"""
Fronta těžkých úloh (app/core/job_queue.py): nejvýš `slots` najednou,
ostatní v pořadí příchodu; dekorátor nemění podpis endpointu.
"""

import asyncio
import inspect

import pytest

from app.core.job_queue import HeavyJobQueue, heavy_job

pytestmark = pytest.mark.unit


async def test_jobs_run_one_at_a_time_in_arrival_order():
    queue = HeavyJobQueue(slots=1)
    order, active, peak = [], 0, 0
    gate = asyncio.Event()              # první úloha drží slot, dokud ji nepustíme

    async def job(i):
        nonlocal active, peak
        async with queue.slot("test"):
            active += 1
            peak = max(peak, active)
            order.append(i)
            await gate.wait()
            active -= 1

    tasks = []
    for i in range(6):
        tasks.append(asyncio.create_task(job(i)))
        await asyncio.sleep(0)          # zaručí pořadí příchodu
    await asyncio.sleep(0)
    assert queue.status() == {"slots": 1, "running": 1, "waiting": 5}
    gate.set()
    await asyncio.gather(*tasks)
    assert order == list(range(6)) and peak == 1
    assert queue.status() == {"slots": 1, "running": 0, "waiting": 0}


async def test_earlier_session_goes_first():
    """Druhá úloha uživatele A (příprava -> topologie) předběhne uživatele B,
    který začal až po A - jinak by se session střídaly a všichni skončili až
    na konci."""
    queue = HeavyJobQueue(slots=1)
    order = []
    gate = asyncio.Event()

    async def job(name, ws, wait=False):
        async with queue.slot("test", ws):
            order.append(name)
            if wait:
                await gate.wait()

    a1 = asyncio.create_task(job("A1", "ws-a", wait=True))
    await asyncio.sleep(0)
    b1 = asyncio.create_task(job("B1", "ws-b"))
    await asyncio.sleep(0)
    a2 = asyncio.create_task(job("A2", "ws-a"))
    await asyncio.sleep(0)
    gate.set()
    await asyncio.gather(a1, b1, a2)
    assert order == ["A1", "A2", "B1"]


async def test_cancelled_waiter_does_not_block_the_queue():
    queue = HeavyJobQueue(slots=1)
    gate = asyncio.Event()
    done = []

    async def job(name, wait=False):
        async with queue.slot("test"):
            if wait:
                await gate.wait()
            done.append(name)

    first = asyncio.create_task(job("first", wait=True))
    await asyncio.sleep(0)
    cancelled = asyncio.create_task(job("cancelled"))
    last = asyncio.create_task(job("last"))
    await asyncio.sleep(0)
    cancelled.cancel()
    gate.set()
    await asyncio.gather(first, last)
    assert done == ["first", "last"] and queue.status()["running"] == 0


async def test_more_slots_run_in_parallel():
    queue = HeavyJobQueue(slots=3)
    active = peak = 0

    async def job():
        nonlocal active, peak
        async with queue.slot("test"):
            active += 1
            peak = max(peak, active)
            await asyncio.sleep(0.01)
            active -= 1

    await asyncio.gather(*(job() for _ in range(7)))
    assert peak == 3


async def test_slot_is_released_after_an_error():
    queue = HeavyJobQueue(slots=1)
    with pytest.raises(ValueError):
        async with queue.slot("test"):
            raise ValueError("boom")
    async with queue.slot("test"):      # nezasekne se
        pass


def test_decorator_keeps_the_endpoint_signature():
    async def endpoint(workspace_id: str, hmr: bool = False):
        return workspace_id

    wrapped = heavy_job("test")(endpoint)
    assert inspect.signature(wrapped) == inspect.signature(endpoint)
    assert asyncio.run(wrapped("ws")) == "ws"


def test_queue_endpoint(client):
    assert client.get("/api/queue").json() == {"slots": 1, "running": 0, "waiting": 0}


async def test_status_reports_how_many_jobs_are_ahead():
    queue = HeavyJobQueue(slots=1)
    gate = asyncio.Event()

    async def job(ws):
        async with queue.slot("test", ws):
            await gate.wait()

    tasks = []
    for ws in ("ws-a", "ws-b", "ws-c"):
        tasks.append(asyncio.create_task(job(ws)))
        await asyncio.sleep(0)
    await asyncio.sleep(0)
    assert queue.status("ws-a")["state"] == "running"
    assert queue.status("ws-b") == {"slots": 1, "running": 1, "waiting": 2, "state": "waiting", "ahead": 1}
    assert queue.status("ws-c")["ahead"] == 2
    assert queue.status("ws-other")["state"] is None
    gate.set()
    await asyncio.gather(*tasks)


async def test_identical_request_shares_the_result_but_a_changed_one_runs_again():
    calls = []
    gate = asyncio.Event()

    @heavy_job("dedupe-test")
    async def endpoint(workspace_id: str, setting: int):
        calls.append(setting)
        await gate.wait()
        return {"setting": setting, "call": len(calls)}

    first = asyncio.create_task(endpoint(workspace_id="ws-x", setting=1))
    await asyncio.sleep(0)
    double_click = asyncio.create_task(endpoint(workspace_id="ws-x", setting=1))
    changed = asyncio.create_task(endpoint(workspace_id="ws-x", setting=2))
    await asyncio.sleep(0)
    gate.set()
    a, b, c = await asyncio.gather(first, double_click, changed)
    assert a == b == {"setting": 1, "call": 1}      # dvojklik: spočítáno jednou
    assert c["setting"] == 2 and calls == [1, 2]     # jiné nastavení: nová úloha


async def test_shared_request_gets_the_same_error():
    @heavy_job("dedupe-error")
    async def endpoint(workspace_id: str):
        await asyncio.sleep(0.01)
        raise ValueError("boom")

    results = await asyncio.gather(endpoint(workspace_id="ws-e"), endpoint(workspace_id="ws-e"),
                                   return_exceptions=True)
    assert all(isinstance(r, ValueError) for r in results)


async def test_one_client_cannot_flood_the_queue(monkeypatch):
    from app.core.config import settings
    from app.core.exceptions import TooManyJobsError
    from app.core.request_context import client_address

    monkeypatch.setattr(settings, "HEAVY_JOBS_PER_CLIENT", 2)
    gate = asyncio.Event()

    @heavy_job("limit-test", dedupe=False)
    async def endpoint(workspace_id: str):
        await gate.wait()
        return workspace_id

    token = client_address.set("10.0.0.1")
    try:
        held = [asyncio.create_task(endpoint(workspace_id=f"ws-{i}")) for i in range(2)]
        await asyncio.sleep(0)
        with pytest.raises(TooManyJobsError):                  # třetí od stejné adresy
            await endpoint(workspace_id="ws-3")
    finally:
        client_address.reset(token)
    other = client_address.set("10.0.0.2")                      # jiný klient projde
    try:
        other_task = asyncio.create_task(endpoint(workspace_id="ws-other"))
    finally:
        client_address.reset(other)
    gate.set()
    assert await asyncio.gather(*held, other_task) == ["ws-0", "ws-1", "ws-other"]
