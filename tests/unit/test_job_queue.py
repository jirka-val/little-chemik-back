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
