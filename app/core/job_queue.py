"""
Fronta těžkých úloh: příprava, topologie, export ... běží v jednom procesu
nejvýš po HEAVY_JOB_SLOTS najednou, ostatní čekají.

Pořadí: podle toho, kdy workspace (uživatel) poslal svou PRVNÍ těžkou úlohu,
ne podle jednotlivých požadavků. Session se skládá z několika těžkých
požadavků (příprava -> topologie -> export); při řazení po požadavcích čekala
topologie prvního uživatele za přípravami všech ostatních a všichni skončili
zase až na konci (změřeno tests/load: p50 = p95 = ~100 s při 16 uživatelích).
Takhle uživatel, který začal dřív, doběhne celý a pak další.

Proč vůbec fronta: výpočty jsou v Pythonu a drží GIL, souběžná vlákna se jen
střídají a všechna drží data v paměti najednou (16 velkých systémů by se
nevešlo do RAM). Méně vláken navíc méně soupeří o GIL s event loopem, takže
lehké požadavky (sekvence, viewer, health) odpovídají rychleji.

Lehké a interaktivní požadavky (slidery bočních řetězců, sekvence, viewer)
do fronty nepatří - uživatel u nich čeká na okamžitou odezvu.

Použití - dekorátor pod @router.post(...):

    @router.post("/{workspace_id}/generate")
    @heavy_job("topology")
    async def generate_topology(workspace_id: str, ...): ...
"""

from __future__ import annotations

import asyncio
import functools
import hashlib
import heapq
import json
import itertools
import logging
import time
from collections import Counter
from contextlib import asynccontextmanager
from typing import Any, AsyncIterator, Awaitable, Callable, Dict, List, Optional, Tuple, TypeVar

from app.core.config import settings
from app.core.exceptions import TooManyJobsError
from app.core.request_context import client_address

logger = logging.getLogger(__name__)

T = TypeVar("T")

# Jak dlouho si pamatovat, kdy workspace začal (pak se bere jako nový uživatel).
_SESSION_MEMORY_SECONDS = 2 * 60 * 60


class HeavyJobQueue:
    def __init__(self, slots: int):
        self.slots = max(1, slots)
        self.running = 0
        self._heap: List[Tuple[int, int, Optional[str], asyncio.Future]] = []
        self._running_ws: Counter = Counter()
        self._seq = itertools.count()
        self._arrival = itertools.count()
        self._first_seen: Dict[str, Tuple[int, float]] = {}  # workspace -> (pořadí, čas)
        self._loop: Optional[asyncio.AbstractEventLoop] = None

    @property
    def waiting(self) -> int:
        return sum(1 for *_, fut in self._heap if not fut.done())

    def _check_loop(self) -> None:
        # Fronta patří k event loopu; testy (TestClient) můžou loop měnit.
        loop = asyncio.get_running_loop()
        if self._loop is not loop:
            self._loop, self.running, self._heap, self._running_ws = loop, 0, [], Counter()

    def _priority(self, workspace_id: Optional[str]) -> int:
        """Pořadí příchodu session (čítač - time.monotonic má na Windows
        rozlišení ~15 ms a dva uživatelé by dostali stejnou hodnotu)."""
        order = next(self._arrival)
        if not workspace_id:
            return order
        now = time.monotonic()
        if len(self._first_seen) > 1000:
            self._first_seen = {k: v for k, v in self._first_seen.items() if now - v[1] < _SESSION_MEMORY_SECONDS}
        first = self._first_seen.get(workspace_id)
        if first is None or now - first[1] > _SESSION_MEMORY_SECONDS:
            self._first_seen[workspace_id] = (order, now)
            return order
        return first[0]

    def _wake_next(self) -> None:
        while self.running < self.slots and self._heap:
            *_, fut = heapq.heappop(self._heap)
            if not fut.done():          # zrušené čekání (klient odešel) přeskočit
                self.running += 1
                fut.set_result(None)

    @asynccontextmanager
    async def slot(self, kind: str, workspace_id: Optional[str] = None) -> AsyncIterator[None]:
        self._check_loop()
        priority = self._priority(workspace_id)
        queued_at = time.perf_counter()
        if self.running < self.slots and not self._heap:
            self.running += 1
        else:
            fut = asyncio.get_running_loop().create_future()
            heapq.heappush(self._heap, (priority, next(self._seq), workspace_id, fut))
            try:
                await fut
            except BaseException:
                if fut.done() and not fut.cancelled():
                    # slot už byl přidělen, ale čekající mezitím zrušen - vrátit ho
                    self.running -= 1
                    self._wake_next()
                raise
        waited = time.perf_counter() - queued_at
        if waited > 1.0:
            logger.info(f"Job queue: {kind} waited {waited:.1f} s ({self.waiting} still waiting).")
        self._running_ws[workspace_id] += 1
        try:
            yield
        finally:
            self._running_ws[workspace_id] -= 1
            if self._running_ws[workspace_id] <= 0:
                del self._running_ws[workspace_id]
            self.running -= 1
            self._wake_next()

    def status(self, workspace_id: Optional[str] = None) -> Dict[str, Any]:
        """Celkový stav; s workspace_id i stav jeho úlohy a kolik úloh je před ní
        (pro hlášku "server je vytížený, před vámi jsou N" ve frontendu)."""
        out: Dict[str, Any] = {"slots": self.slots, "running": self.running, "waiting": self.waiting}
        if workspace_id is None:
            return out
        waiting = sorted((p, seq, ws) for p, seq, ws, fut in self._heap if not fut.done())
        mine = next(((p, seq) for p, seq, ws in waiting if ws == workspace_id), None)
        if mine is not None:
            out.update(state="waiting", ahead=self.running + sum(1 for p, seq, _ in waiting if (p, seq) < mine))
        else:
            out.update(state="running" if self._running_ws.get(workspace_id) else None, ahead=0)
        return out


def _configured_slots() -> int:
    if settings.HEAVY_JOB_SLOTS > 0:
        return settings.HEAVY_JOB_SLOTS
    from app.core.process_pool import configured_workers
    return max(1, configured_workers())


heavy_jobs = HeavyJobQueue(_configured_slots())


def _workspace_of(kwargs: Dict[str, Any]) -> Optional[str]:
    if isinstance(kwargs.get("workspace_id"), str):
        return kwargs["workspace_id"]
    request = kwargs.get("request")
    return getattr(request, "workspace_id", None)


# Kolik těžkých úloh (běžících + čekajících) má který klient - ochrana proti botům.
_client_jobs: Counter = Counter()


def _check_client_limit(client: Optional[str]) -> None:
    limit = settings.HEAVY_JOBS_PER_CLIENT
    if limit > 0 and client and _client_jobs[client] >= limit:
        logger.warning(f"Job queue: client {client} already has {_client_jobs[client]} heavy jobs - rejected.")
        raise TooManyJobsError(
            f"Too many calculations from your address at once ({_client_jobs[client]}). "
            "Please wait until some of them finish and try again.")


# Rozpracované těžké požadavky: (druh, workspace, otisk požadavku) -> výsledek.
_inflight: Dict[Tuple[str, str, str], asyncio.Future] = {}


def _fingerprint(kwargs: Dict[str, Any]) -> str:
    """Otisk parametrů požadavku - stejné nastavení = stejný otisk."""
    def plain(value: Any) -> Any:
        dump = getattr(value, "model_dump", None)
        return dump(mode="json") if callable(dump) else value
    data = json.dumps({k: plain(v) for k, v in kwargs.items()}, sort_keys=True, default=str)
    return hashlib.sha1(data.encode("utf-8")).hexdigest()


def heavy_job(kind: str, dedupe: bool = True) -> Callable[[Callable[..., Awaitable[T]]], Callable[..., Awaitable[T]]]:
    """Async endpoint poběží až po získání slotu ve frontě (podpis zůstává, FastAPI ho vidí).

    dedupe: když stejný workspace pošle PŘESNĚ stejný požadavek (dvojklik,
    znovu odeslání po obnovení stránky) a předchozí ještě běží nebo čeká,
    nepočítá se znovu - oba dostanou stejný výsledek. Jiné nastavení = jiný
    otisk = normální nová úloha (uživatel to chce jinak). Vypnout u
    endpointů, které vrací stream (export), ten nejde poslat dvakrát.
    """
    def decorator(func: Callable[..., Awaitable[T]]) -> Callable[..., Awaitable[T]]:
        @functools.wraps(func)
        async def wrapper(*args: Any, **kwargs: Any) -> T:
            workspace_id = _workspace_of(kwargs)
            key = (kind, workspace_id, _fingerprint(kwargs)) if dedupe and workspace_id else None
            running = _inflight.get(key) if key else None
            if running is not None and running.get_loop() is asyncio.get_running_loop():
                logger.info(f"Job queue: identical {kind} for {workspace_id} already in progress - sharing its result.")
                return await asyncio.shield(running)
            client = client_address.get()
            _check_client_limit(client)
            shared = asyncio.get_running_loop().create_future() if key else None
            if key:
                _inflight[key] = shared
            _client_jobs[client] += 1
            try:
                async with heavy_jobs.slot(kind, workspace_id):
                    result = await func(*args, **kwargs)
                if shared is not None:
                    shared.set_result(result)
                return result
            except BaseException as exc:
                if shared is not None and not shared.done():
                    if isinstance(exc, asyncio.CancelledError):
                        shared.cancel()
                    else:
                        shared.set_exception(exc)
                        shared.exception()     # označit jako převzatou, ať asyncio nevaruje
                raise
            finally:
                _client_jobs[client] -= 1
                if _client_jobs[client] <= 0:
                    del _client_jobs[client]
                if key and _inflight.get(key) is shared:
                    del _inflight[key]
        return wrapper
    return decorator
