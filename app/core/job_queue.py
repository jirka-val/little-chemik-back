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
import heapq
import itertools
import logging
import time
from contextlib import asynccontextmanager
from typing import Any, AsyncIterator, Awaitable, Callable, Dict, List, Optional, Tuple, TypeVar

from app.core.config import settings

logger = logging.getLogger(__name__)

T = TypeVar("T")

# Jak dlouho si pamatovat, kdy workspace začal (pak se bere jako nový uživatel).
_SESSION_MEMORY_SECONDS = 2 * 60 * 60


class HeavyJobQueue:
    def __init__(self, slots: int):
        self.slots = max(1, slots)
        self.running = 0
        self._heap: List[Tuple[int, int, asyncio.Future]] = []
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
            self._loop, self.running, self._heap = loop, 0, []

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
            heapq.heappush(self._heap, (priority, next(self._seq), fut))
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
        try:
            yield
        finally:
            self.running -= 1
            self._wake_next()

    def status(self) -> Dict[str, Any]:
        return {"slots": self.slots, "running": self.running, "waiting": self.waiting}


heavy_jobs = HeavyJobQueue(settings.HEAVY_JOB_SLOTS)


def _workspace_of(kwargs: Dict[str, Any]) -> Optional[str]:
    if isinstance(kwargs.get("workspace_id"), str):
        return kwargs["workspace_id"]
    request = kwargs.get("request")
    return getattr(request, "workspace_id", None)


def heavy_job(kind: str) -> Callable[[Callable[..., Awaitable[T]]], Callable[..., Awaitable[T]]]:
    """Async endpoint poběží až po získání slotu ve frontě (podpis zůstává, FastAPI ho vidí)."""
    def decorator(func: Callable[..., Awaitable[T]]) -> Callable[..., Awaitable[T]]:
        @functools.wraps(func)
        async def wrapper(*args: Any, **kwargs: Any) -> T:
            async with heavy_jobs.slot(kind, _workspace_of(kwargs)):
                return await func(*args, **kwargs)
        return wrapper
    return decorator
