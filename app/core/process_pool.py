"""
Pool výpočetních procesů pro těžké bezstavové úlohy (příprava, topologie,
analýza solvatace, PDB z topologie a souřadnic).

Proč procesy: výpočty jsou v Pythonu a drží GIL, takže v jednom procesu
běží vždy jen jeden a event loop na ně čeká (tests/load: jeden proces = ~1
jádro, lehké požadavky mají při zátěži špičky). V poolu běží každá úloha ve
vlastním procesu - víc jader a webový proces zůstane volný.

Proč ne víc procesů uvicornu: interaktivní dostavování bočních řetězců drží
živé objekty builderu v paměti procesu mezi požadavky (sidechain_service.py).
Tady zůstávají ve webovém procesu a do poolu jdou jen bezstavové výpočty.

Nastavení: HEAVY_JOB_WORKERS (app/core/config.py) - počet procesů; 0 = bez
poolu, úlohy běží ve vlákně webového procesu (tak běží testy). Fronta těžkých
úloh (job_queue.py) pouští najednou tolik úloh, kolik je procesů.

Konzole ve frontendu ukazuje zprávy loggeru "app.console" podle workspace
požadavku (app/core/logging.py) - zprávy z procesu poolu se proto zachytí a
ve webovém procesu se zalogují znovu, v kontextu požadavku.
"""

from __future__ import annotations

import asyncio
import functools
import logging
import multiprocessing
import os
import threading
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from typing import Any, Callable, List, Optional, Tuple, TypeVar

from fastapi.concurrency import run_in_threadpool

from app.core.config import settings
from app.core.exceptions import InternalError

logger = logging.getLogger(__name__)

T = TypeVar("T")

_CONSOLE_LOGGER = "app.console"


def configured_workers() -> int:
    """HEAVY_JOB_WORKERS; záporné = automaticky (jádra - 1, nejvýš 4)."""
    workers = settings.HEAVY_JOB_WORKERS
    if workers < 0:
        workers = min(4, max(1, (os.cpu_count() or 2) - 1))
    return workers


# ---------------------------------------------------------------------------
# Strana procesu v poolu
# ---------------------------------------------------------------------------

class _ListHandler(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.INFO)
        self.records: List[Tuple[int, str]] = []

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.records.append((record.levelno, record.getMessage()))
        except Exception:  # noqa: BLE001 - log nesmí shodit úlohu
            pass


def _exit_with_parent() -> None:
    """Když webový proces skončí natvrdo (pád, kill -9, OOM killer), procesy
    poolu by zůstaly viset a držely paměť - ukončí se samy."""
    parent = multiprocessing.parent_process()
    if parent is None:
        return
    parent.join()          # čeká na konec rodiče (sentinel, funguje na Linuxu i Windows)
    os._exit(0)


def _worker_init() -> None:
    from app.core.logging import setup_logging
    setup_logging()
    threading.Thread(target=_exit_with_parent, name="exit-with-parent", daemon=True).start()


def _portable(exc: BaseException) -> BaseException:
    """Výjimka, která přežije cestu do webového procesu (pickle).

    Aplikační chyby (AppBaseException) se, pokud je nejde přenést celé
    (např. ForgeMissingDOFError se staví z celé molekuly), pošlou jako obecná
    AppBaseException se stejným HTTP stavem, kódem, zprávou a payloadem -
    s tím pracuje globální handler, frontend dostane stejnou odpověď.
    """
    import pickle

    from app.core.exceptions import AppBaseException
    try:
        pickle.loads(pickle.dumps(exc))
        return exc
    except Exception:  # noqa: BLE001
        if isinstance(exc, AppBaseException):
            try:
                pickle.dumps(exc.payload)
                payload = exc.payload
            except Exception:  # noqa: BLE001
                payload = {}
            return AppBaseException(exc.message, status_code=exc.status_code, code=exc.code, payload=payload)
        return RuntimeError(f"{type(exc).__name__}: {exc}")


def _run_captured(func: Callable[..., T], args: tuple, kwargs: dict) -> Tuple[str, Any, List[Tuple[int, str]]]:
    """V procesu poolu: spustí úlohu a vrátí (ok|error, výsledek|výjimka, zprávy konzole)."""
    handler = _ListHandler()
    console = logging.getLogger(_CONSOLE_LOGGER)
    console.addHandler(handler)
    try:
        return "ok", func(*args, **kwargs), handler.records
    except BaseException as exc:  # noqa: BLE001 - výjimka se předá webovému procesu
        logging.getLogger(__name__).debug("Job failed in pool worker", exc_info=True)
        return "error", _portable(exc), handler.records
    finally:
        console.removeHandler(handler)


def _warm_up() -> int:
    # Importy builderu a služeb trvají sekundy - zaplatit je při startu, ne u prvního uživatele.
    import app.api.v1.endpoints.analysis  # noqa: F401
    import app.api.v1.endpoints.molecules  # noqa: F401
    import app.api.v1.endpoints.sidechains  # noqa: F401
    import app.api.v1.endpoints.topology  # noqa: F401
    return os.getpid()


# ---------------------------------------------------------------------------
# Strana webového procesu
# ---------------------------------------------------------------------------

class HeavyProcessPool:
    def __init__(self, workers: int):
        self.workers = workers
        self._executor: Optional[ProcessPoolExecutor] = None
        self._lock = threading.Lock()

    @property
    def enabled(self) -> bool:
        return self.workers > 0

    def _get(self) -> ProcessPoolExecutor:
        with self._lock:
            if self._executor is None:
                # spawn: fork procesu s běžícím event loopem a vlákny není bezpečný
                self._executor = ProcessPoolExecutor(
                    max_workers=self.workers, mp_context=multiprocessing.get_context("spawn"),
                    initializer=_worker_init)
            return self._executor

    def _reset(self) -> None:
        with self._lock:
            executor, self._executor = self._executor, None
        if executor is not None:
            executor.shutdown(wait=False, cancel_futures=True)

    async def warm_up(self) -> None:
        if not self.enabled:
            return
        loop = asyncio.get_running_loop()
        executor = self._get()
        pids = await asyncio.gather(*(loop.run_in_executor(executor, _warm_up) for _ in range(self.workers * 2)))
        logger.info(f"Process pool ready: {self.workers} worker(s), pids {sorted(set(pids))}.")

    async def run(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        if not self.enabled:
            return await run_in_threadpool(func, *args, **kwargs)
        loop = asyncio.get_running_loop()
        try:
            status, value, records = await loop.run_in_executor(
                self._get(), functools.partial(_run_captured, func, args, kwargs))
        except BrokenProcessPool:
            # Proces poolu spadl (typicky došla paměť) - pool se založí znovu.
            logger.error(f"Process pool broke while running {getattr(func, '__name__', func)}; restarting it.")
            self._reset()
            raise InternalError("The computation process stopped unexpectedly (possibly out of memory). "
                                "Please try again; very large systems may exceed the server's memory.")
        console = logging.getLogger(_CONSOLE_LOGGER)
        for level, message in records:
            console.log(level, message)
        if status == "error":
            raise value
        return value

    def shutdown(self) -> None:
        self._reset()


heavy_pool = HeavyProcessPool(configured_workers())


async def run_heavy(func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
    """Spustí bezstavovou těžkou úlohu v poolu procesů (nebo ve vlákně, když je pool vypnutý).

    `func` musí být funkce na úrovni modulu a argumenty i výsledek musí jít
    přenést mezi procesy (pickle).
    """
    return await heavy_pool.run(func, *args, **kwargs)
