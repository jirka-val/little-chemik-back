"""
Pomůcky pro zátěžové testy (tests/load/): klient s měřením odezvy každého
requestu, sledování paměti/CPU serverového procesu a velikosti workspaces,
jedna "uživatelská session" (celá příprava jako ve frontendu) a statistiky.
"""

from __future__ import annotations

import re
import statistics
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional

import httpx
import psutil

from tests.golden.pipeline import PipelineCase, run_pipeline

GOLDEN_PDB_DIR = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "golden"
_UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
_QUERY = re.compile(r"\?.*$")

# Delší než nejpomalejší krok při plném souběhu (příprava 125k atomů za
# frontou dalších uživatelů) - timeout tu nemá nic měřit, jen nezaseknout běh.
REQUEST_TIMEOUT = 900.0


def endpoint(method: str, url: str) -> str:
    """"POST /api/sidechains/start/{ws}" - workspace id a query pryč, ať se dá agregovat."""
    path = _QUERY.sub("", httpx.URL(url).path)
    return f"{method} {_UUID.sub('{ws}', path)}"


# ---------------------------------------------------------------------------
# Server v samostatném procesu
# ---------------------------------------------------------------------------

BACKEND_ROOT = Path(__file__).resolve().parents[2]


class ServerStartError(RuntimeError):
    pass


@dataclass
class LiveServer:
    base_url: str
    pid: int
    workspace_dir: Path
    incidents_dir: Path
    log_path: Path
    workers: int


def ff_data_available() -> bool:
    ff_cache = BACKEND_ROOT / "data" / "ff_cache"
    return ((BACKEND_ROOT / "data" / "ff_catalog.json").exists() and ff_cache.exists()
            and any(ff_cache.iterdir()))


def _free_port() -> int:
    import socket
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@contextmanager
def start_server(root: Path, workers: Optional[int] = None, pool_workers: Optional[int] = None) -> Iterator[LiveServer]:
    """uvicorn app.main:app v samostatném procesu (jako Dockerfile: bez reloadu),
    s vlastními WORKSPACE_DIR a INCIDENTS_DIR pod `root`. pool_workers =
    HEAVY_JOB_WORKERS (None = výchozí nastavení serveru, 0 = bez poolu)."""
    import os
    import subprocess
    import sys

    workspace_dir, incidents_dir = root / "workspaces", root / "incidents"
    workspace_dir.mkdir(parents=True, exist_ok=True)
    incidents_dir.mkdir(parents=True, exist_ok=True)
    port = _free_port()
    env = {**os.environ, "WORKSPACE_DIR": str(workspace_dir), "INCIDENTS_DIR": str(incidents_dir),
           "PYTHONUNBUFFERED": "1"}
    # tests/conftest.py vypíná pool výpočetních procesů pro TestClient - server
    # tady má běžet jako v produkci (výchozí nastavení), pokud ho volající neurčí.
    env.pop("HEAVY_JOB_WORKERS", None)
    # Všichni simulovaní uživatelé chodí z 127.0.0.1 - limit na klienta (ochrana
    # proti botům) by měření kapacity zkreslil.
    env["HEAVY_JOBS_PER_CLIENT"] = "0"
    if pool_workers is not None:
        env["HEAVY_JOB_WORKERS"] = str(pool_workers)
    cmd = [sys.executable, "-m", "uvicorn", "app.main:app", "--host", "127.0.0.1", "--port", str(port),
           "--log-level", "warning"] + (["--workers", str(workers)] if workers else [])
    log_path = root / "server.log"
    log = open(log_path, "w", encoding="utf-8")
    process = subprocess.Popen(cmd, cwd=BACKEND_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    base_url = f"http://127.0.0.1:{port}"

    # tests/golden/pipeline.py čte výstupy (structure.pdb, prmtop ...) přímo z
    # disku přes workspace_manager - v tomto procesu musí mířit na složku
    # serveru (běží na stejném stroji).
    from app.workspaces import manager
    original_dir = manager.WORKSPACE_DIR
    try:
        deadline = time.time() + 90
        while True:
            try:
                if httpx.get(f"{base_url}/api/health", timeout=2).status_code == 200:
                    break
            except httpx.HTTPError:
                pass
            if process.poll() is not None or time.time() > deadline:
                log.flush()
                raise ServerStartError(f"Server did not start:\n{log_path.read_text(encoding='utf-8')[-3000:]}")
            time.sleep(0.3)
        manager.WORKSPACE_DIR = str(workspace_dir)
        yield LiveServer(base_url, process.pid, workspace_dir, incidents_dir, log_path, workers or 1)
    finally:
        manager.WORKSPACE_DIR = original_dir
        # Celý strom procesů: na Windows terminate() ukončí jen server a procesy
        # poolu (app/core/process_pool.py) by přežily a držely otevřený log.
        try:
            children = psutil.Process(process.pid).children(recursive=True)
        except psutil.Error:
            children = []
        process.terminate()
        try:
            process.wait(timeout=20)
        except subprocess.TimeoutExpired:
            process.kill()
        for child in children:
            try:
                child.kill()
            except psutil.Error:
                pass
        psutil.wait_procs(children, timeout=10)
        log.close()


# ---------------------------------------------------------------------------
# Klient s měřením
# ---------------------------------------------------------------------------

@dataclass
class RequestRecord:
    endpoint: str
    status: int
    seconds: float
    started: float


class TimedClient(httpx.Client):
    """httpx.Client, který si pamatuje dobu každého requestu (včetně čtení těla).

    Má stejné rozhraní jako FastAPI TestClient, takže ho bere i
    tests/golden/pipeline.run_pipeline - zátěž jde přesně cestou frontendu.
    """

    def __init__(self, base_url: str, records: Optional[List[RequestRecord]] = None):
        super().__init__(base_url=base_url, timeout=REQUEST_TIMEOUT)
        self.records: List[RequestRecord] = records if records is not None else []
        self._lock = threading.Lock()

    def request(self, method, url, *args, **kwargs):  # type: ignore[override]
        started = time.perf_counter()
        response = super().request(method, url, *args, **kwargs)
        response.read()
        elapsed = time.perf_counter() - started
        with self._lock:
            self.records.append(RequestRecord(endpoint(method, str(response.request.url)),
                                              response.status_code, elapsed, started))
        return response


# ---------------------------------------------------------------------------
# Sledování serveru
# ---------------------------------------------------------------------------

def dir_size_mb(path: Path) -> float:
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:  # soubor mezitím smazal úklid / jiný request
            pass
    return total / 2 ** 20


class ServerMonitor:
    """Vzorkuje RSS a CPU serverového procesu (+ potomků) a velikost workspaces."""

    def __init__(self, pid: int, workspace_dir: Path, interval: float = 0.1, disk_every: int = 10):
        self.process = psutil.Process(pid)
        self.workspace_dir = workspace_dir
        self.interval = interval
        self.disk_every = disk_every
        self.rss_mb: List[float] = []
        self.cpu_percent: List[float] = []
        self.disk_mb: List[float] = []
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _processes(self):
        try:
            return [self.process] + self.process.children(recursive=True)
        except psutil.Error:
            return [self.process]

    def _rss(self) -> float:
        total = 0
        for proc in self._processes():
            try:
                total += proc.memory_info().rss
            except psutil.Error:
                pass
        return total / 2 ** 20

    def _cpu_seconds(self) -> float:
        """Součet CPU času procesu a potomků - na Windows je python.exe z venv
        jen spouštěč a server běží v potomkovi."""
        total = 0.0
        for proc in self._processes():
            try:
                times = proc.cpu_times()
                total += times.user + times.system
            except psutil.Error:
                pass
        return total

    def _run(self):
        last_cpu, last_t = self._cpu_seconds(), time.perf_counter()
        n = 0
        while not self._stop.is_set():
            self.rss_mb.append(self._rss())
            cpu, now = self._cpu_seconds(), time.perf_counter()
            if now > last_t:
                self.cpu_percent.append(100.0 * (cpu - last_cpu) / (now - last_t))
            last_cpu, last_t = cpu, now
            if n % self.disk_every == 0:
                self.disk_mb.append(dir_size_mb(self.workspace_dir))
            n += 1
            self._stop.wait(self.interval)

    def __enter__(self) -> "ServerMonitor":
        self.baseline_rss_mb = self._rss()
        self.baseline_disk_mb = dir_size_mb(self.workspace_dir)
        self.cpu_start = self._cpu_seconds()
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join()
        self.cpu_seconds = self._cpu_seconds() - self.cpu_start
        self.final_disk_mb = dir_size_mb(self.workspace_dir)

    def summary(self) -> Dict[str, float]:
        rss = self.rss_mb or [self.baseline_rss_mb]
        cpu = self.cpu_percent or [0.0]
        return {
            "rss_baseline_mb": round(self.baseline_rss_mb),
            "rss_peak_mb": round(max(rss)),
            "rss_peak_over_baseline_mb": round(max(rss) - self.baseline_rss_mb),
            # 100 % = jedno plně vytížené jádro
            "cpu_mean_percent": round(statistics.fmean(cpu)),
            "cpu_peak_percent": round(max(cpu)),
            "server_cpu_seconds": round(self.cpu_seconds, 1),
            "disk_added_mb": round(self.final_disk_mb - self.baseline_disk_mb, 1),
        }


# ---------------------------------------------------------------------------
# Uživatelská session
# ---------------------------------------------------------------------------

@dataclass
class SessionResult:
    name: str
    ok: bool
    seconds: float
    workspace_id: Optional[str] = None
    error: Optional[str] = None
    atoms: Optional[int] = None
    prmtop_sha: Optional[str] = None
    records: List[RequestRecord] = field(default_factory=list)


def run_session(base_url: str, pdb_name: str, **case_kwargs: Any) -> SessionResult:
    """Jeden uživatel: upload -> sekvence/altloc -> FF -> příprava -> topologie
    -> viewer -> mdin -> export ZIP (tests/golden/pipeline.run_pipeline)."""
    import hashlib

    records: List[RequestRecord] = []
    started = time.perf_counter()
    with TimedClient(base_url, records) as client:
        try:
            path = GOLDEN_PDB_DIR / pdb_name
            if not path.exists():           # další fixtures (např. glu_open_branch.pdb)
                path = GOLDEN_PDB_DIR.parent / pdb_name
            result = run_pipeline(client, PipelineCase(name=pdb_name, pdb_path=path, **case_kwargs))
        except BaseException as exc:  # noqa: BLE001 - chyba jednoho uživatele je výsledek, ne pád běhu
            return SessionResult(pdb_name, False, time.perf_counter() - started, error=f"{type(exc).__name__}: {exc}"[:500],
                                 records=records)
    crd = result.files["structure.crd"].splitlines()
    return SessionResult(
        pdb_name, True, time.perf_counter() - started, workspace_id=result.workspace_id,  # type: ignore[attr-defined]
        atoms=int(crd[1].split()[0]),
        # bez 1. řádku: %VERSION ... DATE = čas zápisu, jinak je výstup deterministický
        prmtop_sha=hashlib.sha256(result.files["structure.prmtop"].split("\n", 1)[1].encode()).hexdigest()[:16],
        records=records,
    )


def run_concurrently(n: int, job: Callable[[int], Any]) -> List[Any]:
    """n uživatelů najednou - všichni startují ve stejnou chvíli (bariéra)."""
    barrier = threading.Barrier(n)

    def start(i: int):
        barrier.wait()
        return job(i)

    with ThreadPoolExecutor(max_workers=n) as pool:
        return list(pool.map(start, range(n)))


# ---------------------------------------------------------------------------
# Statistiky
# ---------------------------------------------------------------------------

def percentile(values: List[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    k = (len(ordered) - 1) * q / 100
    lo, hi = int(k), min(int(k) + 1, len(ordered) - 1)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (k - lo)


def latency_stats(values: List[float]) -> Dict[str, float]:
    return {
        "n": len(values),
        "p50_s": round(percentile(values, 50), 3),
        "p95_s": round(percentile(values, 95), 3),
        "max_s": round(max(values), 3) if values else float("nan"),
    }


def per_endpoint(records: List[RequestRecord], top: int = 8) -> List[Dict[str, Any]]:
    """Nejpomalejší endpointy podle součtu času (kde uživatel čeká nejvíc)."""
    groups: Dict[str, List[float]] = {}
    for rec in records:
        groups.setdefault(rec.endpoint, []).append(rec.seconds)
    rows = [{"endpoint": ep, **latency_stats(v), "total_s": round(sum(v), 1)} for ep, v in groups.items()]
    return sorted(rows, key=lambda r: -r["total_s"])[:top]
