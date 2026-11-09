"""
Benchmark serveru FORGE - opakovatelné měření rychlosti a kapacity, se kterým
se porovnává každá optimalizace (baseline -> změna -> znovu benchmark).

    # změřit (~10 min), uloží tests/load/benchmarks/<datum>_<label>_<commit>.json + .md
    .venv/Scripts/python.exe -m tests.load.benchmark run --label baseline
    .venv/Scripts/python.exe -m tests.load.benchmark run --label 4workers --workers 4
    .venv/Scripts/python.exe -m tests.load.benchmark run --label quick --quick      # ~3 min, orientačně

    # porovnat dva běhy (KPI, změna v %, lepší/horší)
    .venv/Scripts/python.exe -m tests.load.benchmark compare tests/load/benchmarks/A.json tests/load/benchmarks/B.json

Server běží v samostatném procesu (uvicorn jako v Dockeru) s vlastními
složkami v tmp. Scénáře jsou pevné - stejné systémy, stejné počty uživatelů -
aby se běhy daly porovnat. Srovnávat se mají běhy ze STEJNÉHO stroje; mezi
stroji je nejstabilnější KPI server_cpu_s (CPU čas na jednu přípravu).
Během měření na stroji nic dalšího náročného nespouštět.

Správnost pod zátěží hlídají pytest testy (tests/load/test_capacity.py);
benchmark jen měří.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import psutil

from tests.load.helpers import (
    BACKEND_ROOT, GOLDEN_PDB_DIR, LiveServer, ServerMonitor, SessionResult, TimedClient, dir_size_mb,
    ff_data_available, latency_stats, per_endpoint, run_concurrently, run_session, start_server,
)

BENCH_DIR = Path(__file__).resolve().parent / "benchmarks"
MEDIUM = "HIS_1CA2_carbonic_anhydrase_zinc.pdb"
SINGLE_SYSTEMS: List[Tuple[str, str, float]] = [
    ("small", "1CRN.pdb", 1.0),          # ~9k atomů
    ("medium", MEDIUM, 1.0),             # ~30k atomů, protein + Zn
    ("rna", "3SKR.pdb", 1.0),            # ~36k atomů
    ("large", "3SKR.pdb", 3.0),          # ~124k atomů
]
MIX = ["1CRN.pdb", "1CRN.pdb", MEDIUM, "3SKR.pdb", "1CRN.pdb", MEDIUM]

# KPI pro porovnání: (klíč, popis, směr) - "lower" = menší je lepší
KPIS: List[Tuple[str, str, str]] = [
    ("single.medium.session_s", "1 user, medium system: session time (s)", "lower"),
    ("single.medium.server_cpu_s", "1 user, medium system: server CPU time (s)", "lower"),
    ("single.small.session_s", "1 user, small system: session time (s)", "lower"),
    ("single.large.session_s", "1 user, large system (124k atoms): session time (s)", "lower"),
    ("single.large.rss_over_idle_mb", "1 user, large system: extra RAM (MB)", "lower"),
    ("concurrency.4.session_p95_s", "4 users at once: slowest wait (s)", "lower"),
    ("concurrency.8.session_p95_s", "8 users at once: slowest wait (s)", "lower"),
    ("concurrency.16.session_p95_s", "16 users at once: slowest wait (s)", "lower"),
    ("concurrency.16.throughput_per_min", "16 users at once: preparations per minute", "higher"),
    ("concurrency.16.rss_peak_mb", "16 users at once: server RAM peak (MB)", "lower"),
    ("responsiveness.health_p95_s", "8 preparing: health check p95 (s)", "lower"),
    ("responsiveness.sequence_p95_s", "8 preparing: sequence request p95 (s)", "lower"),
    ("arrivals.slowdown_p95", "user every 2 s: p95 slowdown vs alone (x)", "lower"),
    ("light_burst.p95_s", "50 users load at once: request p95 (s)", "lower"),
    ("idle.rss_mb", "idle server RAM (MB)", "lower"),
]


def _median(values: List[float]) -> float:
    return round(statistics.median(values), 3)


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=BACKEND_ROOT, capture_output=True, text=True,
                              timeout=10).stdout.strip()
    except Exception:  # noqa: BLE001 - benchmark nesmí spadnout kvůli gitu
        return ""


def _log(text: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {text}", flush=True)


def _ok(results: List[SessionResult], what: str) -> None:
    failed = [r for r in results if not r.ok]
    if failed:
        raise RuntimeError(f"{what}: {len(failed)} session(s) failed, first: {failed[0].error}")


# ---------------------------------------------------------------------------
# Scénáře
# ---------------------------------------------------------------------------

def bench_single(server: LiveServer, repeats: int) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, pdb, padding in SINGLE_SYSTEMS:
        runs = []
        for _ in range(repeats):
            with ServerMonitor(server.pid, server.workspace_dir) as monitor:
                result = run_session(server.base_url, pdb, prepare_overrides={"box_padding_nm": padding})
            _ok([result], f"single {key}")
            runs.append((result, monitor.summary()))
        out[key] = {
            "pdb": pdb, "padding_nm": padding, "atoms": runs[0][0].atoms, "repeats": repeats,
            "session_s": _median([r.seconds for r, _ in runs]),
            "server_cpu_s": _median([s["server_cpu_seconds"] for _, s in runs]),
            "rss_over_idle_mb": max(s["rss_peak_over_baseline_mb"] for _, s in runs),
            "workspace_mb": round(dir_size_mb(server.workspace_dir / runs[0][0].workspace_id), 1),
            # kde uživatel čeká: medián doby každého endpointu přes opakování
            "endpoints": per_endpoint([rec for r, _ in runs for rec in r.records], top=6),
        }
        _log(f"single {key}: {out[key]['session_s']} s, cpu {out[key]['server_cpu_s']} s")
    return out


def bench_concurrency(server: LiveServer, ladder: List[int]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for users in ladder:
        with ServerMonitor(server.pid, server.workspace_dir) as monitor:
            started = time.perf_counter()
            results = run_concurrently(users, lambda _i: run_session(server.base_url, MEDIUM))
            wall = time.perf_counter() - started
        _ok(results, f"concurrency {users}")
        stats = latency_stats([r.seconds for r in results])
        summary = monitor.summary()
        out[str(users)] = {
            "users": users, "wall_s": round(wall, 2),
            "session_p50_s": stats["p50_s"], "session_p95_s": stats["p95_s"], "session_max_s": stats["max_s"],
            "throughput_per_min": round(users / wall * 60, 1),
            "cpu_mean_cores": round(summary["cpu_mean_percent"] / 100, 2),
            "rss_peak_mb": summary["rss_peak_mb"],
        }
        _log(f"concurrency {users}: p95 {stats['p95_s']} s, {out[str(users)]['throughput_per_min']}/min")
    return out


def bench_responsiveness(server: LiveServer, users: int = 8) -> Dict[str, Any]:
    with TimedClient(server.base_url) as client:
        pdb = (GOLDEN_PDB_DIR / "1CRN.pdb").read_bytes()
        ws = client.post("/api/molecules/upload", files={"file": ("1CRN.pdb", pdb)}).json()["workspace_id"]
        done = threading.Event()
        heavy: List[SessionResult] = []

        def background():
            heavy.extend(run_concurrently(users, lambda _i: run_session(server.base_url, MEDIUM)))
            done.set()

        thread = threading.Thread(target=background)
        thread.start()
        health: List[float] = []
        sequence: List[float] = []
        time.sleep(2.0)
        while not done.is_set():
            t0 = time.perf_counter()
            client.get("/api/health")
            health.append(time.perf_counter() - t0)
            t0 = time.perf_counter()
            client.get(f"/api/analysis/sequence/{ws}?fill_gaps=true")
            sequence.append(time.perf_counter() - t0)
            time.sleep(0.25)
        thread.join()
    _ok(heavy, "responsiveness")
    h, s = latency_stats(health), latency_stats(sequence)
    _log(f"responsiveness: health p95 {h['p95_s']} s, sequence p95 {s['p95_s']} s")
    return {"preparing_users": users, "samples": h["n"],
            "health_p50_s": h["p50_s"], "health_p95_s": h["p95_s"], "health_max_s": h["max_s"],
            "sequence_p50_s": s["p50_s"], "sequence_p95_s": s["p95_s"], "sequence_max_s": s["max_s"]}


def bench_arrivals(server: LiveServer, alone: Dict[str, float], interval_s: float = 2.0,
                   users: int = 15) -> Dict[str, Any]:
    def job(i: int) -> SessionResult:
        time.sleep(i * interval_s)
        return run_session(server.base_url, MIX[i % len(MIX)])

    with ServerMonitor(server.pid, server.workspace_dir) as monitor:
        results = run_concurrently(users, job)
    _ok(results, "arrivals")
    slowdown = [r.seconds / alone.get(r.name, r.seconds) for r in results]
    stats = latency_stats([r.seconds for r in results])
    out = {"interval_s": interval_s, "users": users, "rate_per_hour": round(3600 / interval_s),
           "session_p50_s": stats["p50_s"], "session_p95_s": stats["p95_s"],
           "slowdown_p50": round(latency_stats(slowdown)["p50_s"], 2),
           "slowdown_p95": round(latency_stats(slowdown)["p95_s"], 2),
           "cpu_mean_cores": round(monitor.summary()["cpu_mean_percent"] / 100, 2)}
    _log(f"arrivals every {interval_s} s: p95 slowdown {out['slowdown_p95']}x")
    return out


def bench_light_burst(server: LiveServer, users: int = 50) -> Dict[str, Any]:
    pdb = (GOLDEN_PDB_DIR / "1CRN.pdb").read_bytes()

    def browse(_i: int) -> List[float]:
        with TimedClient(server.base_url) as client:
            ws = client.post("/api/molecules/upload", files={"file": ("1CRN.pdb", pdb)}).json()["workspace_id"]
            for path in (f"/api/analysis/sequence/{ws}?fill_gaps=true", f"/api/analysis/altlocs/{ws}",
                         f"/api/analysis/solvation/{ws}", f"/api/download/{ws}?light=true&filename=structure.pdb"):
                client.get(path).raise_for_status()
            return [r.seconds for r in client.records]

    started = time.perf_counter()
    latencies = [x for per_user in run_concurrently(users, browse) for x in per_user]
    wall = time.perf_counter() - started
    stats = latency_stats(latencies)
    _log(f"light burst: p95 {stats['p95_s']} s")
    return {"users": users, "requests": stats["n"], "wall_s": round(wall, 2), "p50_s": stats["p50_s"],
            "p95_s": stats["p95_s"], "max_s": stats["max_s"], "requests_per_s": round(stats["n"] / wall, 1)}


# ---------------------------------------------------------------------------
# Běh, uložení, porovnání
# ---------------------------------------------------------------------------

def run(label: str, workers: Optional[int], quick: bool) -> Path:
    if not ff_data_available():
        sys.exit("Benchmark needs the downloaded FF catalog and data/ff_cache (see tests/README.md).")
    commit = _git("rev-parse", "--short", "HEAD")
    dirty = bool(_git("status", "--porcelain", "--untracked-files=no"))
    result: Dict[str, Any] = {
        "label": label,
        "meta": {
            "date": datetime.now().isoformat(timespec="seconds"),
            "commit": commit + ("-dirty" if dirty else ""),
            "workers": workers or 1,
            "quick": quick,
            "os": platform.platform(terse=True),
            "cpu": f"{psutil.cpu_count(logical=False)} cores / {psutil.cpu_count()} threads",
            "ram_gb": round(psutil.virtual_memory().total / 2 ** 30, 1),
            "python": platform.python_version(),
        },
    }
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="forge_bench_") as tmp, start_server(Path(tmp), workers) as server:
        _log(f"server up ({server.workers} worker(s)), warming up")
        time.sleep(1.0)
        result["idle"] = {"rss_mb": round(ServerMonitor(server.pid, server.workspace_dir)._rss())}
        _ok([run_session(server.base_url, "1CRN.pdb")], "warm-up")   # FF cache, importy, JIT ...
        result["single"] = bench_single(server, repeats=1 if quick else 3)
        result["concurrency"] = bench_concurrency(server, [1, 4, 8] if quick else [1, 2, 4, 8, 16])
        result["responsiveness"] = bench_responsiveness(server)
        alone = {v["pdb"]: v["session_s"] for k, v in result["single"].items() if v["padding_nm"] == 1.0}
        result["arrivals"] = bench_arrivals(server, alone)
        result["light_burst"] = bench_light_burst(server)
    result["meta"]["duration_s"] = round(time.perf_counter() - started)

    BENCH_DIR.mkdir(exist_ok=True)
    stem = f"{datetime.now():%Y-%m-%d_%H%M}_{label}_{commit or 'nogit'}"
    json_path = BENCH_DIR / f"{stem}.json"
    json_path.write_text(json.dumps(result, indent=1), encoding="utf-8")
    (BENCH_DIR / f"{stem}.md").write_text(render(result), encoding="utf-8")
    _log(f"saved {json_path}")
    return json_path


def _get(data: Dict[str, Any], dotted: str) -> Optional[float]:
    node: Any = data
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node if isinstance(node, (int, float)) else None


def _table(rows: List[List[Any]], header: List[str]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def render(result: Dict[str, Any]) -> str:
    meta = result["meta"]
    parts = [f"# FORGE benchmark: {result['label']}", "",
             ", ".join(f"{k} {v}" for k, v in meta.items()), "", "## KPIs", "",
             _table([[desc, _get(result, key)] for key, desc, _ in KPIS], ["KPI", "value"]),
             "## Single user (median of repeats)", "",
             _table([[k, v["atoms"], v["session_s"], v["server_cpu_s"], v["rss_over_idle_mb"], v["workspace_mb"]]
                     for k, v in result["single"].items()],
                    ["system", "atoms", "session s", "server CPU s", "extra RAM MB", "workspace MB"])]
    for key, v in result["single"].items():
        parts += [f"Slowest steps, {key}: " + "; ".join(f"`{e['endpoint']}` p50 {e['p50_s']} s"
                                                        for e in v["endpoints"][:4]), ""]
    parts += ["## Concurrent users (medium system, all start at once)", "",
              _table([[v["users"], v["session_p50_s"], v["session_p95_s"], v["throughput_per_min"],
                       v["cpu_mean_cores"], v["rss_peak_mb"]] for v in result["concurrency"].values()],
                     ["users", "p50 s", "p95 s", "per min", "busy cores", "RAM peak MB"]),
              "## Other scenarios", "", "```json",
              json.dumps({k: result[k] for k in ("responsiveness", "arrivals", "light_burst", "idle")}, indent=1),
              "```", ""]
    return "\n".join(parts)


def compare(a_path: Path, b_path: Path, threshold: float = 0.10) -> str:
    a, b = (json.loads(Path(p).read_text(encoding="utf-8")) for p in (a_path, b_path))
    rows = []
    for key, desc, direction in KPIS:
        va, vb = _get(a, key), _get(b, key)
        if va is None or vb is None:
            rows.append([desc, va, vb, "", "n/a"])
            continue
        change = (vb - va) / va if va else 0.0
        better = change < 0 if direction == "lower" else change > 0
        verdict = "same" if abs(change) < threshold else ("better" if better else "WORSE")
        rows.append([desc, va, vb, f"{change * 100:+.0f} %", verdict])
    head = (f"A: {a['label']} ({a['meta']['date']}, {a['meta']['commit']}, {a['meta']['workers']} worker(s))\n"
            f"B: {b['label']} ({b['meta']['date']}, {b['meta']['commit']}, {b['meta']['workers']} worker(s))\n")
    if a["meta"].get("cpu") != b["meta"].get("cpu"):
        head += "WARNING: different machines - compare only server CPU time.\n"
    return head + "\n" + _table(rows, ["KPI", "A", "B", "change", f"verdict (±{threshold:.0%})"])


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the benchmark")
    r.add_argument("--label", default="run")
    r.add_argument("--workers", type=int, default=None, help="uvicorn --workers (default: 1, like the Dockerfile)")
    r.add_argument("--quick", action="store_true", help="fewer repeats and users (~3 min, rough numbers)")
    c = sub.add_parser("compare", help="compare two benchmark JSON files")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    args = parser.parse_args(argv)
    if args.cmd == "run":
        print(Path(run(args.label, args.workers, args.quick)).with_suffix(".md").read_text(encoding="utf-8"))
    else:
        print(compare(args.a, args.b))


if __name__ == "__main__":
    main()
