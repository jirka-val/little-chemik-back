"""
Kapacita serveru FORGE: kolik uživatelů najednou, kolik paměti, CPU a disku.

Každý test běží proti skutečnému serveru (fixture live_server) a zapisuje
naměřená čísla do reportu (tests/load/results/latest.md). Asserty hlídají
správnost (žádné chyby, uživatelé si nelezou do dat, server odpovídá);
kapacitní čísla se jen zapisují - závisí na stroji.

    pytest -m load tests/load -s
    LOAD_WORKERS=4 pytest -m load tests/load -s     # uvicorn --workers 4
"""

from __future__ import annotations

import threading
import time
from typing import Dict, List

import httpx
import pytest

from tests.load.helpers import (
    GOLDEN_PDB_DIR, ServerMonitor, SessionResult, TimedClient, dir_size_mb, latency_stats, per_endpoint,
    run_concurrently, run_session,
)

pytestmark = pytest.mark.load

MEDIUM = "HIS_1CA2_carbonic_anhydrase_zinc.pdb"   # protein 260 reziduí + Zn, ~30k atomů po solvataci
SYSTEMS = [
    ("1CRN.pdb", 1.0, "small protein"),
    (MEDIUM, 1.0, "medium protein + Zn"),
    ("3SKR.pdb", 1.0, "RNA"),
    ("3SKR.pdb", 3.0, "RNA, 3 nm water"),
    ("3SKR.pdb", 5.0, "RNA, 5 nm water"),
]

_baseline: Dict[str, float] = {}    # doba session jednoho uživatele (pro zpomalení při souběhu)
_sessions: List[SessionResult] = []  # všechny úspěšné session (disk na session)


def _errors(results: List[SessionResult]) -> List[str]:
    return [f"{r.name}: {r.error}" for r in results if not r.ok]


# ---------------------------------------------------------------------------
# 1. Jeden uživatel - náročnost jedné přípravy podle velikosti systému
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pdb, padding, label", SYSTEMS, ids=[f"{s[0]}-{s[1]}nm" for s in SYSTEMS])
def test_single_user_baseline(live_server, load_report, pdb, padding, label):
    with ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        result = run_session(live_server.base_url, pdb, prepare_overrides={"box_padding_nm": padding})
    assert result.ok, result.error
    _sessions.append(result)
    _baseline.setdefault(f"{pdb}@{padding}", result.seconds)
    slow = per_endpoint(result.records, top=3)
    stats = monitor.summary()
    load_report.add(f"Single user: {pdb} ({label}, {padding} nm water)", [{
        "atoms": result.atoms,
        "session_s": round(result.seconds, 1),
        "server_cpu_s": stats["server_cpu_seconds"],
        "rss_peak_mb": stats["rss_peak_mb"],
        "rss_over_idle_mb": stats["rss_peak_over_baseline_mb"],
        "workspace_mb": round(dir_size_mb(live_server.workspace_dir / result.workspace_id), 1),
        "requests": len(result.records),
    }], note="Slowest steps: " + "; ".join(f"`{r['endpoint']}` {r['total_s']} s" for r in slow))


# ---------------------------------------------------------------------------
# 2. Souběh - N uživatelů klikne "Prepare" ve stejnou chvíli
# ---------------------------------------------------------------------------

_concurrency_rows: List[Dict] = []


@pytest.mark.parametrize("users", [1, 2, 4, 8, 16])
def test_concurrent_users(live_server, load_report, users):
    with ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        started = time.perf_counter()
        results = run_concurrently(users, lambda _i: run_session(live_server.base_url, MEDIUM))
        wall = time.perf_counter() - started
    assert not _errors(results), _errors(results)
    # Stejný vstup -> bajtově stejná topologie u všech (kromě data v %VERSION):
    # session se navzájem neovlivňují.
    assert len({r.atoms for r in results}) == 1
    assert len({r.prmtop_sha for r in results}) == 1
    _sessions.extend(results)

    single = _baseline.get(f"{MEDIUM}@1.0") or min(r.seconds for r in results)
    times = [r.seconds for r in results]
    stats = monitor.summary()
    _concurrency_rows.append({
        "users": users,
        "wall_s": round(wall, 1),
        **{f"session_{k}": v for k, v in latency_stats(times).items() if k != "n"},
        "slowdown_p95": f"{latency_stats(times)['p95_s'] / single:.1f}x",
        "throughput_per_min": round(users / wall * 60, 1),
        "cpu_mean_cores": round(stats["cpu_mean_percent"] / 100, 2),
        "rss_peak_mb": stats["rss_peak_mb"],
        "identical_prmtop": len({r.prmtop_sha for r in results}) == 1,
    })
    if users == 16:
        load_report.add(
            f"Concurrent users ({MEDIUM}, all start at once)", _concurrency_rows,
            note=f"One user alone needs {single:.1f} s. `slowdown_p95` = how much longer the slowest users wait. "
                 "`cpu_mean_cores` = average number of busy cores of the server process.")


# ---------------------------------------------------------------------------
# 3. Odezva lehkých požadavků, zatímco jiní připravují
# ---------------------------------------------------------------------------

def test_light_requests_stay_responsive_under_load(live_server, load_report):
    with TimedClient(live_server.base_url) as client:
        pdb = (GOLDEN_PDB_DIR / "1CRN.pdb").read_bytes()
        ws = client.post("/api/molecules/upload", files={"file": ("1CRN.pdb", pdb)}).json()["workspace_id"]

        done = threading.Event()
        heavy: List[SessionResult] = []

        def background():
            heavy.extend(run_concurrently(8, lambda _i: run_session(live_server.base_url, MEDIUM)))
            done.set()

        thread = threading.Thread(target=background)
        thread.start()
        health: List[float] = []
        sequence: List[float] = []
        time.sleep(2.0)  # ať už všech 8 příprav běží
        while not done.is_set():
            t0 = time.perf_counter()
            client.get("/api/health")
            health.append(time.perf_counter() - t0)
            t0 = time.perf_counter()
            assert client.get(f"/api/analysis/sequence/{ws}?fill_gaps=true").status_code == 200
            sequence.append(time.perf_counter() - t0)
            time.sleep(0.25)
        thread.join()
    assert not _errors(heavy), _errors(heavy)
    _sessions.extend(heavy)
    load_report.add("Light requests while 8 users prepare", [
        {"request": "GET /api/health", **latency_stats(health)},
        {"request": "GET /api/analysis/sequence/{ws}", **latency_stats(sequence)},
    ], note="If the health check is slow too, heavy work blocks the event loop and the whole site freezes.")
    # Health je čistě async - musí odpovědět, i když server počítá.
    assert latency_stats(health)["p95_s"] < 1.0


# ---------------------------------------------------------------------------
# 4. Realistický provoz - uživatelé přicházejí průběžně, různé systémy
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("interval_s, users", [(6.0, 15), (2.0, 15)], ids=["every-6s", "every-2s"])
def test_users_arriving_over_time(live_server, load_report, interval_s, users):
    mix = ["1CRN.pdb", "1CRN.pdb", MEDIUM, "3SKR.pdb", "1CRN.pdb", MEDIUM]

    def job(i: int) -> SessionResult:
        time.sleep(i * interval_s)
        return run_session(live_server.base_url, mix[i % len(mix)])

    with ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        started = time.perf_counter()
        results = run_concurrently(users, job)
        wall = time.perf_counter() - started
    assert not _errors(results), _errors(results)
    _sessions.extend(results)
    slowdown = [r.seconds / _baseline.get(f"{r.name}@1.0", r.seconds) for r in results]
    stats = monitor.summary()
    load_report.add(f"Users arriving every {interval_s:g} s ({users} users, mixed systems)", [{
        "rate_per_hour": round(3600 / interval_s),
        "wall_s": round(wall, 1),
        **{f"session_{k}": v for k, v in latency_stats([r.seconds for r in results]).items() if k != "n"},
        "slowdown_p95": f"{latency_stats(slowdown)['p95_s']:.1f}x",
        "cpu_mean_cores": round(stats["cpu_mean_percent"] / 100, 2),
        "rss_peak_mb": stats["rss_peak_mb"],
    }], note="Mix: 50 % small protein, 33 % medium protein, 17 % RNA. "
             "`slowdown_p95` compares each session with the same system run alone.")


# ---------------------------------------------------------------------------
# 5. Mnoho lehkých návštěv najednou (jen nahrání a prohlížení)
# ---------------------------------------------------------------------------

def test_burst_of_light_sessions(live_server, load_report):
    pdb = (GOLDEN_PDB_DIR / "1CRN.pdb").read_bytes()

    def browse(_i: int) -> List[float]:
        with TimedClient(live_server.base_url) as client:
            ws = client.post("/api/molecules/upload", files={"file": ("1CRN.pdb", pdb)}).json()["workspace_id"]
            for path in (f"/api/analysis/sequence/{ws}?fill_gaps=true", f"/api/analysis/altlocs/{ws}",
                         f"/api/analysis/solvation/{ws}", f"/api/download/{ws}?light=true&filename=structure.pdb"):
                assert client.get(path).status_code == 200, path
            return [r.seconds for r in client.records]

    with ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        started = time.perf_counter()
        latencies = [x for per_user in run_concurrently(50, browse) for x in per_user]
        wall = time.perf_counter() - started
    stats = monitor.summary()
    load_report.add("50 users load a structure at the same moment (upload + sequence + altlocs + viewer)", [{
        "requests": len(latencies),
        "wall_s": round(wall, 1),
        **{k: v for k, v in latency_stats(latencies).items() if k != "n"},
        "requests_per_s": round(len(latencies) / wall, 1),
        "rss_peak_mb": stats["rss_peak_mb"],
    }])


# ---------------------------------------------------------------------------
# 5b. Interaktivní dostavování bočních řetězců (start -> commit) při souběhu
# ---------------------------------------------------------------------------

def test_interactive_sidechain_sessions_under_load(live_server, load_report):
    """Relace drží objekty builderu v paměti webového procesu - výpočet v poolu
    (app/core/process_pool.py) ji vrátí jako needs_session a /start ji
    založí znovu ve webovém procesu. commit pak musí trefit tu samou relaci."""
    results = run_concurrently(6, lambda _i: run_session(live_server.base_url, "glu_open_branch.pdb"))
    assert not _errors(results), _errors(results)
    for result in results:
        commits = [r for r in result.records if r.endpoint == "POST /api/sidechains/commit/{ws}"]
        assert len(commits) == 1 and commits[0].status == 200
    assert len({r.prmtop_sha for r in results}) == 1
    load_report.add("Interactive side-chain sessions (6 users at once, start -> commit -> topology)", [{
        **{f"session_{k}": v for k, v in latency_stats([r.seconds for r in results]).items()},
    }])


# ---------------------------------------------------------------------------
# 6. Velký systém: PDB bez CRYST1 (rekonstrukce boxu) a topologie + souřadnice
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def large_system(live_server):
    result = run_session(live_server.base_url, "3SKR.pdb", prepare_overrides={"box_padding_nm": 3.0})
    assert result.ok, result.error
    ws_dir = live_server.workspace_dir / result.workspace_id
    return {
        "atoms": result.atoms,
        "pdb": (ws_dir / "structure.pdb").read_text(encoding="utf-8"),
        "prmtop": (ws_dir / "structure.prmtop").read_bytes(),
        "crd": (ws_dir / "structure.crd").read_bytes(),
    }


def test_large_pdb_without_cryst1(live_server, load_report, large_system):
    pdb = "".join(line + "\n" for line in large_system["pdb"].splitlines() if not line.startswith("CRYST1"))
    with TimedClient(live_server.base_url) as client, \
            ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        ws = client.post("/api/molecules/upload", files={"file": ("big.pdb", pdb.encode())}).json()["workspace_id"]
        t0 = time.perf_counter()
        analysis = client.get(f"/api/analysis/solvation/{ws}").json()
        seconds = time.perf_counter() - t0
    assert analysis["box_source"] == "reconstructed", analysis
    stats = monitor.summary()
    load_report.add("Large solvated PDB without CRYST1 (box reconstructed from the water)", [{
        "atoms": large_system["atoms"],
        "analysis_s": round(seconds, 1),
        "shape": analysis["box"]["shape"],
        "rss_over_idle_mb": stats["rss_peak_over_baseline_mb"],
    }])


def test_large_topology_and_coordinates_upload(live_server, load_report, large_system):
    with TimedClient(live_server.base_url) as client, \
            ServerMonitor(live_server.pid, live_server.workspace_dir) as monitor:
        t0 = time.perf_counter()
        response = client.post("/api/molecules/upload-inputs", files={
            "topology": ("big.prmtop", large_system["prmtop"]), "coordinates": ("big.rst7", large_system["crd"])})
        seconds = time.perf_counter() - t0
        assert response.status_code == 200, response.text
        ws = response.json()["workspace_id"]
        analysis = client.get(f"/api/analysis/solvation/{ws}").json()
    assert response.json()["built"]["atoms"] == large_system["atoms"]
    assert analysis["box_source"] == "cryst1"
    stats = monitor.summary()
    load_report.add("Large system from topology + coordinates (no PDB)", [{
        "atoms": large_system["atoms"],
        "prmtop_mb": round(len(large_system["prmtop"]) / 2 ** 20, 1),
        "upload_and_build_s": round(seconds, 1),
        "rss_over_idle_mb": stats["rss_peak_over_baseline_mb"],
    }])


# ---------------------------------------------------------------------------
# 7. Disk - kolik zabere jedna session a jestli ho úklid vrátí
# ---------------------------------------------------------------------------

def test_disk_per_session_and_cleanup(live_server, load_report, monkeypatch):
    from app.workspaces.tasks import garbage_collector as gc

    sizes = {}
    for result in _sessions:
        path = live_server.workspace_dir / result.workspace_id
        if path.exists():
            sizes.setdefault(result.name, []).append(dir_size_mb(path))
    total_mb = dir_size_mb(live_server.workspace_dir)
    count = sum(1 for _ in live_server.workspace_dir.iterdir())
    rows = [{"system": name, "sessions": len(v), "mb_per_session": round(sum(v) / len(v), 1)}
            for name, v in sorted(sizes.items())]
    rows.append({"system": "all workspaces now", "sessions": count, "mb_per_session": f"{total_mb:.0f} MB total"})

    # Úklid: workspace starší než MAX_AGE_SECONDS smaže (na serveru běží jednou
    # za CLEANUP_INTERVAL) - tady ručně, se "zestárlými" složkami.
    old = time.time() - gc.MAX_AGE_SECONDS - 60
    import os
    for item in live_server.workspace_dir.iterdir():
        os.utime(item, (old, old))
    monkeypatch.setattr(gc, "WORKSPACE_DIR", str(live_server.workspace_dir))
    deleted = gc._cleanup_old_workspaces_sync()
    assert deleted == count
    assert not any(live_server.workspace_dir.iterdir())
    load_report.add("Disk per session and cleanup", rows,
                    note=f"Workspaces are deleted {gc.MAX_AGE_SECONDS // 3600} h after the last change, the cleanup "
                         f"runs every {gc.CLEANUP_INTERVAL // 60} min - so a workspace can stay up to "
                         f"{(gc.MAX_AGE_SECONDS + gc.CLEANUP_INTERVAL) // 3600} h. Cleanup removed {deleted} "
                         f"workspaces ({total_mb:.0f} MB).")
    # Incidenty: bez nahlášené chyby se nic neukládá mimo workspace.
    assert dir_size_mb(live_server.incidents_dir) < 1.0
    httpx.get(f"{live_server.base_url}/api/health", timeout=5).raise_for_status()  # server přežil celý běh
