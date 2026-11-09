# Load tests (`tests/load/`, marker `load`)

How much the FORGE server needs: concurrent users, memory, CPU and disk.
Unlike `tests/performance/` (services called in-process), these tests start a
**real server** — `uvicorn app.main:app` in its own process, the same way the
Docker image runs it (one process, no reload) — and talk to it over HTTP
exactly like the frontend does (`tests/golden/pipeline.run_pipeline`).

The server gets its own `WORKSPACE_DIR` and `INCIDENTS_DIR` in pytest's tmp
directory, so `temp_workspaces/` and real incident reports are never touched.
The server runs with its default process pool (`HEAVY_JOB_WORKERS`, see
`app/core/process_pool.py`) unless `--pool` says otherwise; the pytest
`TestClient` tests elsewhere run without it.
It needs the downloaded FF catalog and `data/ff_cache` (see `tests/README.md`);
without them the tests are skipped.

## Running

```bash
# all scenarios (~10-20 min), report in tests/load/results/latest.md
.venv/Scripts/python.exe -m pytest -m load tests/load -s

# one scenario
.venv/Scripts/python.exe -m pytest -m load tests/load -s -k concurrent_users

# compare with several uvicorn processes
LOAD_WORKERS=4 .venv/Scripts/python.exe -m pytest -m load tests/load -s
```

Excluded from the default run (`pytest.ini` addopts). Needs `psutil`
(`requirements-dev.txt`).

## Benchmark (`benchmark.py`) - measure every optimisation

The pytest tests check correctness; the benchmark **measures** with fixed
scenarios (warm-up, 3 repeats of single users, concurrency ladder 1-16,
responsiveness under load, users arriving every 2 s, 50-user light burst) and
saves `tests/load/benchmarks/<date>_<label>_<commit>.json` + `.md` (committed,
so the history stays).

```bash
.venv/Scripts/python.exe -m tests.load.benchmark run --label baseline          # ~7 min
.venv/Scripts/python.exe -m tests.load.benchmark run --label x --workers 4     # several uvicorn processes
.venv/Scripts/python.exe -m tests.load.benchmark run --label x --quick         # ~3 min, rough
.venv/Scripts/python.exe -m tests.load.benchmark run --label x --pool 0        # without the process pool
.venv/Scripts/python.exe -m tests.load.benchmark compare A.json B.json         # KPIs, % change, better/WORSE
```

Compare runs from the same machine and keep the machine otherwise idle. The
baseline and the optimisation plan: `docs/performance-plan.md`.

## Scenarios (`test_capacity.py`)

| test | question |
|---|---|
| `test_single_user_baseline` | What does one preparation cost (time, CPU, RAM, disk) for small → very large systems (9k – ~300k atoms)? |
| `test_concurrent_users` | 1, 2, 4, 8, 16 users press *Prepare* at the same moment: how long do they wait, how many cores are used, how much RAM? Do sessions stay isolated (same input → same output)? |
| `test_light_requests_stay_responsive_under_load` | While 8 users prepare, do cheap requests (health, sequence) still answer, or does the whole site freeze? |
| `test_users_arriving_over_time` | Realistic traffic: users arrive every 6 s / 2 s with a mix of systems. |
| `test_burst_of_light_sessions` | 50 users load a structure at once (upload, sequence, altlocs, viewer). |
| `test_large_pdb_without_cryst1` | Box reconstruction from the water on a large solvated PDB. |
| `test_large_topology_and_coordinates_upload` | Building the structure from a large prmtop + rst7 (no PDB). |
| `test_disk_per_session_and_cleanup` | Disk per session, and whether the workspace cleanup frees it. |

Assertions check **correctness** under load (no errors, isolated sessions,
health check answers, cleanup works). Capacity numbers are only **recorded**
— they depend on the machine — in `tests/load/results/latest.md` (+ a
timestamped `.json`), which is git-ignored.
