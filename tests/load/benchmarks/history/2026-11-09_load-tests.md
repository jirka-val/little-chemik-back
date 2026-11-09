# Load tests 2026-11-09 (before the benchmark existed)

Output of tests/load/test_capacity.py on the developer PC (Windows 11, 24 cores / 32 threads, 31.7 GB RAM).
Three runs: production code (commit 75a93ea, 1 worker), then after the fixes of the same day (FF cache race, topology off the event loop) with 1 and 4 uvicorn workers.

---

# Run A: production code, 1 worker

Machine: date 2026-11-09 11:19, os Windows-11, cpu 24 cores / 32 threads, ram_gb 31.7, python 3.12.10, uvicorn_workers 1

## Single user: 1CRN.pdb (small protein, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 1.0 s; `POST /api/topology/{ws}/generate` 0.2 s; `GET /api/analysis/sequence/{ws}` 0.1 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 9442 | 1.7 | 1.6 | 248 | 62 | 4.7 | 14 |

## Single user: HIS_1CA2_carbonic_anhydrase_zinc.pdb (medium protein + Zn, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 3.2 s; `POST /api/topology/{ws}/generate` 0.7 s; `GET /api/analysis/sequence/{ws}` 0.3 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 29675 | 5.0 | 4.7 | 304 | 68 | 11.6 | 14 |

## Single user: 3SKR.pdb (RNA, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 4.3 s; `POST /api/topology/{ws}/generate` 1.0 s; `GET /api/analysis/sequence/{ws}` 0.4 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 36078 | 6.4 | 6.0 | 320 | 68 | 41.6 | 14 |

## Single user: 3SKR.pdb (RNA, 3 nm water, 3.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 7.3 s; `POST /api/topology/{ws}/generate` 2.6 s; `GET /api/analysis/sequence/{ws}` 1.4 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 124261 | 12.8 | 12.2 | 583 | 319 | 67.9 | 14 |

## Single user: 3SKR.pdb (RNA, 5 nm water, 5.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 12.9 s; `POST /api/topology/{ws}/generate` 5.8 s; `GET /api/analysis/sequence/{ws}` 3.1 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 295692 | 25.2 | 24.2 | 1052 | 761 | 118.8 | 14 |

## Concurrent users (HIS_1CA2_carbonic_anhydrase_zinc.pdb, all start at once)

One user alone needs 5.0 s. `slowdown_p95` = how much longer the slowest users wait. `cpu_mean_cores` = average number of busy cores of the server process.

| users | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | throughput_per_min | cpu_mean_cores | rss_peak_mb | identical_prmtop |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 4.8 | 4.826 | 4.826 | 4.826 | 1.0x | 12.4 | 0.91 | 326 | True |
| 2 | 10.6 | 10.543 | 10.563 | 10.565 | 2.1x | 11.3 | 0.99 | 351 | False |
| 4 | 23.4 | 23.307 | 23.339 | 23.339 | 4.6x | 10.3 | 1.06 | 416 | False |
| 8 | 47.1 | 47.041 | 47.063 | 47.066 | 9.3x | 10.2 | 1.08 | 699 | False |
| 16 | 95.8 | 95.708 | 95.759 | 95.781 | 19.0x | 10.0 | 1.09 | 934 | False |

## Light requests while 8 users prepare

If the health check is slow too, heavy work blocks the event loop and the whole site freezes.

| request | n | p50_s | p95_s | max_s |
|---|---|---|---|---|
| GET /api/health | 93 | 0.002 | 0.374 | 2.486 |
| GET /api/analysis/sequence/{ws} | 93 | 0.018 | 0.99 | 3.009 |

## Users arriving every 6 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 600 | 88.9 | 1.613 | 6.568 | 7.2 | 1.0x | 0.52 | 425 |

## Users arriving every 2 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 1800 | 53.3 | 16.069 | 38.914 | 42.963 | 8.8x | 1.08 | 419 |

## 50 users load a structure at the same moment (upload + sequence + altlocs + viewer)

| requests | wall_s | p50_s | p95_s | max_s | requests_per_s | rss_peak_mb |
|---|---|---|---|---|---|---|
| 250 | 3.2 | 0.076 | 0.159 | 0.224 | 77.2 | 367 |

## Large system from topology + coordinates (no PDB)

| atoms | prmtop_mb | upload_and_build_s | rss_over_idle_mb |
|---|---|---|---|
| 124261 | 19.3 | 1.9 | 34 |

## Disk per session and cleanup

Workspaces are deleted 2 h after the last change, the cleanup runs every 60 min - so a workspace can stay up to 3 h. Cleanup removed 128 workspaces (1174 MB).

| system | sessions | mb_per_session |
|---|---|---|
| 1CRN.pdb | 17 | 4.7 |
| 3SKR.pdb | 7 | 56.4 |
| HIS_1CA2_carbonic_anhydrase_zinc.pdb | 50 | 11.6 |
| all workspaces now | 128 | 1174 MB total |

---

# Run B: after fixes, 1 worker (subset)

Machine: date 2026-11-09 11:44, os Windows-11, cpu 24 cores / 32 threads, ram_gb 31.7, python 3.12.10, uvicorn_workers 1

## Single user: 1CRN.pdb (small protein, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 0.9 s; `POST /api/topology/{ws}/generate` 0.2 s; `GET /api/analysis/sequence/{ws}` 0.1 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 9442 | 1.7 | 1.5 | 248 | 61 | 4.7 | 14 |

## Single user: HIS_1CA2_carbonic_anhydrase_zinc.pdb (medium protein + Zn, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 3.2 s; `POST /api/topology/{ws}/generate` 0.7 s; `GET /api/analysis/sequence/{ws}` 0.3 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 29675 | 4.9 | 4.6 | 309 | 72 | 11.6 | 14 |

## Single user: 3SKR.pdb (RNA, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 4.2 s; `POST /api/topology/{ws}/generate` 1.0 s; `GET /api/analysis/sequence/{ws}` 0.4 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 36078 | 6.2 | 5.8 | 324 | 73 | 41.6 | 14 |

## Single user: 3SKR.pdb (RNA, 3 nm water, 3.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 7.2 s; `POST /api/topology/{ws}/generate` 2.5 s; `GET /api/analysis/sequence/{ws}` 1.3 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 124261 | 12.6 | 11.9 | 576 | 312 | 67.9 | 14 |

## Single user: 3SKR.pdb (RNA, 5 nm water, 5.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 12.8 s; `POST /api/topology/{ws}/generate` 5.7 s; `GET /api/analysis/sequence/{ws}` 3.1 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 295692 | 24.8 | 23.8 | 1042 | 752 | 118.8 | 14 |

## Concurrent users (HIS_1CA2_carbonic_anhydrase_zinc.pdb, all start at once)

One user alone needs 4.9 s. `slowdown_p95` = how much longer the slowest users wait. `cpu_mean_cores` = average number of busy cores of the server process.

| users | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | throughput_per_min | cpu_mean_cores | rss_peak_mb | identical_prmtop |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 4.8 | 4.809 | 4.809 | 4.809 | 1.0x | 12.5 | 0.92 | 331 | True |
| 2 | 10.3 | 10.286 | 10.305 | 10.307 | 2.1x | 11.6 | 1.0 | 360 | True |
| 4 | 23.6 | 23.575 | 23.595 | 23.598 | 4.8x | 10.2 | 1.07 | 483 | True |
| 8 | 47.0 | 46.908 | 46.987 | 46.988 | 9.6x | 10.2 | 1.09 | 502 | True |
| 16 | 95.5 | 95.424 | 95.488 | 95.489 | 19.6x | 10.1 | 1.09 | 874 | True |

## Light requests while 8 users prepare

If the health check is slow too, heavy work blocks the event loop and the whole site freezes.

| request | n | p50_s | p95_s | max_s |
|---|---|---|---|---|
| GET /api/health | 102 | 0.002 | 0.212 | 0.872 |
| GET /api/analysis/sequence/{ws} | 102 | 0.019 | 0.775 | 1.825 |

## Users arriving every 6 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 600 | 88.9 | 1.629 | 7.227 | 7.276 | 1.2x | 0.51 | 421 |

## Users arriving every 2 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 1800 | 53.6 | 12.768 | 39.136 | 43.168 | 8.0x | 1.09 | 423 |

## 50 users load a structure at the same moment (upload + sequence + altlocs + viewer)

| requests | wall_s | p50_s | p95_s | max_s | requests_per_s | rss_peak_mb |
|---|---|---|---|---|---|---|
| 250 | 3.3 | 0.07 | 0.17 | 0.258 | 75.2 | 371 |

---

# Run C: after fixes, 4 workers (subset)

Machine: date 2026-11-09 11:48, os Windows-11, cpu 24 cores / 32 threads, ram_gb 31.7, python 3.12.10, uvicorn_workers 4

## Single user: 1CRN.pdb (small protein, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 0.8 s; `POST /api/topology/{ws}/generate` 0.3 s; `GET /api/analysis/sequence/{ws}` 0.1 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 9442 | 1.7 | 1.8 | 851 | 77 | 4.7 | 14 |

## Single user: HIS_1CA2_carbonic_anhydrase_zinc.pdb (medium protein + Zn, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 3.2 s; `POST /api/topology/{ws}/generate` 0.7 s; `GET /api/analysis/sequence/{ws}` 0.3 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 29675 | 5.0 | 4.8 | 938 | 105 | 11.6 | 14 |

## Single user: 3SKR.pdb (RNA, 1.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 4.2 s; `POST /api/topology/{ws}/generate` 0.9 s; `GET /api/analysis/sequence/{ws}` 0.4 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 36078 | 6.3 | 6.0 | 1012 | 127 | 41.6 | 14 |

## Single user: 3SKR.pdb (RNA, 3 nm water, 3.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 7.1 s; `POST /api/topology/{ws}/generate` 2.5 s; `GET /api/analysis/sequence/{ws}` 1.4 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 124261 | 12.6 | 12.0 | 1286 | 337 | 67.9 | 14 |

## Single user: 3SKR.pdb (RNA, 5 nm water, 5.0 nm water)

Slowest steps: `POST /api/sidechains/start/{ws}` 12.7 s; `POST /api/topology/{ws}/generate` 5.7 s; `GET /api/analysis/sequence/{ws}` 3.3 s

| atoms | session_s | server_cpu_s | rss_peak_mb | rss_over_idle_mb | workspace_mb | requests |
|---|---|---|---|---|---|---|
| 295692 | 25.1 | 24.1 | 1862 | 866 | 118.8 | 14 |

## Concurrent users (HIS_1CA2_carbonic_anhydrase_zinc.pdb, all start at once)

One user alone needs 5.0 s. `slowdown_p95` = how much longer the slowest users wait. `cpu_mean_cores` = average number of busy cores of the server process.

| users | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | throughput_per_min | cpu_mean_cores | rss_peak_mb | identical_prmtop |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 4.9 | 4.901 | 4.901 | 4.901 | 1.0x | 12.2 | 0.9 | 1127 | True |
| 2 | 10.2 | 10.143 | 10.161 | 10.163 | 2.0x | 11.8 | 1.02 | 1163 | True |
| 4 | 10.7 | 8.036 | 10.727 | 10.735 | 2.2x | 22.3 | 1.91 | 1207 | True |
| 8 | 18.4 | 11.225 | 18.336 | 18.34 | 3.7x | 26.2 | 2.58 | 1277 | True |
| 16 | 37.8 | 31.783 | 37.796 | 37.798 | 7.6x | 25.4 | 2.84 | 1539 | True |

## Light requests while 8 users prepare

If the health check is slow too, heavy work blocks the event loop and the whole site freezes.

| request | n | p50_s | p95_s | max_s |
|---|---|---|---|---|
| GET /api/health | 59 | 0.001 | 0.05 | 0.204 |
| GET /api/analysis/sequence/{ws} | 59 | 0.005 | 0.045 | 0.282 |

## Users arriving every 6 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 600 | 89.0 | 1.619 | 6.504 | 7.206 | 1.0x | 0.52 | 1328 |

## Users arriving every 2 s (15 users, mixed systems)

Mix: 50 % small protein, 33 % medium protein, 17 % RNA. `slowdown_p95` compares each session with the same system run alone.

| rate_per_hour | wall_s | session_p50_s | session_p95_s | session_max_s | slowdown_p95 | cpu_mean_cores | rss_peak_mb |
|---|---|---|---|---|---|---|---|
| 1800 | 33.0 | 2.159 | 6.984 | 7.866 | 1.3x | 1.42 | 1329 |

## 50 users load a structure at the same moment (upload + sequence + altlocs + viewer)

| requests | wall_s | p50_s | p95_s | max_s | requests_per_s | rss_peak_mb |
|---|---|---|---|---|---|---|
| 250 | 3.0 | 0.019 | 0.056 | 0.073 | 82.5 | 1289 |
