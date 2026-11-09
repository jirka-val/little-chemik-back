"""
Fixtures pro zátěžové testy: skutečný server (uvicorn v samostatném procesu,
jako v Dockeru - jeden proces, bez reloadu) s vlastní složkou workspaces a
incidentů v tmp, a report, do kterého testy zapisují naměřená čísla
(tests/load/results/latest.md + .json).
"""

from __future__ import annotations

import json
import os
import platform
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import psutil
import pytest

from tests.load.helpers import LiveServer, ServerStartError, ff_data_available, start_server

RESULTS_DIR = Path(__file__).resolve().parent / "results"


@pytest.fixture(scope="session")
def live_server(tmp_path_factory) -> LiveServer:
    if not ff_data_available():
        pytest.skip("Zátěžové testy potřebují stažený FF katalog a data/ff_cache (viz tests/README.md).")
    workers = os.environ.get("LOAD_WORKERS")  # volitelně: srovnání s víc procesy uvicornu
    try:
        with start_server(tmp_path_factory.mktemp("load_server"), int(workers) if workers else None) as server:
            yield server
    except ServerStartError as exc:
        pytest.fail(str(exc))


class LoadReport:
    """Sekce s tabulkami a poznámkami; na konci session se zapíše do results/."""

    def __init__(self):
        self.sections: List[Dict[str, Any]] = []

    def add(self, title: str, rows: List[Dict[str, Any]], note: str = "") -> None:
        self.sections.append({"title": title, "rows": rows, "note": note})

    @staticmethod
    def _table(rows: List[Dict[str, Any]]) -> str:
        if not rows:
            return "_(no data)_\n"
        cols = list(dict.fromkeys(k for row in rows for k in row))
        out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        out += ["| " + " | ".join(str(row.get(c, "")) for c in cols) + " |" for row in rows]
        return "\n".join(out) + "\n"

    def write(self, machine: Dict[str, Any]) -> Path:
        RESULTS_DIR.mkdir(exist_ok=True)
        lines = ["# FORGE load test results", "",
                 "Machine: " + ", ".join(f"{k} {v}" for k, v in machine.items()), ""]
        for section in self.sections:
            lines += [f"## {section['title']}", ""]
            if section["note"]:
                lines += [section["note"], ""]
            lines.append(self._table(section["rows"]))
        text = "\n".join(lines)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        (RESULTS_DIR / "latest.md").write_text(text, encoding="utf-8")
        (RESULTS_DIR / f"{stamp}.json").write_text(
            json.dumps({"machine": machine, "sections": self.sections}, indent=1, default=str), encoding="utf-8")
        return RESULTS_DIR / "latest.md"


@pytest.fixture(scope="session")
def load_report():
    report = LoadReport()
    yield report
    if report.sections:
        machine = {
            "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "os": platform.platform(terse=True),
            "cpu": f"{psutil.cpu_count(logical=False)} cores / {psutil.cpu_count()} threads",
            "ram_gb": round(psutil.virtual_memory().total / 2 ** 30, 1),
            "python": platform.python_version(),
            "uvicorn_workers": os.environ.get("LOAD_WORKERS", "1"),
        }
        path = report.write(machine)
        print(f"\nLoad test report: {path}")
