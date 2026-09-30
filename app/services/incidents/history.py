"""
Deník akcí a snímky struktury uvnitř workspace (podklad pro hlášení chyb).

Každý workspace si vede složku `_history/`:
  journal.jsonl   - jeden řádek na API požadavek (čas, metoda, cesta, JSON
                    nastavení z těla požadavku, status, doba trvání)
  errors.jsonl    - serverové chyby včetně tracebacku
  snapshots/      - gzip kopie structure.pdb před každým krokem, který ji mění
                    (000_upload_structure.pdb.gz, 001_prepare_structure.pdb.gz...)

Nic z toho neopouští workspace, dokud uživatel se souhlasem neodešle report
(viz store.py) - a s workspace to smaže i garbage collector.
"""

from __future__ import annotations

import gzip
import json
import logging
import re
import shutil
import threading
import time
from pathlib import Path
from typing import Any, Dict, Optional

from app.core.config import settings
from app.workspaces.manager import WORKSPACE_DIR

logger = logging.getLogger(__name__)

HISTORY_DIRNAME = "_history"
JOURNAL_FILE = "journal.jsonl"
ERRORS_FILE = "errors.jsonl"
SNAPSHOT_DIRNAME = "snapshots"
_SNAPSHOT_STATE_FILE = "snapshot_state.json"

# Jen hlavní struktura - structure_preview.pdb se mění při každém posunu
# slideru v side-chain GUI a do reportu se přidá jen její aktuální stav.
SNAPSHOT_FILES = ("structure.pdb",)

# Velké JSON tělo (celý sidechain payload apod.) se do deníku zkrátí a nad
# tuhle velikost deníku se už těla neukládají vůbec, jen cesta a status.
_MAX_BODY_CHARS = 200_000
_MAX_JOURNAL_BYTES = 20 * 1024 * 1024

_UUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
_lock = threading.Lock()


def enabled() -> bool:
    return bool(settings.INCIDENT_REPORTS_ENABLED)


def workspace_path(workspace_id: str) -> Optional[Path]:
    """Adresář existujícího workspace, nebo None (nikdy ho nezakládá)."""
    if not workspace_id or not _UUID_RE.fullmatch(workspace_id):
        return None
    path = Path(WORKSPACE_DIR) / workspace_id
    return path if path.is_dir() else None


def history_dir(workspace_id: str) -> Optional[Path]:
    ws = workspace_path(workspace_id)
    if ws is None:
        return None
    path = ws / HISTORY_DIRNAME
    path.mkdir(exist_ok=True)
    return path


def workspace_id_from(path: str, body: Any) -> Optional[str]:
    match = _UUID_RE.search(path)
    if match:
        return match.group(0)
    if isinstance(body, dict) and isinstance(body.get("workspace_id"), str):
        return body["workspace_id"]
    return None


def _append(workspace_id: str, filename: str, entry: Dict[str, Any]) -> None:
    hdir = history_dir(workspace_id)
    if hdir is None:
        return
    line = json.dumps(entry, ensure_ascii=False, default=str)
    with _lock, open(hdir / filename, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def compact_body(body: Any, raw_size: int, journal_path: Optional[Path]) -> Any:
    if body is None:
        return None
    if journal_path is not None and journal_path.exists() and journal_path.stat().st_size > _MAX_JOURNAL_BYTES:
        return f"[omitted - journal over {_MAX_JOURNAL_BYTES // (1024 * 1024)} MB]"
    if raw_size > _MAX_BODY_CHARS:
        return f"[omitted - {raw_size} bytes]"
    return body


def record_request(workspace_id: str, entry: Dict[str, Any]) -> None:
    if not enabled():
        return
    hdir = history_dir(workspace_id)
    if hdir is None:
        return
    entry = dict(entry)
    entry["body"] = compact_body(entry.get("body"), entry.pop("body_size", 0), hdir / JOURNAL_FILE)
    _append(workspace_id, JOURNAL_FILE, entry)


def record_error(workspace_id: str, entry: Dict[str, Any]) -> None:
    if not enabled():
        return
    _append(workspace_id, ERRORS_FILE, {"time": _now(), **entry})


def snapshot_structure(workspace_id: str, label: str) -> None:
    """
    Uloží gzip kopii souborů ze SNAPSHOT_FILES, pokud se od posledního snímku
    změnily (velikost/mtime) - volá se PŘED krokem, který je může přepsat.
    """
    if not enabled():
        return
    ws = workspace_path(workspace_id)
    if ws is None:
        return
    hdir = history_dir(workspace_id)
    snap_dir = hdir / SNAPSHOT_DIRNAME
    snap_dir.mkdir(exist_ok=True)
    safe_label = re.sub(r"[^A-Za-z0-9]+", "-", label).strip("-")[:40] or "step"

    with _lock:
        state_path = hdir / _SNAPSHOT_STATE_FILE
        try:
            state = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            state = {"seq": 0, "files": {}}
        changed = False
        for name in SNAPSHOT_FILES:
            source = ws / name
            if not source.is_file():
                continue
            stat = source.stat()
            signature = [stat.st_size, stat.st_mtime_ns]
            if state["files"].get(name) == signature:
                continue
            target = snap_dir / f"{state['seq']:03d}_{safe_label}_{name}.gz"
            try:
                with open(source, "rb") as src, gzip.open(target, "wb", compresslevel=6) as dst:
                    shutil.copyfileobj(src, dst)
            except OSError as exc:
                logger.warning(f"Incident snapshot of {source} failed: {exc}")
                continue
            state["files"][name] = signature
            state["seq"] += 1
            changed = True
        if changed:
            state_path.write_text(json.dumps(state), encoding="utf-8")


def record_workspace_created(workspace_id: str, source: str) -> None:
    """Volá se z endpointů, které workspace zakládají (upload, stažení z RCSB)."""
    if not enabled():
        return
    _append(workspace_id, JOURNAL_FILE, {"time": _now(), "event": "workspace_created", "source": source})
    snapshot_structure(workspace_id, "upload")


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())
