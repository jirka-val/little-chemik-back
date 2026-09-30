"""
Trvalé úložiště nahlášených chyb (incidentů).

Report vzniká jen na výslovnou žádost uživatele se souhlasem (POST
/api/incidents) - buď po serverové chybě, nebo ručně tlačítkem "!" (chyba,
kterou aplikace nepoznala, ale uživatel ji vidí). Každý incident je složka

    INCIDENTS_DIR/<YYYYMMDD-HHMMSS>-<hex>/
        meta.json          - druh, popis od uživatele, chyba, prohlížeč, verze
        client_log.json    - poslední akce v prohlížeči (kliky, API volání, JS chyby)
        screenshot.png     - snímek 3D vieweru (pokud ho uživatel přiložil)
        journal.jsonl      - deník API požadavků workspace (viz history.py)
        errors.jsonl       - serverové chyby s tracebackem
        snapshots/         - struktura před každým krokem, který ji měnil
        current/           - aktuální soubory workspace v okamžiku reportu

Nic se nemaže automaticky: po dosažení INCIDENTS_MAX_TOTAL_MB se nové
reporty odmítnou a v logu je chyba, ať si toho správce všimne.
"""

from __future__ import annotations

import base64
import gzip
import io
import json
import logging
import re
import secrets
import shutil
import time
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.core.config import settings
from app.services.incidents import history

logger = logging.getLogger(__name__)

_INCIDENT_ID_RE = re.compile(r"^\d{8}-\d{6}-[0-9a-f]{6}$")
_MB = 1024 * 1024


class IncidentStorageFull(Exception):
    pass


def incidents_dir() -> Path:
    path = Path(settings.INCIDENTS_DIR)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _dir_size(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def storage_used_bytes() -> int:
    return _dir_size(incidents_dir())


def _decode_screenshot(data_url: Optional[str]) -> Optional[bytes]:
    if not data_url or not data_url.startswith("data:image/png;base64,"):
        return None
    try:
        return base64.b64decode(data_url.split(",", 1)[1], validate=True)
    except (ValueError, IndexError):
        return None


def create_incident(
    *,
    workspace_id: Optional[str],
    kind: str,
    description: str,
    error: Optional[Dict[str, Any]],
    client_log: List[Any],
    page: Dict[str, Any],
    screenshot: Optional[str],
) -> str:
    used = storage_used_bytes()
    if used >= settings.INCIDENTS_MAX_TOTAL_MB * _MB:
        logger.error(
            f"Incident storage full ({used / _MB:.0f} MB >= {settings.INCIDENTS_MAX_TOTAL_MB} MB) - "
            f"new report rejected. Clean up {incidents_dir()}."
        )
        raise IncidentStorageFull()

    incident_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{secrets.token_hex(3)}"
    target = incidents_dir() / incident_id
    target.mkdir()
    budget = settings.INCIDENT_MAX_MB * _MB
    omitted: List[str] = []

    def write(name: str, data: bytes) -> None:
        nonlocal budget
        (target / name).parent.mkdir(parents=True, exist_ok=True)
        (target / name).write_bytes(data)
        budget -= len(data)

    write("client_log.json", json.dumps(client_log, ensure_ascii=False, indent=1, default=str).encode("utf-8"))
    png = _decode_screenshot(screenshot)
    if png:
        write("screenshot.png", png)

    ws = history.workspace_path(workspace_id) if workspace_id else None
    if ws is not None:
        hdir = ws / history.HISTORY_DIRNAME
        for name in (history.JOURNAL_FILE, history.ERRORS_FILE):
            if (hdir / name).is_file():
                write(name, (hdir / name).read_bytes())

        # Aktuální stav workspace (zabalený) - to, na čem chyba nastala.
        for source in sorted(p for p in ws.iterdir() if p.is_file()):
            packed = gzip.compress(source.read_bytes(), compresslevel=6)
            if len(packed) > budget:
                omitted.append(f"current/{source.name}")
                continue
            write(f"current/{source.name}.gz", packed)

        # Snímky: originál (první) má přednost, pak od nejnovějších.
        snaps = sorted((hdir / history.SNAPSHOT_DIRNAME).glob("*.gz")) if (hdir / history.SNAPSHOT_DIRNAME).is_dir() else []
        ordered = snaps[:1] + list(reversed(snaps[1:]))
        for snap in ordered:
            size = snap.stat().st_size
            if size > budget:
                omitted.append(f"snapshots/{snap.name}")
                continue
            (target / "snapshots").mkdir(exist_ok=True)
            shutil.copyfile(snap, target / "snapshots" / snap.name)
            budget -= size

    meta = {
        "id": incident_id,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "kind": kind,
        "description": description,
        "error": error,
        "workspace_id": workspace_id,
        "workspace_found": ws is not None,
        "page": page,
        "backend_version": settings.VERSION,
        "omitted_over_size_limit": omitted,
    }
    (target / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    logger.warning(
        f"Incident {incident_id} reported ({kind}) for workspace {workspace_id}: "
        f"{(description or (error or {}).get('message') or '')[:200]}"
    )
    return incident_id


def list_incidents() -> List[Dict[str, Any]]:
    result = []
    for path in sorted(incidents_dir().iterdir(), reverse=True):
        meta_file = path / "meta.json"
        if not path.is_dir() or not meta_file.is_file():
            continue
        try:
            meta = json.loads(meta_file.read_text(encoding="utf-8"))
        except ValueError:
            continue
        result.append({
            "id": meta.get("id", path.name),
            "created": meta.get("created"),
            "kind": meta.get("kind"),
            "description": (meta.get("description") or "")[:300],
            "error": (meta.get("error") or {}).get("message"),
            "workspace_id": meta.get("workspace_id"),
            "size_mb": round(_dir_size(path) / _MB, 2),
        })
    return result


def incident_zip(incident_id: str) -> Optional[bytes]:
    if not _INCIDENT_ID_RE.match(incident_id):
        return None
    path = incidents_dir() / incident_id
    if not path.is_dir():
        return None
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for file in path.rglob("*"):
            if file.is_file():
                zf.write(file, f"{incident_id}/{file.relative_to(path).as_posix()}")
    return buffer.getvalue()
