"""
Záznamy kroků pro souhrnný report (PDF v Exportu).

Endpointy po úspěšném kroku uloží, co uživatel nastavil, do
`<workspace>/_report/<krok>.json`; report.pdf se z nich skládá až při
exportu. Na rozdíl od deníku pro hlášení chyb (`_history/`) se vede vždy -
je to obsah pro uživatele, ne diagnostika. Smaže se s workspace.
"""

from __future__ import annotations

import json
import logging
from typing import Any, Dict

from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

REPORT_DIRNAME = "_report"

SOURCE = "source"
CONFORMATIONS = "conformations"
PREPARATION = "preparation"
TOPOLOGY = "topology"

# Pozdější kroky stavějí na dřívějších - nový záznam kroku zneplatní ty po něm.
_ORDER = (SOURCE, CONFORMATIONS, PREPARATION, TOPOLOGY)

# Pole FF, která patří do reportu (bez base64 obsahu souborů).
_FF_FIELDS = (
    "ff_name", "display_name", "molecule_type", "doi", "reference_article_doi", "data_publication_time",
)


def ff_summary(ff_data: Dict[str, Any]) -> Dict[str, Any]:
    out = {key: ff_data.get(key) for key in _FF_FIELDS if ff_data.get(key)}
    creators = [c.get("name") for c in ff_data.get("creators") or [] if isinstance(c, dict) and c.get("name")]
    if creators:
        out["creators"] = creators
    for flag in ("tier", "has_ghbfix", "has_nbfix"):
        if flag in ff_data:
            out[flag] = ff_data[flag]
    return out


def save_record(workspace_id: str, step: str, data: Dict[str, Any]) -> None:
    """Uloží záznam kroku a smaže záznamy pozdějších kroků. Chyba reportu nesmí shodit krok."""
    try:
        report_dir = workspace_manager.get_workspace_dir(workspace_id) / REPORT_DIRNAME
        report_dir.mkdir(exist_ok=True)
        (report_dir / f"{step}.json").write_text(json.dumps(data, ensure_ascii=False, default=str), encoding="utf-8")
        for later in _ORDER[_ORDER.index(step) + 1:]:
            (report_dir / f"{later}.json").unlink(missing_ok=True)
    except Exception as e:
        logger.warning(f"Could not save report record '{step}' for workspace {workspace_id}: {e}")


def update_record(workspace_id: str, step: str, data: Dict[str, Any]) -> None:
    """Doplní pole do existujícího záznamu kroku (bez mazání pozdějších)."""
    try:
        path = workspace_manager.get_workspace_dir(workspace_id) / REPORT_DIRNAME / f"{step}.json"
        current = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
        current.update(data)
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps(current, ensure_ascii=False, default=str), encoding="utf-8")
    except Exception as e:
        logger.warning(f"Could not update report record '{step}' for workspace {workspace_id}: {e}")


def load_records(workspace_id: str) -> Dict[str, Dict[str, Any]]:
    report_dir = workspace_manager.get_workspace_dir(workspace_id) / REPORT_DIRNAME
    records: Dict[str, Dict[str, Any]] = {}
    for step in _ORDER:
        path = report_dir / f"{step}.json"
        if not path.exists():
            continue
        try:
            records[step] = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as e:
            logger.warning(f"Unreadable report record {path}: {e}")
    return records
