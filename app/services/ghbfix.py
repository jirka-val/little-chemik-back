"""
gHBfix a NBfix korekce z šestého souboru FF v IDA (`ghbfix_nbfix_ff_file`).

Soubor je tabulka párů typů atomů se sloupci
`Type1 Type2 R eps eta rbeg rend c`. Většina řádků je nulová; nenulové
řádky jsou dvojího druhu:

- NBfix (R, eps): upravené LJ parametry pro daný pár typů - patří do
  topologie (prmtop), ne do mdin. Zatím se nepoužívají
  (viz MEETING_PLAN_2026-10.md, fáze B).
- gHBfix (eta, rbeg, rend, c): korekce vodíkových vazeb, kterou simulace
  načítá jako restrainty - mdin dostane `nmropt=1`, blok `&wt` a
  `DISANG=ghbfix.f` (render_amber_mdin), soubor ghbfix.f se exportuje
  spolu s topologií.

Formát ghbfix.f zatím není ověřený: obsahuje nenulové gHBfix řádky ve
stejném sloupcovém formátu jako FF z IDA. Ke konzultaci, viz plán.
"""

import base64
import logging
from functools import lru_cache
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

logger = logging.getLogger(__name__)

GHBFIX_FILENAME = "ghbfix.f"
GHBFIX_LISTOUT = "ghbfix.out"

_HEADER = ";Type1 Type2   R       eps      eta      rbeg     rend     c"


class CorrectionRows(NamedTuple):
    ghbfix: Tuple[str, ...]
    nbfix: Tuple[str, ...]


@lru_cache(maxsize=128)
def _parse(raw: str) -> CorrectionRows:
    try:
        text = base64.b64decode(raw).decode("utf-8", errors="replace")
    except Exception:
        text = raw
    ghbfix: List[str] = []
    nbfix: List[str] = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 8 or line.lstrip().startswith(";"):
            continue
        try:
            r, eps, eta, rbeg, rend, c = (float(v) for v in parts[2:8])
        except ValueError:
            continue
        if r or eps:
            nbfix.append(line.rstrip())
        if eta or rbeg or rend or c:
            ghbfix.append(line.rstrip())
    return CorrectionRows(tuple(ghbfix), tuple(nbfix))


def correction_rows(ff_data: Dict[str, Any]) -> CorrectionRows:
    """Nenulové gHBfix a NBfix řádky jednoho FF (prázdné, když FF soubor nemá)."""
    raw = ff_data.get("ghbfix_nbfix_ff_file")
    if not raw or not isinstance(raw, str):
        return CorrectionRows((), ())
    return _parse(raw)


def correction_flags(ff_data: Dict[str, Any]) -> Dict[str, bool]:
    """`has_ghbfix`/`has_nbfix` pro odpověď FF endpointu."""
    rows = correction_rows(ff_data)
    return {"has_ghbfix": bool(rows.ghbfix), "has_nbfix": bool(rows.nbfix)}


def render_ghbfix_file(ff_selections: Dict[str, Any]) -> Optional[str]:
    """Obsah ghbfix.f pro vybrané FF, nebo None, když žádný gHBfix nemá."""
    sections = []
    for ff_data in ff_selections.values():
        rows = correction_rows(ff_data)
        if not rows.ghbfix:
            continue
        name = ff_data.get("display_name") or ff_data.get("ff_name") or "unknown_ff"
        sections.append(f"; gHBfix terms of {name} ({len(rows.ghbfix)} type pairs)\n" + "\n".join(rows.ghbfix))
    if not sections:
        return None
    return (
        "; gHBfix restraint definitions exported by FORGE.\n"
        "; Non-zero gHBfix rows of the selected force fields, in the column\n"
        "; format of the force field's ghbfix_nbfix file.\n"
        f"{_HEADER}\n" + "\n".join(sections) + "\n"
    )
