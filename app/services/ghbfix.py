"""
gHBfix a NBfix korekce z šestého souboru FF v IDA (`ghbfix_nbfix_ff_file`).

Soubor je tabulka párů typů atomů se sloupci
`Type1 Type2 R eps eta rbeg rend c`. Většina řádků je nulová; nenulové
řádky jsou dvojího druhu:

- NBfix (R, eps): upravené LJ parametry pro daný pár typů - patří do
  topologie (prmtop), ne do mdin. FF instance je má načtené v `ff.nbfix`
  (FF_IDA.ff, parametr ghbfix_nbfix_file) a create_AMBER_topology je
  promítá do LENNARD_JONES_ACOEF/BCOEF.
- gHBfix (eta, rbeg, rend, c): korekce vodíkových vazeb, kterou simulace
  načítá jako restrainty - mdin dostane `nmropt=1`, blok `&wt` a
  `DISANG=ghbfix.f` (render_amber_mdin). Soubor ghbfix.f (bloky `&rst`
  pro jednotlivé páry atomů) zapisuje create_AMBER_topology
  (write_ghbfix_file) spolu s topologií.

Tady se jen zjišťuje, jestli FF nějaké korekce má (has_ghbfix/has_nbfix
pro FF endpoint a report).
"""

import base64
import logging
from functools import lru_cache
from typing import Any, Dict, List, NamedTuple, Tuple

logger = logging.getLogger(__name__)

GHBFIX_FILENAME = "ghbfix.f"
GHBFIX_LISTOUT = "ghbfix.out"


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
