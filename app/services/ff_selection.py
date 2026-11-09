"""
Výběr silových polí mezi frontendem a backendem.

Dřív GET /api/forcefields/{ws} posílal celé obsahy FF souborů (~7 MB) a
frontend je při přípravě i topologii posílal zpátky v `ff_selections`
(~1,7 MB na požadavek); backend z nich psal sdílenou FF cache a stavěl
topologii. Tedy: zbytečný přenos a hlavně důvěra klientovi - upravený
požadavek mohl podstrčit vlastní "silové pole" do cache všech uživatelů a z
`display_name` se skládá cesta ke složce cache.

Teď:
- public_ff() - odpověď FF endpointu bez obsahů souborů,
- resolve_ff_selections() - na vstupu se každé vybrané FF najde v katalogu
  (podle `id`, jinak podle názvu a skupiny) a obsah souborů, názvy i
  příznaky korekcí se vezmou z katalogu. Od klienta zůstanou jen doplňky
  z FF endpointu (tier, is_default ...). FF, které v katalogu není, se odmítne.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from app.core.exceptions import BadRequestError
from app.services.ff_catalog_service import catalog_service
from app.services.ghbfix import correction_flags

# Obsahy souborů z IDA - frontend je nepotřebuje, backend je bere z katalogu.
FF_FILE_KEYS = (
    "force_field_file", "nonbonded_ff_file", "bonded_ff_file", "atom_type_ff_file",
    "residue_lib_ff_file", "ghbfix_nbfix_ff_file",
)
# Spočítá server z obsahu souborů, klientovi se nevěří.
_DERIVED_KEYS = ("has_ghbfix", "has_nbfix")


def public_ff(ff: Dict[str, Any]) -> Dict[str, Any]:
    """FF objekt pro frontend: vše kromě obsahů souborů."""
    return {k: v for k, v in ff.items() if k not in FF_FILE_KEYS}


def _name(ff: Dict[str, Any]) -> Optional[str]:
    return ff.get("display_name") or ff.get("ff_name")


def _groups_match(ff: Dict[str, Any], group: str) -> bool:
    types = ff.get("molecule_type") or []
    # Obecné "W" (voda) a iontové skupiny se v katalogu vedou jako W3/W4/W5, I1 ...
    return group in types or (group == "W" and any(str(t).startswith("W") for t in types))


def _find_in_catalog(group: str, selected: Dict[str, Any], catalog: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    ff_id = selected.get("id")
    if ff_id is not None:
        return next((ff for ff in catalog if ff.get("id") == ff_id), None)
    name = _name(selected)
    if not name:
        return None
    by_name = [ff for ff in catalog if _name(ff) == name]
    for_group = [ff for ff in by_name if _groups_match(ff, group)]
    candidates = for_group or by_name
    return candidates[0] if len(candidates) == 1 else None


def resolve_ff_selections(ff_selections: Dict[str, Any],
                          catalog: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """{skupina: FF od klienta} -> {skupina: FF z katalogu (+ doplňky klienta)}.

    Vyhodí BadRequestError, když FF v katalogu není (neznámé id / název,
    nebo nejednoznačný název bez id).
    """
    catalog = catalog if catalog is not None else catalog_service.get_forcefields()
    resolved: Dict[str, Any] = {}
    for group, selected in (ff_selections or {}).items():
        if not isinstance(selected, dict):
            raise BadRequestError(f"Force field for group {group} is not an object.", code="unknown_force_field")
        entry = _find_in_catalog(group, selected, catalog)
        if entry is None:
            label = selected.get("id") or _name(selected) or "?"
            raise BadRequestError(
                f"Force field {label} for group {group} is not in the current catalog - reload the page "
                "and select the force fields again.", code="unknown_force_field")
        extras = {k: v for k, v in selected.items()
                  if k not in entry and k not in FF_FILE_KEYS and k not in _DERIVED_KEYS}
        resolved[group] = {**extras, **entry, **correction_flags(entry)}
    return resolved
