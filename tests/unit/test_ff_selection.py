"""
Výběr silových polí (app/services/ff_selection.py): frontend dostává FF bez
obsahů souborů a backend je bere z katalogu, ne z požadavku.
"""

import pytest

from app.core.exceptions import BadRequestError
from app.services.ff_selection import FF_FILE_KEYS, public_ff, resolve_ff_selections

pytestmark = pytest.mark.unit


def _ff(ff_id, name, types, **files):
    return {"id": ff_id, "ff_name": name, "display_name": name, "molecule_type": types,
            **{key: files.get(key, f"{name}:{key}") for key in FF_FILE_KEYS}}


CATALOG = [
    _ff(1, "FF14SB", ["P"]),
    _ff(2, "OL3", ["R"]),
    _ff(3, "FF99", ["P"]),
    _ff(4, "FF99", ["R"]),
    _ff(5, "OPC", ["W4"]),
]


def test_public_ff_drops_file_contents_only():
    ff = {**CATALOG[0], "tier": "recommended"}
    out = public_ff(ff)
    assert not set(FF_FILE_KEYS) & set(out)
    assert out["tier"] == "recommended" and out["id"] == 1


def test_files_come_from_the_catalog_not_from_the_request():
    sent = {"id": 1, "display_name": "../../evil", "bonded_ff_file": "INJECTED", "tier": "recommended",
            "has_ghbfix": True}
    resolved = resolve_ff_selections({"P": sent}, CATALOG)["P"]
    assert resolved["bonded_ff_file"] == "FF14SB:bonded_ff_file"
    assert resolved["display_name"] == "FF14SB"        # žádná cesta z požadavku
    assert resolved["tier"] == "recommended"            # doplněk z FF endpointu zůstává
    assert resolved["has_ghbfix"] is False              # spočítáno z katalogu


def test_stripped_selection_from_the_frontend_is_completed():
    resolved = resolve_ff_selections({"P": public_ff(CATALOG[0]), "W": public_ff(CATALOG[4])}, CATALOG)
    assert resolved["P"]["residue_lib_ff_file"] == "FF14SB:residue_lib_ff_file"
    assert resolved["W"]["ff_name"] == "OPC"


def test_name_without_id_is_resolved_by_group():
    # FF99 je v katalogu dvakrát (protein a RNA) - rozhodne skupina
    assert resolve_ff_selections({"R": {"display_name": "FF99"}}, CATALOG)["R"]["id"] == 4
    assert resolve_ff_selections({"R": {"display_name": "OL3"}}, CATALOG)["R"]["id"] == 2


@pytest.mark.parametrize("selection", [
    {"P": {"id": 999}},                       # neznámé id
    {"P": {"display_name": "NOPE"}},          # neznámý název
    {"D": {"display_name": "FF99"}},          # nejednoznačný název mimo skupiny
    {"P": "FF14SB"},                          # není objekt
])
def test_unknown_force_field_is_rejected(selection):
    with pytest.raises(BadRequestError) as exc:
        resolve_ff_selections(selection, CATALOG)
    assert exc.value.code == "unknown_force_field"
