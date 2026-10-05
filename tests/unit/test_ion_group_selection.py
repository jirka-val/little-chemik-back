"""
FF panel nabízí všechny iontové skupiny (I1/I1+/Im/Im+), aby změna soli v
kroku Structure nevyžadovala nový výběr FF. Builder a topologie ale dostanou
jen skupiny, které struktura a soli opravdu použijí.
"""

import pytest

from app.api.v1.endpoints.forcefields import _with_all_ion_groups
from app.services.structure.groups import drop_unused_ion_groups

pytestmark = pytest.mark.unit


def test_all_ion_groups_offered_needed_first():
    required = {"P": {"reason": "protein"}, "W": {"reason": "solvation"}, "Im": {"reason": "ions: Mg2+ (x3)"}}
    out = _with_all_ion_groups(required)
    assert list(out) == ["P", "W", "Im", "I1", "I1+", "Im+"]
    assert out["Im"] == {"reason": "ions: Mg2+ (x3)", "needed": True}
    assert out["P"]["needed"] is True
    assert out["I1"]["needed"] is False
    assert out["Im+"]["needed"] is False


def test_no_ion_groups_without_solvent():
    out = _with_all_ion_groups({"P": {"reason": "protein"}})
    assert list(out) == ["P"]


def test_drop_unused_ion_groups_keeps_non_ion_entries():
    selections = {"P": "ff14sb", "W": "spce", "I1": "jc", "I1+": "jc+", "Im": "lm", "Im+": "lm+", "I": "legacy"}
    out = drop_unused_ion_groups(selections, {"P", "W3", "I1", "Im"})
    assert out == {"P": "ff14sb", "W": "spce", "I1": "jc", "Im": "lm", "I": "legacy"}
