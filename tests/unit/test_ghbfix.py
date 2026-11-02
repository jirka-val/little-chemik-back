"""
gHBfix/NBfix řádky z `ghbfix_nbfix_ff_file` a jejich promítnutí do
mdin (app/services/ghbfix.py, simulation.render_amber_mdin).
"""

import base64

import pytest

from app.api.v1.endpoints.simulation import AmberMdinRequest, render_amber_mdin
from app.services.ghbfix import correction_flags, correction_rows

pytestmark = pytest.mark.unit

_TABLE = """;Type1 Type2   R       eps      eta      rbeg     rend     c
RAC2  RAC2  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
RAH2  RBO3  2.880800 0.050500 0.000000 0.000000 0.000000 0.000000
RAH6  RAN1  0.000000 0.000000 0.300000 2.000000 3.000000 0.800000
"""


def _ff(name: str, table: str = _TABLE) -> dict:
    return {"ff_name": name, "ghbfix_nbfix_ff_file": base64.b64encode(table.encode()).decode()}


def test_rows_split_into_ghbfix_and_nbfix():
    rows = correction_rows(_ff("OL3CP_gHBfix21_NBfix-0BPh"))
    assert rows.ghbfix == ("RAH6  RAN1  0.000000 0.000000 0.300000 2.000000 3.000000 0.800000",)
    assert rows.nbfix == ("RAH2  RBO3  2.880800 0.050500 0.000000 0.000000 0.000000 0.000000",)
    assert correction_flags(_ff("x")) == {"has_ghbfix": True, "has_nbfix": True}


def test_all_zero_table_and_missing_file():
    zeros = _ff("FF14SB", ";header\nRAC2  RAC2  0 0 0 0 0 0\n")
    assert correction_flags(zeros) == {"has_ghbfix": False, "has_nbfix": False}
    assert correction_flags({"ff_name": "none"}) == {"has_ghbfix": False, "has_nbfix": False}


def test_mdin_ghbfix_block():
    mdin = render_amber_mdin(AmberMdinRequest(duration_ns=1.0, dt_fs=2.0, ghbfix=True))
    assert "  nmropt=1,\n/\n&wt\n" in mdin
    assert "type='REST', istep1=1, istep2=500000, value1=1.0, value2=1.0," in mdin
    assert mdin.endswith("&wt\n  type='END',\n/\nLISTOUT=ghbfix.out\nDISANG=ghbfix.f\n")
    plain = render_amber_mdin(AmberMdinRequest(duration_ns=1.0, dt_fs=2.0))
    assert "nmropt" not in plain and "DISANG" not in plain
