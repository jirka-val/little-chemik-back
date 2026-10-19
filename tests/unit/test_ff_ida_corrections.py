"""
Šestý FF soubor (gHBfix/NBfix) načtený do FF_IDA.ff - vstup pro promítnutí
NBfix do prmtop v create_AMBER_topology.
"""

import pytest

from app.utils.adams4sims_processing_library import FF_IDA

pytestmark = pytest.mark.unit

_RTP = """[ RA ]
 [ atoms ]
  H2  RAH2  0.0  1
  O3  RBO3  0.0  2
"""
_NB = """[ atomtypes ]
RAH2  1  1.008  0.0  A  2.42146e-01  6.27600e-02
RBO3  8  16.00  0.0  A  3.00001e-01  7.11280e-01
"""
_B = "[ bondtypes ]\n"
_CORRECTIONS = """;Type1 Type2   R       eps      eta      rbeg     rend     c
RAH2  RAH2  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
RBO3  RAH2  2.880800 0.050500 0.000000 0.000000 0.000000 0.000000
RAH6  RAN1  0.000000 0.000000 0.300000 2.000000 3.000000 0.800000
"""


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def _ff(tmp_path, corrections=None):
    return FF_IDA.ff(
        _write(tmp_path, "ff.rtp", _RTP),
        _write(tmp_path, "nb.itp", _NB),
        _write(tmp_path, "b.itp", _B),
        None,
        _write(tmp_path, "corr.itp", corrections) if corrections is not None else None,
    )


def test_corrections_split_into_nbfix_and_ghbfix(tmp_path):
    ff = _ff(tmp_path, _CORRECTIONS)
    # klíč je seřazená dvojice typů, stejně jako u bondtypes
    assert ff.nbfix == {("RAH2", "RBO3"): [2.8808, 0.0505]}
    assert ff.ghbfix == {("RAH6", "RAN1"): [0.3, 2.0, 3.0, 0.8]}


def test_without_corrections_file(tmp_path):
    ff = _ff(tmp_path)
    assert ff.nbfix == {} and ff.ghbfix == {}
    assert ff.units["RA"]["atoms"]["type"] == ["RAH2", "RBO3"]


def test_malformed_row_is_rejected(tmp_path):
    with pytest.raises(Exception, match="ghbfix_nbfix"):
        _ff(tmp_path, "RAH2 RBO3 2.88 0.05\n")
