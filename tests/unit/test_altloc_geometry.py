"""
Unit testy pro geometrické doporučení AltLoc variant a detekci zlomů
(analysis_service.analyze_pdb_altlocs / find_altloc_selection_breaks).

Regrese z 2TRA: doporučení vybíralo podle occupancy/B-faktoru bez ohledu na
to, jestli zvolená konformace navazuje na sousedy (G34-C36 = B při 50/50
occupancy -> zbytečný zlom U33|G34), a varování ve frontendu porovnávalo jen
písmena. Pravidlo: nejdřív co nejméně zlomů (O3'-P / C-N), pak occupancy a
B-faktor sečtené přes celý kus reziduí propojených altloc atomy.

Stačí syntetické řetězce jen s vazebnými atomy O3' a P - geometrie se
posuzuje jen přes ně.
"""

import pytest

from app.services.analysis_service import analyze_pdb_altlocs, find_altloc_selection_breaks

pytestmark = pytest.mark.unit


def _atom(serial, name, altloc, resname, resseq, x, occ=1.0, b=20.0):
    return (
        f"ATOM  {serial:5d} {name:<4s}{altloc:1s}{resname:>3s} A{resseq:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{b:6.2f}           {name[0]}"
    )


def _pdb(atoms):
    return "\n".join(_atom(i + 1, *a) for i, a in enumerate(atoms)) + "\n"


def _recommended(pdb_text):
    return {r["resseq"]: r["recommended_alt"] for r in analyze_pdb_altlocs(pdb_text)["residues"]}


def test_geometry_beats_higher_occupancy():
    # G2: B má vyšší occupancy, ale jeho P je od O3' U1 3.4 A daleko.
    pdb = _pdb([
        ("O3'", " ", "U", 1, 0.0),
        ("P", "A", "G", 2, 1.6, 0.4),
        ("P", "B", "G", 2, 3.4, 0.6),
        ("O3'", " ", "G", 2, 6.0),
        ("P", " ", "C", 3, 7.6),
    ])
    assert _recommended(pdb) == {2: "A"}


def test_occupancy_decides_when_geometry_is_equal():
    pdb = _pdb([
        ("O3'", " ", "U", 1, 0.0),
        ("P", " ", "G", 2, 1.6),
        ("C1'", "A", "G", 2, 3.0, 0.3),
        ("C1'", "B", "G", 2, 3.2, 0.7),
    ])
    assert _recommended(pdb) == {2: "B"}


def test_connected_piece_is_decided_as_a_whole():
    # G2 a C3 jsou propojené altloc atomy O3'(2)-P(3): A-A i B-B navazují,
    # křížové kombinace ne. Po reziduích by vyšlo A (0.6) + B (0.6) = zlom;
    # přes celý kus je occupancy shodná a rozhodne nižší B-faktor větve B.
    pdb = _pdb([
        ("O3'", "A", "G", 2, 0.0, 0.6, 30.0),
        ("O3'", "B", "G", 2, 5.0, 0.4, 20.0),
        ("P", "A", "C", 3, 1.6, 0.4, 30.0),
        ("P", "B", "C", 3, 6.6, 0.6, 20.0),
    ])
    assert _recommended(pdb) == {2: "B", 3: "B"}
    assert find_altloc_selection_breaks(pdb, {"A_2_G": "B", "A_3_C": "B"}) == []


def test_selection_breaks_report_only_avoidable_breaks():
    pdb = _pdb([
        ("O3'", " ", "U", 1, 0.0),
        ("P", "A", "G", 2, 1.6, 0.5),
        ("P", "B", "G", 2, 3.4, 0.5),
        ("O3'", " ", "G", 2, 6.0),
        # C3 je od G2 daleko při jakékoli volbě - skutečná mezera, nehlásí se.
        ("P", "A", "C", 3, 12.0, 0.5),
        ("P", "B", "C", 3, 13.0, 0.5),
    ])
    breaks = find_altloc_selection_breaks(pdb, {"A_2_G": "B", "A_3_C": "A"})
    assert len(breaks) == 1
    assert breaks[0]["prev"]["resseq"] == 1 and breaks[0]["prev"]["altloc"] is None
    assert breaks[0]["next"]["resseq"] == 2 and breaks[0]["next"]["altloc"] == "B"
    assert find_altloc_selection_breaks(pdb, {"A_2_G": "A", "A_3_C": "A"}) == []


def test_occupancy_and_bfactor_are_averaged_over_variant_atoms():
    pdb = _pdb([
        ("C1'", "A", "G", 2, 3.0, 0.2, 10.0),
        ("C2'", "A", "G", 2, 4.0, 0.4, 30.0),
        ("C1'", "B", "G", 2, 3.2, 0.8, 50.0),
        ("C2'", "B", "G", 2, 4.2, 0.6, 70.0),
    ])
    alt = analyze_pdb_altlocs(pdb)["residues"][0]["altLocs"]
    assert alt["A"] == {"occupancy": 30.0, "bFactor": 20.0}
    assert alt["B"] == {"occupancy": 70.0, "bFactor": 60.0}
