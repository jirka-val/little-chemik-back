"""
Unit testy pro kontroly kroku 3.5 "Structure Check" (structure_review.py a
find_removed_heterogens): nulová obsazenost, otočení amidu ASN/GLN,
odstraněné ligandy a aplikace rozhodnutí uživatele na PDB text.
"""

import pytest

from app.services.structure.forge_service import find_removed_heterogens
from app.services.structure.structure_review import (
    apply_structure_edits,
    find_amide_flips,
    find_zero_occupancy,
)

pytestmark = pytest.mark.unit


def _atom(serial, name, resname, resseq, x, y, z, occ=1.0, element=None, record="ATOM", chain="A"):
    element = element or name[0]
    atom_field = name.ljust(4) if len(name) == 4 else f" {name}".ljust(4)
    return (
        f"{record:<6}{serial:>5} {atom_field} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


# ASN10: OD1 na (0,0,0), ND2 na (2,0,0). Backbone O rezidua 50 je 2.9 Å od
# ND2 (dobře: N-H...O) a backbone N rezidua 60 je 2.9 Å od OD1 (dobře).
# Otočení to celé pokazí -> žádný návrh.
def _asn(o_partner_x, n_partner_x):
    return "\n".join([
        _atom(1, "N", "ASN", 10, -3.0, 2.0, 0.0),
        _atom(2, "CA", "ASN", 10, -2.0, 1.5, 0.0),
        _atom(3, "C", "ASN", 10, -2.0, 3.0, 0.0),
        _atom(4, "O", "ASN", 10, -1.0, 3.5, 0.0),
        _atom(5, "CB", "ASN", 10, -1.0, 0.5, 0.0),
        _atom(6, "CG", "ASN", 10, 1.0, 1.0, 0.0),
        _atom(7, "OD1", "ASN", 10, 0.0, 0.0, 0.0),
        _atom(8, "ND2", "ASN", 10, 2.0, 0.0, 0.0),
        _atom(9, "HD21", "ASN", 10, 2.5, -0.8, 0.0, element="H"),
        _atom(10, "O", "GLY", 50, n_partner_x, -2.9, 0.0),
        _atom(11, "N", "GLY", 60, o_partner_x, -2.9, 0.0),
    ])


def test_amide_in_good_orientation_is_not_flagged():
    assert find_amide_flips(_asn(o_partner_x=0.0, n_partner_x=2.0)) == []


def test_amide_in_wrong_orientation_is_suggested_for_flip():
    # Partneři prohození: backbone O míří na OD1 a backbone N na ND2.
    [flip] = find_amide_flips(_asn(o_partner_x=2.0, n_partner_x=0.0))
    assert flip["residue"] == "A10"
    assert flip["score_flipped"] > flip["score_current"]


def test_decided_amide_is_skipped():
    pdb = _asn(o_partner_x=2.0, n_partner_x=0.0)
    assert find_amide_flips(pdb, decided=[("A", 10, "")]) == []


def test_flip_swaps_amide_atoms_and_drops_their_hydrogens():
    edited = apply_structure_edits(_asn(0.0, 2.0), flip_residues=[("A", 10, "")])
    lines = {line[12:16].strip(): line for line in edited.splitlines() if " ASN " in line}
    assert "HD21" not in lines
    assert lines["ND2"][30:38].strip() == "0.000"   # ND2 je teď na místě původního OD1
    assert lines["OD1"][30:38].strip() == "2.000"
    assert lines["ND2"][76:78].strip() == "N"


def test_zero_occupancy_side_chain_is_rebuildable():
    pdb = "\n".join([
        _atom(1, "N", "LYS", 222, 0, 0, 0),
        _atom(2, "CA", "LYS", 222, 1.5, 0, 0),
        _atom(3, "CD", "LYS", 222, 3, 0, 0, occ=0.0),
        _atom(4, "NZ", "LYS", 222, 4, 0, 0, occ=0.0),
        _atom(5, "O", "GLY", 223, 6, 0, 0, occ=0.0),
    ])
    rows = {r["residue"]: r for r in find_zero_occupancy(pdb)}
    assert rows["A222"]["atoms"] == ["CD", "NZ"] and rows["A222"]["rebuildable"]
    assert rows["A223"]["rebuildable"] is False   # páteř nejde postavit znovu

    edited = apply_structure_edits(pdb, rebuild_residues=[("A", 222, "")])
    names = [line[12:16].strip() for line in edited.splitlines()]
    assert names == ["N", "CA", "O"]


def test_removed_heterogens_skip_water_and_ions():
    pdb = "\n".join([
        _atom(1, "C1", "BEN", 1, 0, 0, 0, record="HETATM"),
        _atom(2, "C2", "BEN", 1, 1, 0, 0, record="HETATM"),
        _atom(3, "CA", "CA", 480, 5, 0, 0, record="HETATM", element="CA"),
        _atom(4, "O", "HOH", 500, 9, 0, 0, record="HETATM"),
    ])
    assert find_removed_heterogens(pdb, "keep_water") == [
        {"key": "A:1::BEN", "residue": "A1", "resname": "BEN", "atoms": 2}
    ]
    assert find_removed_heterogens(pdb, "keep_all") == []
