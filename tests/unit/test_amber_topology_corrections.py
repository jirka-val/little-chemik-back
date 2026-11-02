"""
HMR, NBfix, gHBfix a 1-4 interakce v AMBER_topology.create_AMBER_topology
nad skutečným RNA FF s korekcemi (OL3CP_gHBfix21_NBfix-0BPh z data/ff_cache).
Molekula je krátký řetězec RG5-RA-RU-RC3 bez souřadnic - topologie je
potřebuje jen pro box.
"""

import math
import re
from pathlib import Path

import pytest

from app.utils.adams4sims_processing_library import FF_IDA
from app.utils.adams4sims_processing_library.utils import AMBER_topology

FF_NAME = "OL3CP_gHBfix21_NBfix-0BPh"
FF_DIR = Path(__file__).resolve().parents[2] / "data" / "ff_cache" / FF_NAME

pytestmark = [
    pytest.mark.unit,
    pytest.mark.skipif(not FF_DIR.exists(), reason=f"data/ff_cache/{FF_NAME} not present locally"),
]


@pytest.fixture(scope="module")
def rna_ff():
    path = FF_DIR
    return FF_IDA.ff(
        str(path / f"{FF_NAME}.rtp"),
        str(path / f"nonbonded_{FF_NAME}.itp"),
        str(path / f"bonded_{FF_NAME}.itp"),
        str(path / f"{FF_NAME}.atp"),
        str(path / f"ghbfix_nbfix_{FF_NAME}.itp"),
    )


def _mol(ff):
    residues = [{"resn": resn, "mol_type": "R"} for resn in ("RG5", "RA", "RU", "RC3")]
    return {"residues": residues, "force_field_data": {"R": ff}}


@pytest.fixture
def build(rna_ff, tmp_path):
    def _build(**kwargs):
        ghbfix = tmp_path / "ghbfix.f"
        top = AMBER_topology.create_AMBER_topology(_mol(rna_ff), ghbfix_file=ghbfix, **kwargs)
        return top, ghbfix
    return _build


def _pair_lj(top, i, j):
    """Rmin (A) a eps (kcal/mol) páru atomů i, j (0-based) z A/B koeficientů."""
    ntypes = top["POINTERS"][1]
    ti, tj = top["ATOM_TYPE_INDEX"][i] - 1, top["ATOM_TYPE_INDEX"][j] - 1
    index = top["NONBONDED_PARM_INDEX"][ti * ntypes + tj] - 1
    a, b = top["LENNARD_JONES_ACOEF"][index], top["LENNARD_JONES_BCOEF"][index]
    return (2 * a / b) ** (1 / 6), b * b / (4 * a)


def _atom(top, residue, name):
    start = top["RESIDUE_POINTER"][residue] - 1
    end = top["RESIDUE_POINTER"][residue + 1] - 1 if residue + 1 < len(top["RESIDUE_POINTER"]) else len(top["ATOM_NAME"])
    return start + top["ATOM_NAME"][start:end].index(name)


def test_hmr_sets_hydrogens_to_3_024_da_and_keeps_total_mass(build):
    plain, _ = build()
    hmr, _ = build(HMR=True)
    hydrogens = [i for i, z in enumerate(hmr["ATOMIC_NUMBER"]) if z == 1]
    assert hydrogens and all(hmr["MASS"][i] == pytest.approx(3.024) for i in hydrogens)
    assert math.fsum(hmr["MASS"]) == pytest.approx(math.fsum(plain["MASS"]))
    heavy = [i for i, z in enumerate(hmr["ATOMIC_NUMBER"]) if z > 1]
    assert any(hmr["MASS"][i] < plain["MASS"][i] for i in heavy)


def test_ghbfix_file_has_two_rst_terms_per_pair(build):
    _, ghbfix = build()
    text = ghbfix.read_text(encoding="utf-8")
    pairs = re.findall(r"&rst iat=(\d+), (\d+),", text)
    assert pairs and len(pairs) % 2 == 0
    assert all(int(i) < int(j) for i, j in pairs)
    assert "r1=0.0, r2=2, r3=2, r4=3," in text
    assert "r1=0.0, r2=2, r3=2, r4=2.8," in text


def test_ghbfix_skips_base_pairs_within_one_residue(build, rna_ff):
    top, ghbfix = build()
    residue_of = {}
    for residue, start in enumerate(top["RESIDUE_POINTER"]):
        end = top["RESIDUE_POINTER"][residue + 1] - 1 if residue + 1 < len(top["RESIDUE_POINTER"]) else len(top["ATOM_NAME"])
        residue_of.update({atom: residue for atom in range(start - 1, end)})
    base_type = re.compile(r"^R[GACU][NOH]\d$")
    pairs = {(int(i) - 1, int(j) - 1) for i, j in re.findall(r"&rst iat=(\d+), (\d+),", ghbfix.read_text(encoding="utf-8"))}
    types = [atom_type for res in _mol(rna_ff)["residues"] for atom_type in rna_ff.units[res["resn"]]["atoms"]["type"]]
    intra_base = [(i, j) for i, j in pairs if residue_of[i] == residue_of[j]
                  and base_type.match(types[i]) and base_type.match(types[j])]
    assert intra_base == []


def test_end_group_14_pairs_are_counted_once(build):
    top, _ = build()
    bonded = {}
    for key in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        bonds = top[key]
        for k in range(0, len(bonds), 3):
            a, b = bonds[k] // 3, bonds[k + 1] // 3
            bonded.setdefault(a, set()).add(b)
            bonded.setdefault(b, set()).add(a)
    counted = []
    for key in ("DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"):
        dihedrals = top[key]
        for k in range(0, len(dihedrals), 5):
            a, _, c, d = (dihedrals[k + n] for n in range(4))
            if c >= 0 and d >= 0:
                a, d = a // 3, abs(d) // 3
                counted.append(tuple(sorted((a, d))))
    assert len(counted) == len(set(counted))
    # Ribóza je pětičlenný kruh: kraje dihedrálu C1'-C2'-C3'-C4' jsou zároveň
    # 1-3 přes O4', takže 1-4 se pro ně počítat nesmí.
    for a, d in counted:
        assert d not in bonded[a] and bonded[a].isdisjoint(bonded[d])


@pytest.mark.xfail(strict=True, reason="NBfix v IDA je Rmin (A) a eps (kcal/mol), kód ho přepočítává jako nm/kJ - čeká na potvrzení")
def test_nbfix_overrides_base_hydrogen_phosphate_oxygen_pair(build):
    top, _ = build()
    h8 = _atom(top, 0, "H8")       # RGH8 (RG5)
    o5 = _atom(top, 1, "O5'")      # RBO5 (RA)
    rmin, eps = _pair_lj(top, h8, o5)
    assert rmin == pytest.approx(2.8808, abs=1e-4)
    assert eps == pytest.approx(0.0505, abs=1e-4)
