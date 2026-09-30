"""
Unit testy pro detekci identických kopií molekuly v jednom PDB
(analysis_service.find_identical_chain_copies) a jejich odstranění přes
process_structure(remove_chains=...). Vzor: 3SKR - řetězce A a B jsou dvě
kopie téhož riboswitche, REMARK 350 je uvádí jako dvě monomerní jednotky.
"""

import pytest

from app.services.analysis_service import (
    analyze_pdb_altlocs,
    find_identical_chain_copies,
    process_structure,
)

pytestmark = pytest.mark.unit


def _atom(record, serial, name, resname, chain, resseq, b):
    return (
        f"{record:<6s}{serial:5d} {name:<4s} {resname:>3s} {chain}{resseq:4d}    "
        f"{float(resseq):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{b:6.2f}           {name[0]}"
    )


def _copies_pdb(with_seqres=True, remark350=True, b_a=40.0, b_b=60.0, drop_b_residue=False):
    lines = []
    if remark350:
        lines += [
            "REMARK 350 BIOMOLECULE: 1",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A",
            "REMARK 350 BIOMOLECULE: 2",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: B",
        ]
    if with_seqres:
        lines += [
            "SEQRES   1 A    3    G   C   U",
            "SEQRES   1 B    3    G   C   U",
        ]
    serial = 1
    for chain, b in (("A", b_a), ("B", b_b)):
        for resseq, resname in ((1, "G"), (2, "C"), (3, "U")):
            if drop_b_residue and chain == "B" and resseq == 3:
                continue
            lines.append(_atom("ATOM", serial, "P", resname, chain, resseq, b))
            serial += 1
        lines.append(_atom("HETATM", serial, "MG", "MG", chain, 101, 30.0))
        serial += 1
    return "\n".join(lines) + "\n"


def test_detects_copies_and_recommends_lower_bfactor():
    groups = find_identical_chain_copies(_copies_pdb())
    assert len(groups) == 1
    group = groups[0]
    assert [c["chain"] for c in group["chains"]] == ["A", "B"]
    assert group["recommended"] == "A"
    assert group["separateBiologicalUnits"] is True
    assert group["chains"][0]["heteroGroups"] == 1


def test_more_complete_copy_wins_over_bfactor():
    groups = find_identical_chain_copies(_copies_pdb(b_a=60.0, b_b=40.0, drop_b_residue=True))
    assert groups[0]["recommended"] == "A"


def test_detects_copies_without_seqres():
    groups = find_identical_chain_copies(_copies_pdb(with_seqres=False, remark350=False))
    assert groups[0]["recommended"] == "A"
    assert groups[0]["fromSeqres"] is False
    assert groups[0]["separateBiologicalUnits"] is False


def test_single_chain_has_no_copies():
    pdb = "\n".join(l for l in _copies_pdb().splitlines() if l[21:22] != "B" and "SEQRES   1 B" not in l)
    assert find_identical_chain_copies(pdb) == []


def test_analysis_reports_copy_groups():
    assert len(analyze_pdb_altlocs(_copies_pdb())["copyGroups"]) == 1


def test_remove_chains_drops_copy_with_its_hetero_groups():
    out = process_structure(_copies_pdb(), target_model=1, apply_symmetry=False, selection={}, remove_chains=["B"])
    chains = {l[21] for l in out.splitlines() if l.startswith(("ATOM", "HETATM"))}
    assert chains == {"A"}
    assert sum(1 for l in out.splitlines() if l.startswith("HETATM")) == 1


def test_remove_nothing_keeps_all_copies():
    out = process_structure(_copies_pdb(), target_model=1, apply_symmetry=False, selection={})
    chains = {l[21] for l in out.splitlines() if l.startswith(("ATOM", "HETATM"))}
    assert chains == {"A", "B"}
