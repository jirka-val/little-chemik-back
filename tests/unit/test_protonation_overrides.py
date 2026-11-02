"""
Unit testy pro ručně vnucené protonační stavy (forced_states) v řešiči
builderu a pro serializaci HIS reportu (titratable_residue_report).

Řešič se volá přímo nad rodinou HIS z data/protonation_states_v1.json a
syntetickou evidencí vodíkové vazby - bez parseru, šablon a silových polí.
"""

import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

BACKEND_DIR = Path(__file__).resolve().parents[2]
BUILDER_DIR = BACKEND_DIR / "app" / "builder"
if str(BUILDER_DIR) not in sys.path:
    sys.path.insert(0, str(BUILDER_DIR))

from forge_molecule_parser import Residue  # noqa: E402
from forge_molecule_state_assignment import (  # noqa: E402
    FixedSiteEvidence,
    ProtonationAssignmentReport,
    ProtonationResidueAssignment,
    _load_protonation_families,
    _solve_hydrogen_bond_network_states,
    _StateVariable,
)

from app.services.structure.reports import titratable_residue_report  # noqa: E402

STATE_DATA = json.loads((BACKEND_DIR / "data" / "protonation_states_v1.json").read_text(encoding="utf-8"))


def _his(resname="HIS"):
    families, by_state = _load_protonation_families(STATE_DATA, 7.4)
    family = by_state[("P", "HIE")] if resname == "HIS" else by_state[("P", resname)]
    residue = Residue(chain_id="A", resseq=24, icode="", ff_resname="HIE", atoms={}, index_in_chain=0)
    item = _StateVariable(index=0, residue=residue, kind="protonation", options=family.options,
                          variable_sites=family.variable_sites, default_option=family.default_option,
                          default_level_index=family.default_level_index)
    return item, family


def _option(family, resname):
    return next(o for o in family.options if o.resname == resname)


def _nd1_must_donate():
    # Partner (např. ASP OD1) je jednoznačný akceptor -> ND1 musí mít vodík.
    # Při pH 7.4 je výchozí HIE, takže evidence musí řešič přetlačit na HID.
    return {
        ("A", 24, "", "ND1"): FixedSiteEvidence(
            variable_site=("A", 24, "", "ND1"),
            required_role="donor",
            partner_site=("A", 45, "", "OD1"),
            distance_angstrom=2.8,
            donor_deviation_deg=10.0,
        )
    }


def test_solver_follows_hbond_evidence_without_override():
    item, family = _his()
    report = ProtonationAssignmentReport(pH=7.4)
    chosen = _solve_hydrogen_bond_network_states([item], _nd1_must_donate(), [], report)
    assert chosen[0].resname == "HID"
    assert report.conflicts == []


def test_forced_state_wins_and_contradiction_is_reported():
    item, family = _his()
    report = ProtonationAssignmentReport(pH=7.4)
    chosen = _solve_hydrogen_bond_network_states(
        [item], _nd1_must_donate(), [], report, forced={0: _option(family, "HIE")}
    )
    assert chosen[0].resname == "HIE"
    # HIE nemá vodík na ND1 - uživatel volí proti geometrii, musí to být vidět.
    assert [c.kind for c in report.conflicts] == ["unsatisfied_fixed_evidence"]


def test_report_serializes_source_evidence_and_base_names():
    prot = ProtonationAssignmentReport(pH=7.4)
    prot.family_defaults.append((("NHIP", "NHIE", "NHID"), "NHIE"))
    prot.fixed_evidence.extend(_nd1_must_donate().values())
    prot.assignments.append(ProtonationResidueAssignment(
        residue_key=("A", 24, ""), old_mol_type="P", old_resname="NHIS",
        new_mol_type="P", new_resname="NHID", default_resname="NHIE",
        is_default=False, is_forced=True,
    ))

    [row] = titratable_residue_report(prot)

    assert row["residue"] == "A24"
    assert row["assigned"] == "HID"
    assert row["default"] == "HIE"
    assert row["options"] == ["HIP", "HIE", "HID"]
    assert row["source"] == "user"
    assert row["evidence"] == [{
        "kind": "fixed", "site": "ND1", "role": "donor",
        "partner": "A45:OD1", "distance_angstrom": 2.8,
    }]


def _his_assignment(resseq, forced=False):
    return ProtonationResidueAssignment(
        residue_key=("A", resseq, ""), old_mol_type="P", old_resname="HIS",
        new_mol_type="P", new_resname="HIE", default_resname="HIE",
        is_default=True, is_forced=forced,
    )


def test_review_flags_only_uncertain_residues():
    from forge_molecule_state_assignment import ProtonationConflict
    from app.services.structure.reports import build_histidine_review

    prot = ProtonationAssignmentReport(pH=7.4)
    prot.family_defaults.append((("HIP", "HIE", "HID"), "HIE"))
    prot.fixed_evidence.extend(_nd1_must_donate().values())      # A24: jednoznačná H-vazba
    prot.assignments.append(_his_assignment(24))
    prot.assignments.append(_his_assignment(57))                 # A57: žádná H-vazba
    prot.assignments.append(_his_assignment(90, forced=True))    # A90: rozhodl uživatel
    prot.conflicts.append(ProtonationConflict(kind="x", message="general problem"))

    review = build_histidine_review(prot)

    assert [r["residue"] for r in review["histidines"]] == ["A57"]
    assert review["histidines"][0]["review_reasons"] == ["no_hbond"]
    assert review["total_histidines"] == 3
    assert review["general_issues"] == ["general problem"]


def test_review_is_empty_when_everything_is_decided():
    from app.services.structure.reports import build_histidine_review

    prot = ProtonationAssignmentReport(pH=7.4)
    prot.family_defaults.append((("HIP", "HIE", "HID"), "HIE"))
    prot.fixed_evidence.extend(_nd1_must_donate().values())
    prot.assignments.append(_his_assignment(24))
    assert build_histidine_review(prot)["histidines"] == []


def test_ion_placement_keeps_polymer_numbering():
    # Odebrání vod pro ionty dřív přečíslovalo VŠECHNY řetězce od 1 (D24 -> D5),
    # takže ruční volby i čísla reziduí po přípravě přestaly sedět.
    import numpy as np
    from forge_molecule_ions import _WaterCandidate, _remove_selected_waters
    from forge_molecule_parser import Chain, Molecule

    def res(chain, resseq, name, group):
        return Residue(chain_id=chain, resseq=resseq, icode="", ff_resname=name,
                       atoms={}, index_in_chain=0, group=group)

    protein = [res("D", n, "ALA", "P") for n in (20, 21, 24)]
    waters = [res("W", n, "WAT", "W3") for n in (1, 2, 3)]
    mol = Molecule(chains={"D": Chain("D", protein), "W": Chain("W", waters)})

    _remove_selected_waters(mol, [_WaterCandidate("W", waters[0], np.zeros(3), active=False)])

    assert [r.resseq for r in mol.chains["D"].residues] == [20, 21, 24]
    assert [r.resseq for r in mol.chains["W"].residues] == [1, 2]


# ---------------------------------------------------------------------------
# Orientace ASN/GLN v síti vodíkových vazeb (orientation_families) a jejich
# vynucení uživatelem (forced_orientations)
# ---------------------------------------------------------------------------

def _ubiquitin_states(**kwargs):
    from app.services.analysis import build_sequence_tokens
    from app.services.structure import forge_service
    from forge_molecule_parser import build_molecule_from_forge_json
    from forge_molecule_state_assignment import assign_molecule_states

    resources = forge_service._static_resources()
    pdb = (BACKEND_DIR / "tests" / "fixtures" / "pdb" / "golden" / "HIS_1UBQ_ubiquitin_his68.pdb").read_text()
    pdb = forge_service._strip_unrecognized_heterogens(pdb, "remove_all")
    molecule = build_molecule_from_forge_json(
        {"pdb_text": pdb, "missing_atoms": build_sequence_tokens(pdb, chain=None, fill_gaps=True)},
        resources.converting_dictionary,
    )
    return assign_molecule_states(
        molecule, resources.converting_dictionary, resources.building_template, resources.state_definitions,
        pH=7.4, **kwargs,
    )


def _orientation(report, resseq):
    return next(a for a in report.orientation.assignments if a.residue_key[1] == resseq)


def test_solver_flips_gln49_of_ubiquitin():
    molecule, report = _ubiquitin_states()
    gln49 = _orientation(report, 49)
    assert gln49.flipped and gln49.evaluated and not gln49.is_forced
    assert all(a.evaluated for a in report.orientation.assignments)


def test_forced_orientation_keeps_the_input_amide():
    _molecule, free = _ubiquitin_states()
    _molecule, report = _ubiquitin_states(forced_orientations={("A", 49, ""): "as_is"})
    gln49 = _orientation(report, 49)
    assert not gln49.flipped and gln49.is_forced
    # ostatní amidy rozhoduje solver dál stejně
    others = {a.residue_key: a.flipped for a in free.orientation.assignments if a.residue_key[1] != 49}
    assert {a.residue_key: a.flipped for a in report.orientation.assignments if a.residue_key[1] != 49} == others


def test_forced_orientation_and_state_that_do_not_apply_are_reported():
    _molecule, report = _ubiquitin_states(
        forced_orientations={("A", 1, ""): "flipped"},       # MET1 nemá amid
        forced_protonation_states={("A", 49, ""): "HIE"},    # GLN49 nemá protonační varianty
    )
    assert [c.kind for c in report.orientation.conflicts] == ["forced_orientation_not_applicable"]
    assert "forced_state_not_applicable" in [c.kind for c in report.protonation.conflicts]


def test_metal_coordinated_histidines_keep_the_lone_pair_for_zinc():
    """1CA2: HIS94/96 koordinují zinek přes NE2 (HID), HIS119 přes ND1 (HIE);
    GLN92 dělá H-vazbu na ND1-H HIS94 a nemá se otáčet."""
    from app.services.analysis import build_sequence_tokens
    from app.services.structure import forge_service
    from forge_molecule_parser import build_molecule_from_forge_json
    from forge_molecule_state_assignment import assign_molecule_states

    resources = forge_service._static_resources()
    pdb = (BACKEND_DIR / "tests" / "fixtures" / "pdb" / "golden" / "HIS_1CA2_carbonic_anhydrase_zinc.pdb").read_text()
    pdb = forge_service._strip_unrecognized_heterogens(pdb, "remove_all")
    molecule = build_molecule_from_forge_json(
        {"pdb_text": pdb, "missing_atoms": build_sequence_tokens(pdb, chain=None, fill_gaps=True)},
        resources.converting_dictionary,
    )
    _molecule, report = assign_molecule_states(
        molecule, resources.converting_dictionary, resources.building_template, resources.state_definitions,
        pH=7.4,
    )
    states = {a.residue_key[1]: a.new_resname for a in report.protonation.assignments}
    assert (states[94], states[96], states[119]) == ("HID", "HID", "HIE")
    assert not _orientation(report, 92).flipped
    assert report.protonation.conflicts == []
    metal = [e for e in report.protonation.fixed_evidence if e.partner_site[3] == "ZN"]
    assert {e.variable_site[1] for e in metal} == {94, 96, 119}
