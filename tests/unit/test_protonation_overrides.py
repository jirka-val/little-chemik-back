"""
Unit testy pro ručně vnucené protonační stavy (forced_states) v řešiči
builderu a pro serializaci HIS reportu (_titratable_residue_report).

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
    _solve_protonation_states,
    _TitratableResidue,
)

from app.services.structure.forge_service import _titratable_residue_report  # noqa: E402

STATE_DATA = json.loads((BACKEND_DIR / "data" / "protonation_states_v1.json").read_text(encoding="utf-8"))


def _his(resname="HIS"):
    families, by_state = _load_protonation_families(STATE_DATA, 7.4)
    family = by_state[("P", "HIE")] if resname == "HIS" else by_state[("P", resname)]
    residue = Residue(chain_id="A", resseq=24, icode="", ff_resname="HIE", atoms={}, index_in_chain=0)
    return _TitratableResidue(index=0, residue=residue, family=family), family


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
    chosen = _solve_protonation_states([item], _nd1_must_donate(), [], report)
    assert chosen[0].resname == "HID"
    assert report.conflicts == []


def test_forced_state_wins_and_contradiction_is_reported():
    item, family = _his()
    report = ProtonationAssignmentReport(pH=7.4)
    chosen = _solve_protonation_states(
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

    [row] = _titratable_residue_report(prot)

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
    from app.services.structure.forge_service import build_protonation_review

    prot = ProtonationAssignmentReport(pH=7.4)
    prot.family_defaults.append((("HIP", "HIE", "HID"), "HIE"))
    prot.fixed_evidence.extend(_nd1_must_donate().values())      # A24: jednoznačná H-vazba
    prot.assignments.append(_his_assignment(24))
    prot.assignments.append(_his_assignment(57))                 # A57: žádná H-vazba
    prot.assignments.append(_his_assignment(90, forced=True))    # A90: rozhodl uživatel
    prot.conflicts.append(ProtonationConflict(kind="x", message="general problem"))

    review = build_protonation_review(prot)

    assert [r["residue"] for r in review["residues"]] == ["A57"]
    assert review["residues"][0]["review_reasons"] == ["no_hbond"]
    assert review["total_titratable"] == 3
    assert review["general_issues"] == ["general problem"]


def test_review_is_none_when_everything_is_decided():
    from app.services.structure.forge_service import build_protonation_review

    prot = ProtonationAssignmentReport(pH=7.4)
    prot.family_defaults.append((("HIP", "HIE", "HID"), "HIE"))
    prot.fixed_evidence.extend(_nd1_must_donate().values())
    prot.assignments.append(_his_assignment(24))
    assert build_protonation_review(prot) is None


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
