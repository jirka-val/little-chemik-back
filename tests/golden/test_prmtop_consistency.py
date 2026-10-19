"""
Vygenerovaný prmtop musí projít kontrolou délek sekcí z knihovny
(AMBER_topology.check_AMBER_topology) - POINTERS počítají vazby/úhly/
dihedrály, ne čísla, a RADII/SCREEN mají hodnotu i pro ionty. Zároveň ho
musí umět přečíst souhrn nahraných vstupů (app/services/input_files.py).
"""

from __future__ import annotations

import pytest

from app.services.input_files import summarize_topology
from app.utils.adams4sims_processing_library.utils import AMBER_topology
from app.workspaces.manager import workspace_manager

from .pipeline import run_pipeline
from .test_golden_pipeline import CASES, frozen_force_fields  # noqa: F401 (fixture)

pytestmark = pytest.mark.golden

_SOLVATED_PROTEIN = next(c for c in CASES if c.name == "1crn_protein")


def test_generated_prmtop_is_consistent(client, frozen_force_fields):  # noqa: F811
    result = run_pipeline(client, _SOLVATED_PROTEIN)
    path = workspace_manager.get_file_path(result.workspace_id, "structure.prmtop")

    top = AMBER_topology.read_AMBER_topology(str(path))
    assert top, "read_AMBER_topology rejected the generated prmtop (see its printed error)"
    assert top["POINTERS"][2] * 3 == len(top["BONDS_INC_HYDROGEN"])
    assert len(top["RADII"]) == top["POINTERS"][0]

    summary = summarize_topology(path)
    assert summary["ok"]
    assert summary["composition"]["protein"] == 46
    assert summary["composition"]["water"] > 0
    assert set(summary["composition"]["ions"]) == {"Na+", "Cl-"}
    assert summary["total_charge"] == 0.0
    # box z CRYST1 truncated octahedronu, ne 90 stupňů
    assert top["POINTERS"][27] == 2
    assert summary["box"]["beta"] == pytest.approx(109.4712190)
    assert summary["box"]["shape"] == "truncated octahedron"
