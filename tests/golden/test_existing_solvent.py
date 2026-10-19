"""
Nahraný hotový (solvatovaný) systém: rozpozná se, s volbou "keep" projde
přípravou i topologií beze změny vody, iontů a boxu, bez volby a s
ponecháním krystalové vody se příprava odmítne (jinak by se box obalil
dalším boxem).
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from app.services.input_files import summarize_topology
from app.workspaces.manager import workspace_manager

from .pipeline import run_pipeline
from .test_golden_pipeline import CASES, frozen_force_fields  # noqa: F401 (fixture)

pytestmark = pytest.mark.golden

_SOLVATED_PROTEIN = next(c for c in CASES if c.name == "1crn_protein")


def _counts(pdb_text: str):
    waters = {(l[21], l[22:27]) for l in pdb_text.splitlines() if l[17:20] in ("WAT", "HOH")}
    ions = sum(1 for l in pdb_text.splitlines() if l.startswith("HETATM") and l[17:20].strip() in ("NA", "CL"))
    cryst1 = next((l for l in pdb_text.splitlines() if l.startswith("CRYST1")), None)
    return len(waters), ions, cryst1


def test_uploaded_solvated_system(client, frozen_force_fields, tmp_path):  # noqa: F811
    first = run_pipeline(client, _SOLVATED_PROTEIN)
    solvated = tmp_path / "solvated.pdb"
    solvated.write_text(first.files["structure.pdb"], encoding="utf-8")
    before = _counts(solvated.read_text())

    kept = run_pipeline(client, replace(_SOLVATED_PROTEIN, pdb_path=solvated,
                                        prepare_overrides={"existing_solvent": "keep"}))
    detection = client.get(f"/api/analysis/solvation/{kept.workspace_id}")
    assert detection.status_code == 200
    assert _counts(kept.files["structure.pdb"]) == before

    top = summarize_topology(workspace_manager.get_file_path(kept.workspace_id, "structure.prmtop"))
    assert top["total_charge"] == 0.0
    assert top["box"]["shape"] == "truncated octahedron"
    assert kept.responses["topology"]["total_charge"] == 0.0

    ws = client.post("/api/molecules/upload", files={"file": ("s.pdb", solvated.read_bytes())}).json()["workspace_id"]
    analysis = client.get(f"/api/analysis/solvation/{ws}").json()
    assert analysis["solvated"] is True
    assert all(c["status"] == "ok" for c in analysis["checks"]), analysis["checks"]

    with pytest.raises(AssertionError, match="solvated_system_choice_required"):
        run_pipeline(client, replace(_SOLVATED_PROTEIN, pdb_path=solvated, crystal_water_mode="keep_water"))
