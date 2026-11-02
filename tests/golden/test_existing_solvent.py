"""
Nahraný hotový (solvatovaný) systém: rozpozná se, s volbou "keep" projde
přípravou i topologií beze změny vody, iontů a boxu, bez volby a s
ponecháním krystalové vody se příprava odmítne (jinak by se box obalil
dalším boxem).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
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


def _strip_cryst1_and_rotate(pdb_text: str) -> str:
    """Jako export z programu, který box nezapíše: bez CRYST1, otočené a posunuté."""
    angle = np.radians(33.0)
    rot = np.array([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("CRYST1"):
            continue
        if line.startswith(("ATOM", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ rot.T + (7.0, -4.0, 11.0)
            line = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def test_solvated_system_without_cryst1_keeps_a_reconstructed_box(client, frozen_force_fields, tmp_path):  # noqa: F811
    first = run_pipeline(client, _SOLVATED_PROTEIN)
    original = first.files["structure.pdb"]
    stripped = tmp_path / "no_cryst1.pdb"
    stripped.write_text(_strip_cryst1_and_rotate(original), encoding="utf-8")
    waters, ions, cryst1 = _counts(original)

    ws = client.post("/api/molecules/upload", files={"file": ("s.pdb", stripped.read_bytes())}).json()["workspace_id"]
    analysis = client.get(f"/api/analysis/solvation/{ws}").json()
    assert analysis["solvated"] is True and analysis["box_source"] == "reconstructed"
    assert analysis["box"]["shape"] == "truncated octahedron"
    assert analysis["box"]["a"] == pytest.approx(float(cryst1[6:15]), abs=0.5)
    assert analysis["checks"][0]["id"] == "box_source"

    kept = run_pipeline(client, replace(_SOLVATED_PROTEIN, pdb_path=stripped,
                                        prepare_overrides={"existing_solvent": "keep"}))
    kept_waters, kept_ions, kept_cryst1 = _counts(kept.files["structure.pdb"])
    assert (kept_waters, kept_ions) == (waters, ions)
    assert kept_cryst1 is not None
    top = summarize_topology(workspace_manager.get_file_path(kept.workspace_id, "structure.prmtop"))
    assert top["box"]["shape"] == "truncated octahedron"
    assert top["total_charge"] == 0.0


def test_removing_the_solvent_right_after_upload(client, frozen_force_fields, tmp_path):  # noqa: F811
    """Okno po načtení, "solvatovat znovu": voda a ionty pryč hned, pak běžná příprava."""
    first = run_pipeline(client, _SOLVATED_PROTEIN)
    solvated = first.files["structure.pdb"]
    waters, _ions, _cryst1 = _counts(solvated)

    ws = client.post("/api/molecules/upload", files={"file": ("s.pdb", solvated.encode())}).json()["workspace_id"]
    response = client.post(f"/api/analysis/solvation/{ws}/remove")
    assert response.status_code == 200
    removed = response.json()["removed"]
    assert removed["waters"] == waters and removed["ions"]
    analysis = client.get(f"/api/analysis/solvation/{ws}").json()
    assert analysis["solvated"] is False and analysis["waters"] == 0

    stripped = tmp_path / "stripped.pdb"
    stripped.write_text(workspace_manager.get_file_path(ws, "structure.pdb").read_text(encoding="utf-8"),
                        encoding="utf-8")
    assert _counts(stripped.read_text()) == (0, 0, None)
    # bez vody se struktura připraví jako každá jiná (default add_solvent)
    again = run_pipeline(client, replace(_SOLVATED_PROTEIN, pdb_path=stripped))
    new_waters, new_ions, new_cryst1 = _counts(again.files["structure.pdb"])
    assert new_waters > 0 and new_ions > 0 and new_cryst1 is not None
