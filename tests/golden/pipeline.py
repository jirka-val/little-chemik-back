"""
Drives the whole preparation pipeline through the HTTP API, the same way the
frontend does in Guided mode, and collects every user-visible output.

Guided mode, step by step (see little-chemik-front/src):
  1. upload                    molecule-loader.ts -> POST /molecules/upload
  2. sequence + altloc scan    GET /analysis/sequence, GET /analysis/altlocs
  3. auto-apply recommended    altloc-panel.ts::autoApplyRecommended -> POST /analysis/clean-altlocs
  4. force fields              GET /forcefields?mode=guided, take every `is_default` entry
  5. validation                POST /validation/check
  6. prepare                   hydrogen-panel.ts -> POST /validation/prepare
  7. topology                  export-modal.ts -> POST /topology/{ws}/generate
  8. export                    POST /download/{ws}/export (pdb + prmtop + crd + mdin as ZIP)

The returned `PipelineResult` holds JSON responses and files keyed by a
stable name, so the golden test can store and compare them.
"""

from __future__ import annotations

import io
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

from fastapi.testclient import TestClient


@dataclass
class PipelineCase:
    name: str
    pdb_path: Path
    add_solvent: bool = True
    crystal_water_mode: str = "remove_all"
    ionic_strength: float = 0.15
    hmr: bool = False
    # Extra PreparationRequest fields (positive_ion, box_shape, ...).
    prepare_overrides: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineResult:
    responses: Dict[str, Any] = field(default_factory=dict)
    files: Dict[str, str] = field(default_factory=dict)


def _ok(response, step: str) -> Any:
    assert response.status_code == 200, f"{step}: HTTP {response.status_code}: {response.text[:2000]}"
    return response.json()


def _recommended_altloc_request(altlocs: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Mirrors altloc-panel.ts: what Guided mode sends to /clean-altlocs, or None if nothing to resolve."""
    models = altlocs.get("models") or []
    has_models = len(models) > 1
    has_symmetry = bool(altlocs.get("hasSymmetry"))
    residues = altlocs.get("residues") or []
    has_altlocs = bool(altlocs.get("hasAltLocs")) and len(residues) > 0
    copy_groups = altlocs.get("copyGroups") or []
    if not (has_models or has_symmetry or has_altlocs or copy_groups):
        return None

    selection = {}
    for res in residues:
        key = f"{res['chain']}_{res['resseq']}_{res['resname']}"
        selection[key] = res.get("recommended_alt") or next(iter(res["altLocs"]))

    remove_chains = []
    for group in copy_groups:
        keep = group.get("recommended")
        if keep and keep != "all":
            remove_chains += [c["chain"] for c in group["chains"] if c["chain"] != keep]

    return {
        "model": models[0] if has_models else 1,
        "apply_symmetry": has_symmetry,
        "selection": selection,
        "remove_chains": remove_chains,
    }


def _guided_ff_selections(ff_response: Dict[str, Any]) -> Dict[str, Any]:
    """Mirrors forcefield-panel.ts confirm handler with the Guided preselection (`is_default`)."""
    selections: Dict[str, Any] = {}
    for group_key in ff_response["required_groups"]:
        if group_key == "W":
            profile = next(p for p in ff_response["water_profiles"] if p.get("is_default"))
            selections["W"] = profile["water_ff"]
            continue
        candidates = ff_response["forcefields_by_group"].get(group_key, [])
        selections[group_key] = next(ff for ff in candidates if ff.get("is_default"))
    return selections


def run_pipeline(client: TestClient, case: PipelineCase) -> PipelineResult:
    result = PipelineResult()
    r = result.responses

    with open(case.pdb_path, "rb") as fh:
        upload = _ok(client.post("/api/molecules/upload", files={"file": (case.pdb_path.name, fh, "chemical/x-pdb")}), "upload")
    ws = upload["workspace_id"]
    r["upload"] = {k: v for k, v in upload.items() if k != "workspace_id"}

    r["sequence_initial"] = _ok(client.get(f"/api/analysis/sequence/{ws}?fill_gaps=true"), "sequence")
    altlocs = _ok(client.get(f"/api/analysis/altlocs/{ws}"), "altlocs")
    r["altlocs"] = altlocs

    clean_request = _recommended_altloc_request(altlocs)
    if clean_request is not None:
        r["clean_altlocs_request"] = clean_request
        r["clean_altlocs"] = _ok(client.post(f"/api/analysis/clean-altlocs/{ws}", json=clean_request), "clean-altlocs")
        r["sequence_after_clean"] = _ok(client.get(f"/api/analysis/sequence/{ws}?fill_gaps=true"), "sequence")

    positive_ion = case.prepare_overrides.get("positive_ion", "Na+")
    negative_ion = case.prepare_overrides.get("negative_ion", "Cl-")
    ff_response = _ok(client.get(
        f"/api/forcefields/{ws}",
        params={"mode": "guided", "positive_ion": positive_ion, "negative_ion": negative_ion},
    ), "forcefields")
    r["forcefields"] = _summarize_ff_response(ff_response)
    ff_selections = _guided_ff_selections(ff_response)
    r["ff_selections"] = {k: v.get("display_name") or v.get("ff_name") for k, v in ff_selections.items()}

    r["validation_check"] = _ok(client.post("/api/validation/check", json={"workspace_id": ws}), "validation/check")

    prepare_request = {
        "workspace_id": ws,
        "ff_selections": ff_selections,
        "ph": 7.4,
        "crystal_water_mode": case.crystal_water_mode,
        "add_solvent": case.add_solvent,
        "box_padding_nm": 1.0,
        "box_shape": "octahedron",  # hydrogen-panel.ts default
        "ionic_strength": case.ionic_strength,
        "positive_ion": positive_ion,
        "negative_ion": negative_ion,
        "clean_crystal_ions": True,
        "replace_structural_multivalent_with_mg": False,
        "concentration_mode": "water_ratio",  # hydrogen-panel.ts default
        "additional_salts": [],
        "protonation_overrides": [],
        "structure_decisions": {"amide_flips": [], "zero_occupancy": [], "acknowledged_heterogens": []},
        "review_structure": False,
        "auto_amide_flips": True,
        **case.prepare_overrides,
    }
    # hydrogen-panel.ts prepares through /sidechains/start; in Guided mode a
    # missing_dof outcome is committed right away with the FF-optimal values.
    start = _ok(client.post(f"/api/sidechains/start/{ws}", json=prepare_request), "sidechains/start")
    if start["status"] == "missing_dof":
        r["sidechains_start"] = start
        r["prepare"] = _ok(client.post(f"/api/sidechains/commit/{ws}"), "sidechains/commit")
    else:
        r["prepare"] = start
    result.files["prepared.pdb"] = _read_ws(client, ws, "structure.pdb")
    result.files["structure.forge_meta.json"] = _read_ws(client, ws, "structure.forge_meta.json")
    r["sequence_prepared"] = _ok(client.get(f"/api/analysis/sequence/{ws}?fill_gaps=true"), "sequence")

    r["topology"] = _ok(client.post(f"/api/topology/{ws}/generate", json={
        "pdb_filename": "structure.pdb",
        "ff_selections": ff_selections,
        "hmr": case.hmr,
    }), "topology")
    for name in ("structure.pdb", "structure.prmtop", "structure.crd"):
        result.files[name] = _read_ws(client, ws, name)

    # viewer.ts loads structures through the GET download endpoint (light = without bulk water).
    for light in (True, False):
        viewer = client.get(f"/api/download/{ws}", params={"light": str(light).lower(), "filename": "structure.pdb"})
        assert viewer.status_code == 200, viewer.text
        result.files[f"viewer_light_{str(light).lower()}.pdb"] = viewer.text

    mdin_request = {"hmr": case.hmr, "dt_fs": 4.0 if case.hmr else 2.0}
    mdin = client.post(f"/api/simulation/{ws}/amber-mdin", json=mdin_request)
    assert mdin.status_code == 200, mdin.text
    result.files["production.mdin"] = mdin.text
    r["mdin_headers"] = {"content-disposition": mdin.headers.get("content-disposition")}

    export = client.post(f"/api/download/{ws}/export", json={
        "wants_pdb": True, "wants_top": True, "wants_crd": True, "mdin": mdin_request, "as_zip": True,
    })
    assert export.status_code == 200, export.text
    with zipfile.ZipFile(io.BytesIO(export.content)) as zf:
        names = sorted(zf.namelist())
        r["export_zip"] = {"files": names}
        for name in names:
            result.files[f"export/{name}"] = zf.read(name).decode("utf-8")

    result.responses = r
    result.workspace_id = ws  # type: ignore[attr-defined]
    return result


def _read_ws(client: TestClient, ws: str, filename: str) -> str:
    from app.workspaces.manager import workspace_manager
    return workspace_manager.get_file_path(ws, filename).read_text(encoding="utf-8")


def _summarize_ff_response(ff_response: Dict[str, Any]) -> Dict[str, Any]:
    """The FF payloads carry whole force-field files; keep only what identifies the choice."""
    def brief(ff: Dict[str, Any]) -> Dict[str, Any]:
        return {k: ff.get(k) for k in ("id", "display_name", "molecule_type", "tier", "is_default")}

    return {
        "mode": ff_response["mode"],
        "detected_types": ff_response["detected_types"],
        "required_groups": ff_response["required_groups"],
        "forcefields_by_group": {g: [brief(f) for f in ffs] for g, ffs in ff_response["forcefields_by_group"].items()},
        "water_profiles": [
            {k: v for k, v in p.items() if k not in ("water_ff", "ion_ffs")} | {"water_ff": brief(p["water_ff"])}
            for p in ff_response["water_profiles"]
        ],
        "forcefields_count": len(ff_response["forcefields"]),
    }
