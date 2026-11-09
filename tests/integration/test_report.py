"""
Souhrnný report přípravy (PDF v Exportu): záznamy kroků v `_report/` a jejich
vykreslení (app/services/report/).
"""

import io
import zipfile
from pathlib import Path

import pytest

from app.services.report.records import PREPARATION, SOURCE, load_records

pytestmark = pytest.mark.integration


def _prepare(client, ws_id):
    response = client.post(f"/api/sidechains/start/{ws_id}", json={
        "workspace_id": ws_id,
        "ff_selections": {"R": {"display_name": "OL3CP", "has_ghbfix": False, "reference_article_doi": "10.1/x"}},
        "ph": 6.5,
        "add_solvent": False,
    })
    assert response.status_code == 200, response.text


def test_preparation_is_recorded(client, offline_forge_ff, make_workspace, pdb_rna_gap):
    ws_id = make_workspace(pdb_rna_gap)
    _prepare(client, ws_id)

    prep = load_records(ws_id)[PREPARATION]
    assert prep["settings"]["ph"] == 6.5
    assert prep["force_fields"]["R"]["display_name"] == "OL3CP"
    assert "summary" in prep
    assert "ff_selections" not in prep["settings"]


def test_report_only_export_is_pdf(client, offline_forge_ff, make_workspace, pdb_rna_gap):
    ws_id = make_workspace(pdb_rna_gap)
    _prepare(client, ws_id)

    response = client.post(f"/api/download/{ws_id}/export", json={"report": True, "as_zip": False})
    assert response.status_code == 200, response.text
    assert response.headers["content-type"] == "application/pdf"
    assert response.content.startswith(b"%PDF")


def test_report_goes_into_zip_with_other_files(client, make_workspace, pdb_alanine_single):
    ws_id = make_workspace(pdb_alanine_single)
    response = client.post(f"/api/download/{ws_id}/export", json={
        "wants_pdb": True, "report": True, "mdin": {"duration_ns": 1.0}, "detail_level": "expert",
    })
    assert response.status_code == 200, response.text
    names = zipfile.ZipFile(io.BytesIO(response.content)).namelist()
    assert "structure.pdb" in names and "forge_report.pdf" in names
    assert any(n.endswith(".mdin") for n in names)


def test_procedure_description_export(client, offline_forge_ff, make_workspace, pdb_rna_gap):
    ws_id = make_workspace(pdb_rna_gap)
    _prepare(client, ws_id)

    single = client.post(f"/api/download/{ws_id}/export", json={"procedure": True, "as_zip": False})
    assert single.status_code == 200, single.text
    assert single.headers["content-type"] == "application/pdf"
    assert "forge_procedure.pdf" in single.headers["content-disposition"]

    both = client.post(f"/api/download/{ws_id}/export", json={"procedure": True, "report": True})
    names = zipfile.ZipFile(io.BytesIO(both.content)).namelist()
    assert {"forge_procedure.pdf", "forge_report.pdf"} <= set(names)


def test_export_errors_name_what_is_missing(client, make_workspace, pdb_alanine_single):
    ws_id = make_workspace(pdb_alanine_single)

    nothing = client.post(f"/api/download/{ws_id}/export", json={})
    assert nothing.status_code == 400
    assert nothing.json()["code"] == "nothing_selected"

    no_topology = client.post(f"/api/download/{ws_id}/export", json={"wants_top": True, "wants_crd": True})
    assert no_topology.status_code == 404
    message = no_topology.json()["message"]
    assert "structure.prmtop, structure.crd" in message
    assert "Export step" in message


def test_upload_records_structure_citation(client):
    pdb = (Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "golden" / "1RNA.pdb").read_text()
    response = client.post("/api/molecules/upload", files={"file": ("1RNA.pdb", pdb.encode(), "chemical/x-pdb")})
    assert response.status_code == 200, response.text

    citation = load_records(response.json()["workspace_id"])[SOURCE]["citation"]
    assert citation["pdb_id"] == "1RNA"
    assert citation["doi"] == "10.1016/0022-2836(89)90010-7"


def test_amide_orientations_are_described():
    from app.services.report.procedure import _amide_orientation_sentence

    text = _amide_orientation_sentence({
        "evaluated": 8,
        "flipped": [{"residue": "A49", "resname": "GLN", "forced": False},
                    {"residue": "A60", "resname": "ASN", "forced": True}],
        "kept_by_user": ["A62"],
    })
    assert "8 ASN/GLN amide groups were optimized together with the protonation states" in text
    assert "2 amides were flipped (GLN A49 and ASN A60), 1 of them as chosen in the structure review" in text
    assert "1 amide was kept as deposited by choice" in text
    assert _amide_orientation_sentence({"evaluated": 0}) == ""
