"""
Integrační testy pro základní životní cyklus workspace: upload, stažení,
sekvenční analýza. Žádná síť - vše z lokálních fixture PDB souborů.

Síťové varianty (fetch-pdb z RCSB) jsou v tests/network/test_rcsb_smoke.py.
"""

import pytest

pytestmark = pytest.mark.integration


class TestUpload:
    def test_upload_valid_pdb_returns_workspace_id(self, client, pdb_alanine_single):
        files = {"file": ("test_molecule.pdb", pdb_alanine_single.encode(), "chemical/x-pdb")}
        response = client.post("/api/molecules/upload", files=files)

        assert response.status_code == 200
        data = response.json()
        assert "workspace_id" in data
        assert data["filename"] == "test_molecule.pdb"

    def test_upload_rejects_non_pdb_extension(self, client):
        files = {"file": ("invalid.txt", b"not a pdb", "text/plain")}
        response = client.post("/api/molecules/upload", files=files)

        assert response.status_code == 400
        body = response.json()
        assert body["code"] == "bad_request"
        assert "pdb" in body["message"].lower()


class TestDownload:
    def test_download_returns_uploaded_content(self, client, pdb_alanine_single):
        files = {"file": ("m.pdb", pdb_alanine_single.encode(), "chemical/x-pdb")}
        upload = client.post("/api/molecules/upload", files=files)
        ws_id = upload.json()["workspace_id"]

        response = client.get(f"/api/download/{ws_id}")
        assert response.status_code == 200
        assert response.headers["content-type"] == "chemical/x-pdb"
        assert "ALA" in response.text

    def test_download_nonexistent_workspace_returns_404(self, client):
        response = client.get("/api/download/non-existent-uuid-12345")
        assert response.status_code == 404


class TestSequenceAnalysis:
    def test_sequence_endpoint_returns_tokens(self, client, pdb_alanine_single):
        files = {"file": ("alanine.pdb", pdb_alanine_single.encode(), "chemical/x-pdb")}
        upload = client.post("/api/molecules/upload", files=files)
        ws_id = upload.json()["workspace_id"]

        response = client.get(f"/api/analysis/sequence/{ws_id}")
        assert response.status_code == 200
        data = response.json()
        tokens = data["sequence"]["chains"]["A"]["tokens"]
        assert tokens[0]["pdb_resname"] == "ALA"

    def test_sequence_analysis_detects_gap_terminus(self, client, pdb_protein_gap):
        ws_id = client.post(
            "/api/molecules/upload",
            files={"file": ("gap.pdb", pdb_protein_gap.encode(), "chemical/x-pdb")},
        ).json()["workspace_id"]

        response = client.get(f"/api/analysis/sequence/{ws_id}")
        assert response.status_code == 200
        tokens = response.json()["sequence"]["chains"]["K"]["tokens"]

        glu83 = next(t for t in tokens if t.get("resseq") == 83)
        phe89 = next(t for t in tokens if t.get("resseq") == 89)
        assert glu83["ff_resname"] == "CGLU"
        assert phe89["ff_resname"] == "NPHE"

    def test_sequence_analysis_nonexistent_workspace_returns_404(self, client):
        response = client.get("/api/analysis/sequence/does-not-exist")
        assert response.status_code == 404


class TestAltlocsEndpoint:
    def test_altlocs_endpoint_reports_variants(self, client, pdb_altloc_sample):
        ws_id = client.post(
            "/api/molecules/upload",
            files={"file": ("altloc.pdb", pdb_altloc_sample.encode(), "chemical/x-pdb")},
        ).json()["workspace_id"]

        response = client.get(f"/api/analysis/altlocs/{ws_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["hasAltLocs"] is True
        assert data["residues"][0]["resseq"] == 42
