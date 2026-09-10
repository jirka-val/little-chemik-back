"""
Integrační testy pro /api/validation/check, /preview-selection, /apply-selections.
"""

import pytest

pytestmark = pytest.mark.integration


def _upload(client, pdb_text: str, filename: str = "structure.pdb") -> str:
    files = {"file": (filename, pdb_text.encode(), "chemical/x-pdb")}
    return client.post("/api/molecules/upload", files=files).json()["workspace_id"]


class TestCheck:
    def test_check_reports_altlocs(self, client, pdb_altloc_sample):
        ws_id = _upload(client, pdb_altloc_sample)
        response = client.post("/api/validation/check", json={"workspace_id": ws_id})
        assert response.status_code == 200
        data = response.json()
        assert data["analysis"]["metadata"]["has_alt_locs"] is True
        assert any(e["issue"] == "alt_locs_detected" for e in data["analysis"]["errors"])

    def test_check_nonexistent_workspace_returns_404(self, client):
        response = client.post("/api/validation/check", json={"workspace_id": "nope"})
        assert response.status_code == 404


class TestApplySelections:
    def test_apply_selections_persists_choice_to_disk(self, client, pdb_altloc_sample):
        ws_id = _upload(client, pdb_altloc_sample)

        response = client.post("/api/validation/apply-selections", json={
            "workspace_id": ws_id,
            "selections": {"A_42_SER": "A"},
        })
        assert response.status_code == 200

        download = client.get(f"/api/download/{ws_id}")
        og_lines = [l for l in download.text.splitlines() if l[12:16].strip() == "OG"]
        assert len(og_lines) == 1
        assert og_lines[0][16] == " "  # altloc indikátor odstraněn po zápisu

    def test_apply_selections_nonexistent_workspace_returns_404(self, client):
        response = client.post("/api/validation/apply-selections", json={
            "workspace_id": "nope",
            "selections": {},
        })
        assert response.status_code == 404
