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


def _link_pdb() -> str:
    # U1 O3' -> G2 P: konformace A navazuje (1.6 A), B ne (3.4 A).
    rows = [
        ("O3'", " ", "U", 1, 0.0, 1.00),
        ("P", "A", "G", 2, 1.6, 0.40),
        ("P", "B", "G", 2, 3.4, 0.60),
    ]
    return "".join(
        f"ATOM  {i + 1:5d} {name:<4s}{alt:1s}{res:>3s} A{seq:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f} 20.00           {name[0]}\n"
        for i, (name, alt, res, seq, x, occ) in enumerate(rows)
    )


class TestAltlocBreaks:
    def test_reports_break_created_by_selection(self, client):
        ws_id = _upload(client, _link_pdb())
        response = client.post(f"/api/analysis/altloc-breaks/{ws_id}", json={"selection": {"A_2_G": "B"}})
        assert response.status_code == 200
        breaks = response.json()["breaks"]
        assert len(breaks) == 1
        assert breaks[0]["prev"]["resseq"] == 1
        assert breaks[0]["next"] == {"resseq": 2, "icode": "", "resname": "G", "altloc": "B"}

    def test_connected_selection_has_no_breaks_and_is_recommended(self, client):
        ws_id = _upload(client, _link_pdb())
        response = client.post(f"/api/analysis/altloc-breaks/{ws_id}", json={"selection": {"A_2_G": "A"}})
        assert response.json()["breaks"] == []
        analysis = client.get(f"/api/analysis/altlocs/{ws_id}").json()
        assert analysis["residues"][0]["recommended_alt"] == "A"
