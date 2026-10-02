"""
Integration tests for endpoints that had no coverage before the 2026-10
refactor: the sequence editor (rename/mutate/remove), residue removal, ion
options and the console log feed. They pin today's behaviour so the
refactor cannot drop or change these endpoints unnoticed.
"""

from pathlib import Path

import pytest

from app.core.logging import console_logger
from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.integration

CRAMBIN = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "golden" / "1CRN.pdb"


@pytest.fixture
def crambin_ws(make_workspace) -> str:
    return make_workspace(CRAMBIN.read_text(encoding="utf-8"))


def _residue_atoms(ws: str, chain: str, resseq: int) -> list[tuple[str, str]]:
    text = workspace_manager.get_file_path(ws).read_text(encoding="utf-8")
    return [
        (line[17:20].strip(), line[12:16].strip())
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain and int(line[22:26]) == resseq
    ]


class TestEditor:
    def test_rename_residue(self, client, crambin_ws):
        response = client.post("/api/editor/rename-residue", json={
            "workspace_id": crambin_ws, "chain_id": "A", "residue_number": 1, "new_res_name": "SER",
        })
        assert response.status_code == 200, response.text
        assert response.json() == {"message": "Success"}
        assert {resname for resname, _ in _residue_atoms(crambin_ws, "A", 1)} == {"SER"}

    def test_rename_atom(self, client, crambin_ws):
        response = client.post("/api/editor/rename-atom", json={
            "workspace_id": crambin_ws, "chain_id": "A", "residue_number": 1,
            "old_atom_name": "OG1", "new_atom_name": "OG",
        })
        assert response.status_code == 200, response.text
        names = [name for _, name in _residue_atoms(crambin_ws, "A", 1)]
        assert "OG" in names and "OG1" not in names

    def test_remove_atom(self, client, crambin_ws):
        response = client.post("/api/editor/remove-atom", json={
            "workspace_id": crambin_ws, "chain_id": "A", "residue_number": 1, "atom_name": "CG2",
        })
        assert response.status_code == 200, response.text
        names = [name for _, name in _residue_atoms(crambin_ws, "A", 1)]
        assert names == ["N", "CA", "C", "O", "CB", "OG1"]

    def test_mutate(self, client, crambin_ws):
        response = client.post("/api/editor/mutate", json={
            "workspace_id": crambin_ws, "chain_id": "A", "residue_number": 2, "mutate_to": "ALA",
        })
        assert response.status_code == 200, response.text
        atoms = _residue_atoms(crambin_ws, "A", 2)
        assert {resname for resname, _ in atoms} == {"ALA"}
        assert {"N", "CA", "C", "O", "CB"} <= {name for _, name in atoms}

    def test_rename_residue_rejects_long_name(self, client, crambin_ws):
        response = client.post("/api/editor/rename-residue", json={
            "workspace_id": crambin_ws, "chain_id": "A", "residue_number": 1, "new_res_name": "LONG",
        })
        assert response.status_code == 422

    def test_unknown_workspace_returns_404(self, client):
        response = client.post("/api/editor/rename-residue", json={
            "workspace_id": "00000000-0000-0000-0000-000000000000",
            "chain_id": "A", "residue_number": 1, "new_res_name": "SER",
        })
        assert response.status_code == 404


class TestRemoveResidue:
    def test_removes_residue(self, client, crambin_ws):
        response = client.post(f"/api/molecules/remove-residue/{crambin_ws}", json={"chain": "A", "resseq": 3})
        assert response.status_code == 200, response.text
        assert _residue_atoms(crambin_ws, "A", 3) == []
        assert _residue_atoms(crambin_ws, "A", 4) != []

    def test_missing_residue_returns_404(self, client, crambin_ws):
        response = client.post(f"/api/molecules/remove-residue/{crambin_ws}", json={"chain": "A", "resseq": 999})
        assert response.status_code == 404


class TestIonOptions:
    def test_returns_buildable_ions_per_group(self, client):
        response = client.get("/api/validation/ions")
        assert response.status_code == 200, response.text
        body = response.json()
        assert set(body) >= {"I1"}
        assert "Na+" in body["I1"] and "Cl-" in body["I1"]


class TestConsoleLogs:
    def test_returns_new_entries_since_id(self, client):
        first = client.get("/api/system/logs").json()
        console_logger.info("refactor-check message")
        second = client.get("/api/system/logs", params={"since_id": first["last_id"]}).json()
        assert [e["message"] for e in second["entries"]][-1] == "refactor-check message"
        assert second["last_id"] > first["last_id"]
