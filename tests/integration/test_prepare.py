"""
Integration tests for structure preparation through POST /api/sidechains/start,
the endpoint the frontend uses (hydrogen-panel.ts). It runs the FORGE builder
(app/builder) through ForgeStructureService; structures that need no side-chain
decision finish in one call with status "complete".

Uses the `offline_forge_ff` fixture (see conftest.py), so it touches neither
the network nor the real data/ff_cache/; the tests skip themselves when the
OL3 (RNA) force field is not cached on this machine.
"""

import json

import pytest

from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.integration


def _start(client, ws_id, **overrides):
    return client.post(f"/api/sidechains/start/{ws_id}", json={
        "workspace_id": ws_id,
        "ff_selections": {"R": {"display_name": "OL3"}},
        "ph": 7.0,
        "add_solvent": False,
        **overrides,
    })


class TestPrepareHappyPath:
    def test_prepare_writes_forge_meta_sidecar(self, client, offline_forge_ff, make_workspace, pdb_rna_gap):
        ws_id = make_workspace(pdb_rna_gap)

        response = _start(client, ws_id)

        assert response.status_code == 200, response.text
        assert response.json()["status"] == "complete"
        meta_path = workspace_manager.get_file_path(ws_id, "structure.forge_meta.json")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        # Terminals at the gap boundary come from the builder, not from the residue name.
        assert meta["A:3:"]["ff_resname"] == "RA3"
        assert meta["A:6:"]["ff_resname"] == "RU5"

    def test_prepare_updates_structure_pdb_with_built_atoms(
        self, client, offline_forge_ff, make_workspace, pdb_rna_gap
    ):
        ws_id = make_workspace(pdb_rna_gap)
        original = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text(encoding="utf-8")

        response = _start(client, ws_id)
        assert response.status_code == 200, response.text

        updated = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text(encoding="utf-8")
        original_atom_count = sum(1 for l in original.splitlines() if l.startswith("ATOM"))
        updated_atom_count = sum(1 for l in updated.splitlines() if l.startswith("ATOM"))
        assert updated_atom_count > original_atom_count  # the builder added the missing (hydrogen) atoms

    def test_response_carries_warnings_validation_and_summary(
        self, client, offline_forge_ff, make_workspace, pdb_rna_gap
    ):
        ws_id = make_workspace(pdb_rna_gap)
        body = _start(client, ws_id).json()
        assert body["message"] == "Structure successfully prepared."
        assert {"warnings", "validation", "preparation_summary"} <= set(body)


class TestPrepareValidation:
    def test_missing_ff_selections_returns_422(self, client, make_workspace, pdb_alanine_single):
        ws_id = make_workspace(pdb_alanine_single)
        response = client.post(f"/api/sidechains/start/{ws_id}", json={"workspace_id": ws_id, "ph": 7.0})
        assert response.status_code == 422

    def test_nonexistent_workspace_returns_404(self, client):
        ws_id = "00000000-0000-4000-8000-000000000000"
        response = _start(client, ws_id)
        assert response.status_code == 404
