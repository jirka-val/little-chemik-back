"""
Workspace IDs and file names coming from requests must never turn into paths
outside the workspace directory, and admin endpoints must reject wrong tokens.
"""

import os
from pathlib import Path

import pytest

from app.core.config import settings
from app.core.exceptions import BadRequestError, WorkspaceNotFoundError
from app.workspaces.manager import WORKSPACE_DIR, is_valid_workspace_id, workspace_manager

pytestmark = pytest.mark.integration

CRAMBIN = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "golden" / "1CRN.pdb"
MISSING_WS = "00000000-0000-4000-8000-000000000000"


def _workspace_entries() -> set[str]:
    return set(os.listdir(WORKSPACE_DIR))


class TestWorkspaceIds:
    @pytest.mark.parametrize("bad", ["", "..", "../x", "..\\x", "does-not-exist", "A" * 36, "ABCDEF00-0000-4000-8000-000000000000"])
    def test_invalid_ids_are_rejected(self, bad):
        assert not is_valid_workspace_id(bad)
        assert not workspace_manager.workspace_exists(bad)
        with pytest.raises(WorkspaceNotFoundError):
            workspace_manager.get_workspace_dir(bad)

    def test_created_ids_are_valid(self, make_workspace):
        assert is_valid_workspace_id(make_workspace("ATOM"))

    @pytest.mark.parametrize("path", [
        "/api/analysis/sequence/..%2F..%2Fescape",
        "/api/analysis/altlocs/not-a-workspace",
        f"/api/forcefields/{MISSING_WS}",
    ])
    def test_get_with_bad_workspace_is_404_and_creates_nothing(self, client, path):
        before = _workspace_entries()
        response = client.get(path)
        assert response.status_code == 404
        assert _workspace_entries() == before

    def test_remove_residue_on_unknown_workspace_is_404_and_creates_nothing(self, client):
        before = _workspace_entries()
        response = client.post(f"/api/molecules/remove-residue/{MISSING_WS}", json={"chain": "A", "resseq": 1})
        assert response.status_code == 404
        assert _workspace_entries() == before


class TestFileNames:
    @pytest.mark.parametrize("bad", ["", ".", "..", "../structure.pdb", "sub/structure.pdb", "..\\structure.pdb"])
    def test_path_like_names_are_rejected(self, make_workspace, bad):
        ws = make_workspace("ATOM")
        with pytest.raises(BadRequestError):
            workspace_manager.get_file_path(ws, bad)

    def test_sequence_rejects_path_in_filename(self, client, make_workspace):
        ws = make_workspace(CRAMBIN.read_text(encoding="utf-8"))
        response = client.get(f"/api/analysis/sequence/{ws}", params={"filename": "../../escape.pdb"})
        assert response.status_code == 400

    def test_topology_rejects_path_in_pdb_filename(self, client, make_workspace, tmp_path):
        ws = make_workspace(CRAMBIN.read_text(encoding="utf-8"))
        response = client.post(f"/api/topology/{ws}/generate", json={
            "pdb_filename": "../escape.pdb", "ff_selections": {}, "hmr": False,
        })
        assert response.status_code == 400
        assert not (Path(WORKSPACE_DIR) / "escape.prmtop").exists()


class TestAdminToken:
    @pytest.fixture
    def admin_token(self, monkeypatch):
        monkeypatch.setattr(settings, "ADMIN_TOKEN", "s3cret-token")
        return "s3cret-token"

    def test_wrong_or_missing_token_is_403(self, client, admin_token):
        assert client.get("/api/incidents").status_code == 403
        assert client.get("/api/incidents", headers={"X-Admin-Token": "s3cret-tokeN"}).status_code == 403
        assert client.patch("/api/forcefields/classification", json={}, headers={"X-Admin-Token": "x"}).status_code == 403

    def test_right_token_is_accepted(self, client, admin_token):
        assert client.get("/api/incidents", headers={"X-Admin-Token": admin_token}).status_code == 200

    def test_empty_admin_token_locks_endpoints(self, client, monkeypatch):
        monkeypatch.setattr(settings, "ADMIN_TOKEN", "")
        assert client.get("/api/incidents", headers={"X-Admin-Token": ""}).status_code == 403
