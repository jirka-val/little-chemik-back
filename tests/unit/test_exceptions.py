"""
Unit testy pro jednotnou chybovou obálku (app/core/exceptions.py) a
WorkspaceManager.require_workspace - viz frontend refaktor chyb v konverzaci
(backend mělo 3 různé tvary chybové odpovědi současně, frontend je neuměl
spolehlivě parsovat).
"""

import pytest

from app.core.exceptions import (
    AppBaseException,
    BadRequestError,
    ExternalServiceError,
    InternalError,
    NotFoundError,
    RemoteMoleculeNotFoundError,
    WorkspaceNotFoundError,
)
from app.workspaces.manager import WorkspaceManager

pytestmark = pytest.mark.unit


class TestExceptionDefaults:
    @pytest.mark.parametrize("exc_cls,expected_status,expected_code", [
        (WorkspaceNotFoundError, 404, "workspace_not_found"),
        (NotFoundError, 404, "not_found"),
        (ExternalServiceError, 502, "external_service_error"),
        (BadRequestError, 400, "bad_request"),
        (InternalError, 500, "internal_error"),
    ])
    def test_class_level_status_and_code_defaults(self, exc_cls, expected_status, expected_code):
        exc = exc_cls("test message")
        assert exc.status_code == expected_status
        assert exc.code == expected_code
        assert exc.message == "test message"

    def test_workspace_not_found_has_sensible_default_message(self):
        exc = WorkspaceNotFoundError()
        assert "workspace" in exc.message.lower()

    def test_remote_molecule_not_found_carries_pdb_code_in_payload(self):
        exc = RemoteMoleculeNotFoundError("9XYZ")
        assert exc.status_code == 404
        assert exc.code == "remote_molecule_not_found"
        assert exc.payload["pdb_code"] == "9XYZ"
        assert "9XYZ" in exc.message

    def test_instance_can_override_class_defaults(self):
        exc = ExternalServiceError("IDA API down", status_code=503)
        assert exc.status_code == 503  # instance override vyhrává nad class-level 502

    def test_payload_merges_into_json_response_alongside_error_and_code(self):
        """
        Ověřuje přesně ten kontrakt, na který se bude spoléhat frontend:
        {"error": <třída>, "code": <stabilní_identifikátor>, "message": <text>, **payload}
        """
        exc = AppBaseException("boom", status_code=409, code="custom_conflict", payload={"chain": "A", "resseq": 5})
        content = {"error": exc.__class__.__name__, "code": exc.code, "message": exc.message, **exc.payload}
        assert content == {
            "error": "AppBaseException",
            "code": "custom_conflict",
            "message": "boom",
            "chain": "A",
            "resseq": 5,
        }


class TestRequireWorkspace:
    def test_raises_workspace_not_found_for_missing_workspace(self, tmp_path, monkeypatch):
        import app.workspaces.manager as wm_module
        monkeypatch.setattr(wm_module, "WORKSPACE_DIR", str(tmp_path))
        manager = WorkspaceManager()

        with pytest.raises(WorkspaceNotFoundError):
            manager.require_workspace("does-not-exist")

    def test_does_not_raise_for_existing_workspace(self, tmp_path, monkeypatch):
        import app.workspaces.manager as wm_module
        monkeypatch.setattr(wm_module, "WORKSPACE_DIR", str(tmp_path))
        manager = WorkspaceManager()
        ws_id = manager.create_from_string("ATOM ...")

        manager.require_workspace(ws_id)  # nesmí vyhodit
