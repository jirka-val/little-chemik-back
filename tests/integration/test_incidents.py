"""
Integrační testy hlášení chyb: deník akcí + snímky ve workspace
(app/core/incident_middleware.py, app/services/incidents/history.py) a
ukládání reportů (app/services/incidents/store.py, /api/incidents).
"""

import gzip
import json
from pathlib import Path

import pytest

from app.core.config import settings
from app.services import analysis_service
from app.services.incidents import history
from app.workspaces.manager import WORKSPACE_DIR

pytestmark = pytest.mark.integration


@pytest.fixture
def incidents_tmp(monkeypatch, tmp_path):
    monkeypatch.setattr(settings, "INCIDENTS_DIR", tmp_path / "incidents")
    monkeypatch.setattr(settings, "INCIDENT_REPORTS_ENABLED", True)
    monkeypatch.setattr(settings, "ADMIN_TOKEN", "secret")
    return tmp_path / "incidents"


def _upload(client, pdb_text: str) -> str:
    files = {"file": ("structure.pdb", pdb_text.encode(), "chemical/x-pdb")}
    return client.post("/api/molecules/upload", files=files).json()["workspace_id"]


def _history(ws_id: str) -> Path:
    return Path(WORKSPACE_DIR) / ws_id / history.HISTORY_DIRNAME


def _journal(ws_id: str):
    return [json.loads(l) for l in (_history(ws_id) / history.JOURNAL_FILE).read_text(encoding="utf-8").splitlines()]


def test_journal_records_settings_and_snapshots_structure(client, incidents_tmp, pdb_altloc_sample):
    ws_id = _upload(client, pdb_altloc_sample)
    response = client.post(f"/api/analysis/clean-altlocs/{ws_id}", json={"selection": {"A_42_SER": "B"}})
    assert response.status_code == 200

    journal = _journal(ws_id)
    assert journal[0]["event"] == "workspace_created"
    step = next(e for e in journal if e.get("path", "").endswith(f"/clean-altlocs/{ws_id}"))
    assert step["body"]["selection"] == {"A_42_SER": "B"}
    assert step["status"] == 200

    # Originál po uploadu; clean-altlocs strukturu přepsal, ale snímek se
    # dělá PŘED krokem - změna se zachytí až před dalším krokem.
    snaps = sorted(p.name for p in (_history(ws_id) / history.SNAPSHOT_DIRNAME).iterdir())
    assert snaps == ["000_upload_structure.pdb.gz"]
    client.post(f"/api/analysis/clean-altlocs/{ws_id}", json={"selection": {}})
    snaps = sorted(p.name for p in (_history(ws_id) / history.SNAPSHOT_DIRNAME).iterdir())
    assert snaps == ["000_upload_structure.pdb.gz", "001_analysis-clean-altlocs_structure.pdb.gz"]
    original = gzip.decompress((_history(ws_id) / history.SNAPSHOT_DIRNAME / snaps[0]).read_bytes()).decode()
    assert original == pdb_altloc_sample


def test_server_error_is_recorded_with_original_traceback(client, incidents_tmp, pdb_altloc_sample, monkeypatch):
    ws_id = _upload(client, pdb_altloc_sample)

    def boom(**kwargs):
        raise ValueError("Fixed system charge -56.691900 is not sufficiently close to an integer")

    monkeypatch.setattr("app.api.v1.endpoints.analysis.process_structure", boom)
    response = client.post(f"/api/analysis/clean-altlocs/{ws_id}", json={"selection": {}})
    assert response.status_code == 500

    errors = [json.loads(l) for l in (_history(ws_id) / history.ERRORS_FILE).read_text(encoding="utf-8").splitlines()]
    assert len(errors) == 1
    assert errors[0]["status"] == 500
    assert "not sufficiently close to an integer" in errors[0]["traceback"]


def test_report_requires_consent(client, incidents_tmp):
    response = client.post("/api/incidents", json={"kind": "manual", "description": "x", "consent": False})
    assert response.status_code == 400
    assert not incidents_tmp.exists() or not any(incidents_tmp.iterdir())


def test_manual_report_bundles_workspace_history(client, incidents_tmp, pdb_altloc_sample):
    ws_id = _upload(client, pdb_altloc_sample)
    client.post(f"/api/analysis/clean-altlocs/{ws_id}", json={"selection": {"A_42_SER": "B"}})

    response = client.post("/api/incidents", json={
        "workspace_id": ws_id,
        "kind": "manual",
        "description": "Waters form a weird cluster at the box edge.",
        "consent": True,
        "client_log": [{"type": "click", "label": "Apply & Clean Structure"}],
        "page": {"userAgent": "pytest"},
    })
    assert response.status_code == 200
    incident_id = response.json()["incident_id"]

    folder = incidents_tmp / incident_id
    meta = json.loads((folder / "meta.json").read_text(encoding="utf-8"))
    assert meta["kind"] == "manual"
    assert meta["workspace_found"] is True
    assert (folder / "journal.jsonl").is_file()
    assert (folder / "client_log.json").is_file()
    assert (folder / "current" / "structure.pdb.gz").is_file()
    assert (folder / "snapshots" / "000_upload_structure.pdb.gz").is_file()


def test_admin_list_and_download(client, incidents_tmp):
    incident_id = client.post("/api/incidents", json={"kind": "manual", "description": "x", "consent": True}).json()["incident_id"]

    assert client.get("/api/incidents").status_code == 403
    listing = client.get("/api/incidents", headers={"X-Admin-Token": "secret"}).json()
    assert [i["id"] for i in listing["incidents"]] == [incident_id]

    download = client.get(f"/api/incidents/{incident_id}/download", headers={"X-Admin-Token": "secret"})
    assert download.status_code == 200
    assert download.headers["content-type"] == "application/zip"


def test_admin_delete(client, incidents_tmp):
    incident_id = client.post("/api/incidents", json={"kind": "manual", "description": "x", "consent": True}).json()["incident_id"]
    admin = {"X-Admin-Token": "secret"}

    assert client.delete(f"/api/incidents/{incident_id}").status_code == 403
    assert (incidents_tmp / incident_id).is_dir()

    response = client.delete(f"/api/incidents/{incident_id}", headers=admin)
    assert response.status_code == 200
    assert response.json() == {"deleted": incident_id}
    assert not (incidents_tmp / incident_id).exists()
    assert client.get("/api/incidents", headers=admin).json()["incidents"] == []

    assert client.delete(f"/api/incidents/{incident_id}", headers=admin).status_code == 404
    # ID mimo formát (např. pokus o cestu ven ze složky) se nesmaže
    assert client.delete("/api/incidents/..", headers=admin).status_code in (404, 405)


def test_storage_limit_rejects_new_reports(client, incidents_tmp, monkeypatch):
    monkeypatch.setattr(settings, "INCIDENTS_MAX_TOTAL_MB", 0)
    response = client.post("/api/incidents", json={"kind": "manual", "description": "x", "consent": True})
    assert response.status_code == 507
    assert response.json()["code"] == "incident_storage_full"


def test_disabled_turns_everything_off(client, incidents_tmp, monkeypatch, pdb_altloc_sample):
    monkeypatch.setattr(settings, "INCIDENT_REPORTS_ENABLED", False)
    assert client.get("/api/incidents/config").json() == {"enabled": False}
    response = client.post("/api/incidents", json={"kind": "manual", "description": "x", "consent": True})
    assert response.status_code == 404

    ws_id = _upload(client, pdb_altloc_sample)
    client.post(f"/api/analysis/clean-altlocs/{ws_id}", json={"selection": {}})
    assert not _history(ws_id).exists()
