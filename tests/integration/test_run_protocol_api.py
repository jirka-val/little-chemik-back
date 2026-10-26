"""
Protokol přípravy přes API: kontrola souboru a přibalení k exportu.
"""

import io
import json
import zipfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

_RNA = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "rna_gap.pdb"
_PROTOCOL = {"format": "forge-protocol", "format_version": 1, "name": "RNA test", "extra_key": 1}


def test_check_returns_full_protocol_and_warnings(client):
    response = client.post("/api/protocols/check", json=_PROTOCOL)
    assert response.status_code == 200
    body = response.json()
    assert body["protocol"]["simulation"]["temp0"] == 298.16
    assert body["warnings"] == ["Unknown setting 'extra_key' was ignored."]
    assert body["filename"] == "RNA_test.forge-protocol.json"
    assert len(body["fingerprint"]) == 12


def test_check_rejects_invalid_protocol_with_code(client):
    response = client.post("/api/protocols/check", json={**_PROTOCOL, "preparation": {"ph": 99}})
    assert response.status_code == 400
    assert response.json()["code"] == "protocol_invalid"
    assert "preparation.ph" in response.json()["message"]


def test_export_bundles_protocol(client):
    files = {"file": ("structure.pdb", _RNA.read_bytes(), "chemical/x-pdb")}
    ws = client.post("/api/molecules/upload", files=files).json()["workspace_id"]
    protocol = client.post("/api/protocols/check", json=_PROTOCOL).json()["protocol"]

    response = client.post(f"/api/download/{ws}/export", json={
        "wants_pdb": True, "procedure": True, "protocol": protocol, "as_zip": True,
    })
    assert response.status_code == 200
    archive = zipfile.ZipFile(io.BytesIO(response.content))
    assert "RNA_test.forge-protocol.json" in archive.namelist()
    assert json.loads(archive.read("RNA_test.forge-protocol.json"))["name"] == "RNA test"
    assert "forge_procedure.pdf" in archive.namelist()
