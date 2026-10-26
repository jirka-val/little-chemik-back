"""
Náhled celého MD běhu po krocích (/api/simulation/{ws}/amber-protocol).
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

_RNA = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "rna_gap.pdb"


def _upload(client) -> str:
    files = {"file": ("structure.pdb", _RNA.read_bytes(), "chemical/x-pdb")}
    return client.post("/api/molecules/upload", files=files).json()["workspace_id"]


def test_preview_lists_every_file_of_the_run(client):
    ws = _upload(client)
    response = client.post(f"/api/simulation/{ws}/amber-protocol", json={"protocol": False, "steps": [
        {"name": "mini", "kind": "minimization", "restraint": "solute"},
        {"name": "heat1", "kind": "heating", "temp_target": 100},
    ]})
    assert response.status_code == 200
    files = {f["name"]: f["content"] for f in response.json()}
    assert list(files) == ["mini.in", "heat1.in", "prod.in", "run_md.sh"]
    # maska solutu ze struktury workspace, ne náhradní "vše kromě vody"
    assert "restraintmask=':1-" in files["mini.in"]
    assert "irest=1, ntx=5," in files["prod.in"]  # protocol se zapne vždy


def test_preview_reports_step_errors(client):
    ws = _upload(client)
    response = client.post(f"/api/simulation/{ws}/amber-protocol", json={"steps": [
        {"name": "h", "kind": "heating", "extra": "nstlim=5"},
    ]})
    assert response.status_code == 400
    assert response.json()["code"] == "protocol_step_extra_clash"


def test_preview_needs_a_workspace(client):
    assert client.post("/api/simulation/nope/amber-protocol", json={}).status_code == 404
