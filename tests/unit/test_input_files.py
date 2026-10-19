"""
Souřadnice nahrané v upload dialogu a souhrn vstupů (app/services/input_files.py).
Čtení prmtopu pokrývá tests/golden/test_prmtop_consistency.py.
"""

import struct

import pytest

from app.services import input_files

pytestmark = pytest.mark.unit


def _inpcrd(natom: int, box: bool, velocities: bool = False) -> str:
    coords = [0.0] * (natom * 3)
    lines = ["test system", f"{natom:6d}"]

    def block(values):
        return ["".join(f"{v:12.7f}" for v in values[i:i + 6]) for i in range(0, len(values), 6)]

    lines += block(coords)
    if velocities:
        lines += block(coords)
    if box:
        lines.append("".join(f"{v:12.7f}" for v in (40.0, 40.0, 40.0, 109.4712190, 109.4712190, 109.4712190)))
    return "\n".join(lines) + "\n"


def _netcdf_header(natom: int) -> bytes:
    def name(text: str) -> bytes:
        raw = text.encode()
        return struct.pack(">i", len(raw)) + raw + b"\0" * (-len(raw) % 4)

    dims = [("spatial", 3), ("atom", natom)]
    body = struct.pack(">ii", 0x0A, len(dims)) + b"".join(name(n) + struct.pack(">i", s) for n, s in dims)
    return b"CDF\x02" + struct.pack(">i", 0) + body


def test_ascii_coordinates_with_box(tmp_path):
    path = tmp_path / "c.inpcrd"
    path.write_text(_inpcrd(5, box=True))
    summary = input_files.summarize_coordinates(path)
    assert summary["atoms"] == 5 and summary["format"] == "ASCII coordinates"
    assert summary["box"]["beta"] == pytest.approx(109.471219)


def test_ascii_restart_with_velocities(tmp_path):
    path = tmp_path / "c.rst7"
    path.write_text(_inpcrd(4, box=True, velocities=True))
    summary = input_files.summarize_coordinates(path)
    assert summary["velocities"] is True and summary["format"] == "ASCII restart"
    assert summary["box"]["a"] == 40.0


def test_netcdf_restart_atom_count(tmp_path):
    path = tmp_path / "c.ncrst"
    path.write_bytes(_netcdf_header(1234) + b"\0" * 64)
    summary = input_files.summarize_coordinates(path)
    assert summary == {"ok": True, "format": "NetCDF restart", "atoms": 1234, "box": None}


def test_not_a_coordinate_file(tmp_path):
    path = tmp_path / "c.crd"
    path.write_text("hello\nworld\n")
    assert input_files.summarize_coordinates(path)["ok"] is False


def test_summary_warns_on_atom_count_mismatch(tmp_path):
    (tmp_path / "structure.pdb").write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\n" * 3)
    (tmp_path / "uploaded_coordinates.inpcrd").write_text(_inpcrd(5, box=False))
    input_files.save_manifest(tmp_path, {"coordinates": "sys.inpcrd"})
    summary = input_files.summarize_inputs(tmp_path)
    assert summary["pdb"]["atoms"] == 3
    assert summary["coordinates"]["filename"] == "sys.inpcrd"
    assert summary["topology"] is None


def test_upload_endpoint_rejects_unknown_extension(client):
    pdb = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\nEND\n"
    ws = client.post("/api/molecules/upload", files={"file": ("a.pdb", pdb)}).json()["workspace_id"]
    bad = client.post(f"/api/molecules/{ws}/inputs", files={"topology": ("a.txt", b"x")})
    assert bad.status_code == 400 and bad.json()["code"] == "unsupported_topology_format"

    ok = client.post(f"/api/molecules/{ws}/inputs", files={"coordinates": ("a.inpcrd", _inpcrd(1, box=False))})
    assert ok.status_code == 200
    assert ok.json()["coordinates"]["atoms"] == 1
    assert client.get(f"/api/molecules/{ws}/inputs").json()["coordinates"]["filename"] == "a.inpcrd"
