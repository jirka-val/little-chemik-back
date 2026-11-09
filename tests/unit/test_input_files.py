"""
Souřadnice nahrané v upload dialogu a souhrn vstupů (app/services/input_files.py).
Čtení prmtopu pokrývá tests/golden/test_prmtop_consistency.py.
"""

import struct
from pathlib import Path

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


def _netcdf_restart(xyz, box=None, version: int = 2) -> bytes:
    """Minimální AMBER NetCDF3 restart (jako pmemd ncrst): coordinates
    (atom, spatial), volitelně cell_lengths / cell_angles, globální atribut."""
    def name(text: str) -> bytes:
        raw = text.encode()
        return struct.pack(">i", len(raw)) + raw + b"\0" * (-len(raw) % 4)

    dims = [("spatial", 3), ("atom", len(xyz)), ("cell_spatial", 3), ("cell_angular", 3), ("label", 5)]
    variables = [("coordinates", (1, 0), [v for row in xyz for v in row])]
    if box:
        variables += [("cell_lengths", (2,), box[:3]), ("cell_angles", (3,), box[3:])]
    off = ">i" if version == 1 else ">q"

    def header(begins):
        out = b"CDF" + bytes([version]) + struct.pack(">i", 0)
        out += struct.pack(">ii", 0x0A, len(dims)) + b"".join(name(n) + struct.pack(">i", s) for n, s in dims)
        title = b"AMBERRESTART"
        out += struct.pack(">ii", 0x0C, 1) + name("Conventions") + struct.pack(">ii", 2, len(title)) + title
        out += struct.pack(">ii", 0x0B, len(variables))
        for (var, dim_ids, values), begin in zip(variables, begins):
            out += name(var) + struct.pack(">i", len(dim_ids)) + struct.pack(f">{len(dim_ids)}i", *dim_ids)
            out += struct.pack(">ii", 0, 0) + struct.pack(">ii", 6, 8 * len(values)) + struct.pack(off, begin)
        return out

    start = len(header([0] * len(variables)))
    begins, data = [], b""
    for _, _, values in variables:
        begins.append(start + len(data))
        data += struct.pack(f">{len(values)}d", *values)
    return header(begins) + data


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


@pytest.mark.parametrize("version", [1, 2])
def test_netcdf_restart_coordinates_and_box(tmp_path, version):
    xyz = [(1.0, 2.0, 3.0), (-4.5, 5.25, 6.0)]
    path = tmp_path / "c.ncrst"
    path.write_bytes(_netcdf_restart(xyz, box=[50.0, 51.0, 52.0, 70.53, 109.47, 70.53], version=version))
    coords = input_files.read_coordinates(path)
    assert coords["format"] == "NetCDF restart" and coords["atoms"] == 2
    assert coords["xyz"].tolist() == [list(row) for row in xyz]
    assert coords["box"] == {"a": 50.0, "b": 51.0, "c": 52.0, "alpha": 70.53, "beta": 109.47, "gamma": 70.53}
    assert input_files.summarize_coordinates(path)["box"]["shape"] == "truncated octahedron"


def test_netcdf_restart_without_box(tmp_path):
    path = tmp_path / "c.ncrst"
    path.write_bytes(_netcdf_restart([(0.0, 0.0, 0.0)] * 1234))
    summary = input_files.summarize_coordinates(path)
    assert summary["ok"] and summary["atoms"] == 1234 and summary["box"] is None


def test_netcdf4_is_rejected_with_a_hint(tmp_path):
    path = tmp_path / "c.ncrst"
    path.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\0" * 64)
    summary = input_files.summarize_coordinates(path)
    assert not summary["ok"] and "cpptraj" in summary["error"]


def test_ascii_coordinates_are_read_with_the_full_box(tmp_path):
    path = tmp_path / "c.rst7"
    text = _inpcrd(3, box=False).splitlines()
    text[2] = f"{1.5:12.7f}{-2.25:12.7f}{3.0:12.7f}{0.0:12.7f}{0.0:12.7f}{0.0:12.7f}"
    text.append("".join(f"{v:12.7f}" for v in (51.551, 51.551, 51.551, 70.53, 109.47, 70.53)))
    path.write_text("\n".join(text) + "\n")
    coords = input_files.read_coordinates(path)
    assert coords["xyz"][0].tolist() == [1.5, -2.25, 3.0]
    assert (coords["box"]["alpha"], coords["box"]["beta"], coords["box"]["gamma"]) == (70.53, 109.47, 70.53)


def test_ascii_coordinates_with_wrong_value_count(tmp_path):
    path = tmp_path / "c.inpcrd"
    path.write_text("title\n    5\n   1.0000000   2.0000000\n")
    with pytest.raises(input_files.CoordinateError, match="should hold 5 atoms"):
        input_files.read_coordinates(path)


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


# ---------------------------------------------------------------------------
# PDB z topologie a souřadnic (bez PDB)
# ---------------------------------------------------------------------------

def _topology(box_beta=109.4712190):
    """Dipeptid s AMBER terminály (NALA, CGLY), Na+ a voda TIP3."""
    pointers = [0] * 31
    pointers[0], pointers[27] = 10, 1
    return {
        "POINTERS": pointers,
        "BOX_DIMENSIONS": [box_beta, 40.0, 40.0, 40.0],
        "ATOM_NAME": ["N", "CA", "C", "N", "CA", "C", "Na+", "O", "H1", "H2"],
        "RESIDUE_LABEL": ["NALA", "CGLY", "Na+", "TIP3"],
        "RESIDUE_POINTER": [1, 4, 7, 8],
        "ATOMIC_NUMBER": [7, 6, 6, 7, 6, 6, 11, 8, 1, 1],
        "BONDS_WITHOUT_HYDROGEN": [0, 3, 1, 3, 6, 1, 6, 9, 1, 9, 12, 1, 12, 15, 1],
        "BONDS_INC_HYDROGEN": [21, 24, 1, 21, 27, 1],
    }


def _rst7(natom: int, box=None) -> str:
    values = [float(i) for i in range(natom * 3)]
    lines = ["sim", f"{natom:6d}"]
    lines += ["".join(f"{v:12.7f}" for v in values[i:i + 6]) for i in range(0, len(values), 6)]
    if box:
        lines.append("".join(f"{v:12.7f}" for v in box))
    return "\n".join(lines) + "\n"


@pytest.fixture
def topology(monkeypatch):
    top = _topology()
    monkeypatch.setattr(input_files.AMBER_topology, "read_AMBER_topology", lambda _path: dict(top))
    return top


def _write(tmp_path, coordinates: str):
    (tmp_path / "c.rst7").write_text(coordinates)
    pdb = tmp_path / "structure.pdb"
    info = input_files.write_pdb_from_inputs(tmp_path / "t.prmtop", tmp_path / "c.rst7", pdb)
    return info, pdb.read_text().splitlines()


def test_pdb_from_inputs_takes_the_box_from_the_coordinates(tmp_path, topology):
    # Po NPT je box v rst7 jiný než v prmtopu - a umí všechny tři úhly.
    info, lines = _write(tmp_path, _rst7(10, box=(41.5, 41.5, 41.5, 70.53, 109.47, 70.53)))
    assert info["box_from"] == "coordinates"
    assert lines[0].startswith("CRYST1   41.500   41.500   41.500  70.53 109.47  70.53 P 1")
    atoms = [line for line in lines if line.startswith("ATOM")]
    assert len(atoms) == 10
    # terminály a voda pod 3znakovými jmény, souřadnice beze změny
    assert [atoms[i][17:20] for i in (0, 3, 6, 7)] == ["ALA", "GLY", "Na+", "WAT"]
    assert float(atoms[1][30:38]) == 3.0
    assert info["renamed"] == {"CGLY -> GLY": 1, "NALA -> ALA": 1, "TIP3 -> WAT": 1}
    assert sum(line.startswith("CRYST1") for line in lines) == 1


def test_pdb_from_inputs_falls_back_to_the_topology_box(tmp_path, topology):
    info, lines = _write(tmp_path, _rst7(10))
    assert info["box_from"] == "topology"
    assert lines[0].startswith("CRYST1   40.000   40.000   40.000 109.47 109.47 109.47")


def test_pdb_from_inputs_without_any_box(tmp_path, topology):
    topology["POINTERS"][27] = 0
    info, lines = _write(tmp_path, _rst7(10))
    assert info["box_from"] is None and info["box"] is None
    assert not any(line.startswith("CRYST1") for line in lines)


def test_pdb_from_inputs_rejects_mismatched_atom_count(tmp_path, topology):
    with pytest.raises(input_files.InputsError, match="10 atoms but the coordinates have 9"):
        _write(tmp_path, _rst7(9))


def test_pdb_from_inputs_rejects_long_residue_names(tmp_path, topology):
    topology["RESIDUE_LABEL"][0] = "LIG1"
    with pytest.raises(input_files.InputsError, match="LIG1"):
        _write(tmp_path, _rst7(10))


def test_upload_inputs_endpoint_creates_a_workspace(client, topology):
    response = client.post("/api/molecules/upload-inputs", files={
        "topology": ("sim.prmtop", b"ignored - read_AMBER_topology is patched"),
        "coordinates": ("sim.rst7", _rst7(10, box=(41.5, 41.5, 41.5, 90.0, 90.0, 90.0))),
    })
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["built"]["box_from"] == "coordinates" and data["built"]["atoms"] == 10
    assert data["summary"]["pdb"]["atoms"] == 10
    assert data["summary"]["coordinates"]["filename"] == "sim.rst7"
    from app.workspaces.manager import workspace_manager
    pdb = workspace_manager.get_file_path(data["workspace_id"]).read_text()
    assert pdb.startswith("CRYST1   41.500   41.500   41.500  90.00  90.00  90.00")


def test_upload_inputs_endpoint_reports_errors_and_cleans_up(client, topology, monkeypatch):
    from app.workspaces import manager
    before = set(Path(manager.WORKSPACE_DIR).iterdir())
    response = client.post("/api/molecules/upload-inputs", files={
        "topology": ("sim.prmtop", b"x"), "coordinates": ("sim.rst7", _rst7(9)),
    })
    assert response.status_code == 400 and response.json()["code"] == "inputs_not_convertible"
    assert set(Path(manager.WORKSPACE_DIR).iterdir()) == before

    missing = client.post("/api/molecules/upload-inputs", files={"topology": ("sim.prmtop", b"x")})
    assert missing.status_code in (400, 422)
