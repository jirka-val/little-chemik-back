"""
AMBER_topology.getBox (box z prmtopu po cpptraj) a savePdb (PDB z topologie
a souřadnic) - malé ručně sestavené topologie.
"""

import pytest

from app.services.structure.solvated_system import box_shape
from app.utils.adams4sims_processing_library.utils import AMBER_topology

pytestmark = pytest.mark.unit


def _pointers(ifbox: int) -> list:
    pointers = [0] * 31
    pointers[27] = ifbox
    return pointers


@pytest.mark.parametrize("beta, expected_angles, shape", [
    (90.0, (90.0, 90.0, 90.0), "cubic"),
    (109.4712190, (109.4712190,) * 3, "truncated octahedron"),
    (60.0, (60.0, 90.0, 60.0), "rhombic dodecahedron"),
])
def test_get_box_decodes_box_dimensions(beta, expected_angles, shape):
    top = {"POINTERS": _pointers(1), "BOX_DIMENSIONS": [beta, 50.0, 50.0, 50.0]}
    box = AMBER_topology.getBox(top)
    assert box[:3] == [50.0, 50.0, 50.0]
    assert tuple(box[3:]) == pytest.approx(expected_angles)
    assert box_shape(dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), box))) == shape


def test_get_box_without_periodicity():
    assert AMBER_topology.getBox({"POINTERS": _pointers(0), "BOX_DIMENSIONS": [90.0, 1, 1, 1]}) == []
    assert AMBER_topology.getBox({"POINTERS": _pointers(1)}) == []


def test_get_box_rejects_invalid_dimensions():
    with pytest.raises(ValueError):
        AMBER_topology.getBox({"POINTERS": _pointers(1), "BOX_DIMENSIONS": [90.0, 50.0, -1.0, 50.0]})
    with pytest.raises(ValueError):
        AMBER_topology.getBox({"POINTERS": _pointers(1), "BOX_DIMENSIONS": [90.0, 50.0]})


def _system():
    """Dipeptid (2 rezidua spojená vazbou), Na+ a jedna voda."""
    top = {
        "POINTERS": _pointers(1),
        "BOX_DIMENSIONS": [109.4712190, 40.0, 40.0, 40.0],
        "ATOM_NAME": ["N", "CA", "C", "N", "CA", "C", "Na+", "O", "H1", "H2"],
        "RESIDUE_LABEL": ["ALA", "GLY", "Na+", "WAT"],
        "RESIDUE_POINTER": [1, 4, 7, 8],
        "ATOMIC_NUMBER": [7, 6, 6, 7, 6, 6, 11, 8, 1, 1],
        # AMBER: indexy atomů krát 3, trojice (i, j, typ)
        "BONDS_WITHOUT_HYDROGEN": [0, 3, 1, 3, 6, 1, 6, 9, 1, 9, 12, 1, 12, 15, 1],
        "BONDS_INC_HYDROGEN": [21, 24, 1, 21, 27, 1],
    }
    coords = [(float(i), 2.0 * i, -1.5 * i) for i in range(10)]
    return top, coords


def test_save_pdb_writes_atoms_ter_and_cryst1(tmp_path):
    top, coords = _system()
    path = tmp_path / "out.pdb"
    AMBER_topology.savePdb(str(path), top, coords)
    lines = path.read_text().splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    assert len(atoms) == 10
    assert lines[0].startswith("CRYST1") and "109.47" in lines[0]
    # solut, pak ionty, pak voda - každá část končí TER
    assert [line[17:20].strip() for line in atoms] == ["ALA"] * 3 + ["GLY"] * 3 + ["Na+", "WAT", "WAT", "WAT"]
    assert sum(1 for line in lines if line.startswith("TER")) == 3
    x, y, z = float(atoms[4][30:38]), float(atoms[4][38:46]), float(atoms[4][46:54])
    assert (x, y, z) == pytest.approx(coords[4])
    # vstup se nemění
    assert top["ATOM_NAME"][0] == "N" and coords[0] == (0.0, 0.0, -0.0)


def test_save_pdb_rejects_wrong_coordinate_count(tmp_path):
    top, coords = _system()
    with pytest.raises(ValueError):
        AMBER_topology.savePdb(str(tmp_path / "out.pdb"), top, coords[:-1])
