"""
Regresní testy pro systémy s více než 99 999 atomy (incident 20261001-061334-ce249e).

2TRA v kubickém boxu měla po solvataci 120 096 atomů. PDBService.reorder_and_inject_eps
psal serial jako 5místné číslo bez přetočení, od atomu 100 000 se tak všechny
další sloupce posunuly o znak. parse_pdb_to_topology_dict i crd writer pak
takové řádky tiše zahodily - prmtop i crd měly jen 99 999 atomů a u jedné
stěny boxu zůstal ~18 Å pás bez vody.
"""

import pytest

from app.services.analysis import process_structure
from app.services.pdb_service import PDBService, parse_pdb_to_topology_dict
from app.services.topology_service import TopologyService

pytestmark = pytest.mark.unit

N_WATERS = 34_000  # 102 000 atomů - přes hranici 99 999
WATERS_PER_CHAIN = 9_000  # resseq má v PDB jen 4 sloupce


def _water_line(serial: int, name: str, chain: str, resseq: int, x: float, element: str) -> str:
    return (
        f"HETATM{serial % 100000:5d} {name:<4} WAT {chain}{resseq:4d}    "
        f"{x:8.3f}{1.0:8.3f}{2.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2}"
    )


def _water_pdb(n_waters: int, header: tuple[str, ...] = ()) -> str:
    lines = list(header)
    lines.append("CRYST1  108.239  108.239  108.239  90.00  90.00  90.00 P 1           1")
    serial = 1
    for index in range(n_waters):
        chain = "WXYZ"[index // WATERS_PER_CHAIN]
        resseq = index % WATERS_PER_CHAIN + 1
        x = float(index % 1000) / 10.0
        for name, element in (("O", "O"), ("H1", "H"), ("H2", "H")):
            lines.append(_water_line(serial, name, chain, resseq, x, element))
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.fixture(scope="module")
def large_water_pdb() -> str:
    return _water_pdb(N_WATERS)


@pytest.fixture(scope="module")
def reordered_pdb(large_water_pdb) -> str:
    return PDBService().reorder_and_inject_eps(large_water_pdb, None)


def _atom_lines(pdb: str) -> list[str]:
    return [line for line in pdb.splitlines() if line.startswith(("ATOM", "HETATM"))]


class TestReorderKeepsColumnsAboveSerial99999:
    def test_all_atoms_kept_with_fixed_columns(self, reordered_pdb):
        lines = _atom_lines(reordered_pdb)

        assert len(lines) == 3 * N_WATERS
        assert {line[17:20] for line in lines} == {"WAT"}
        assert {len(line) for line in lines} == {78}

    def test_serial_wraps_like_amber(self, reordered_pdb):
        lines = _atom_lines(reordered_pdb)

        assert lines[99_998][6:11] == "99999"
        assert lines[99_999][6:11] == "    0"
        assert lines[100_000][6:11] == "    1"

    def test_coordinates_survive_after_wrap(self, reordered_pdb):
        line = _atom_lines(reordered_pdb)[100_000]
        assert float(line[30:38]) == pytest.approx(33.3)
        assert float(line[38:46]) == pytest.approx(1.0)
        assert float(line[46:54]) == pytest.approx(2.0)

    def test_topology_dict_sees_every_water(self, reordered_pdb):
        mol = parse_pdb_to_topology_dict(reordered_pdb, {"W": "TIP3P"})

        waters = [res for res in mol["residues"] if res["mol_type"] == "W"]
        assert len(waters) == N_WATERS


class TestAmberCrd:
    def _crd_lines(self, pdb: str, tmp_path) -> list[str]:
        path = tmp_path / "structure.crd"
        # _generate_amber_crd nepotřebuje FF ani workspace - nekonstruujeme je
        TopologyService.__new__(TopologyService)._generate_amber_crd(pdb, path)
        return path.read_text(encoding="utf-8").splitlines()

    def test_has_title_and_all_atoms(self, reordered_pdb, tmp_path):
        lines = self._crd_lines(reordered_pdb, tmp_path)

        assert not lines[0].strip().isdigit()  # 1. řádek je titulek
        assert int(lines[1]) == 3 * N_WATERS
        assert len(lines) == 2 + (3 * 3 * N_WATERS + 5) // 6 + 1
        assert lines[-1].split() == ["108.2390000"] * 3 + ["90.0000000"] * 3

    def test_unparsable_atom_line_is_an_error(self, tmp_path):
        broken = (
            "HETATM100000  O   WAT Z2773      89.826  50.545  12.481  1.00  0.00           O\n"
        )
        with pytest.raises(ValueError, match="Cannot parse coordinates"):
            self._crd_lines(broken, tmp_path)


_BIOMT_IDENTITY_AND_SHIFT = (
    "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
    "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
    "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
    "REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000      200.00000",
    "REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000",
    "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000",
)


class TestProcessStructureAboveSerial99999:
    """Čištění altlocs / BIOMT (první krok po nahrání) přečíslovává atomy taky."""

    def _check(self, cleaned: str, expected_atoms: int) -> list[str]:
        lines = _atom_lines(cleaned)
        assert len(lines) == expected_atoms
        assert {line[17:20] for line in lines} == {"WAT"}
        assert {len(line) for line in lines} == {78}
        assert lines[99_999][6:11] == "    0"
        return lines

    def test_renumbering_without_symmetry(self, large_water_pdb):
        cleaned = process_structure(large_water_pdb, target_model=1, apply_symmetry=False, selection={})

        lines = self._check(cleaned, 3 * N_WATERS)
        assert float(lines[100_000][30:38]) == pytest.approx(33.3)

    def test_symmetry_expansion(self):
        half = N_WATERS // 2  # 2 BIOMT kopie -> 102 000 atomů
        pdb = _water_pdb(half, header=_BIOMT_IDENTITY_AND_SHIFT)

        cleaned = process_structure(pdb, target_model=1, apply_symmetry=True, selection={})

        lines = self._check(cleaned, 2 * 3 * half)
        # Atom 100 001 patří druhé (posunuté) kopii: voda 33333 - 17000 = 16333
        assert float(lines[100_000][30:38]) == pytest.approx(33.3 + 200.0)
        assert float(lines[100_000][38:46]) == pytest.approx(1.0)
