"""
Unit testy pro app/services/structure/forge_service.py - PDB writer.

Nejdůležitější pokrytí: bezpečnostní vrstva pro délku resname
(_pdb_safe_resname), kterou odhalil návrhový agent při plánování Fáze 4 -
format_pdb_atom_line (v app/builder) NEOŘEZÁVÁ resname na 3 znaky, jen ho
tiše rozšíří a posune všechny další sloupce na řádku. Testy tady staví
Molecule/Residue/Atom objekty přímo (bez spouštění run_forge_workflow),
takže jsou rychlé a bez závislosti na síti/silových polích.
"""

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

BUILDER_DIR = Path(__file__).resolve().parents[2] / "app" / "builder"
if str(BUILDER_DIR) not in sys.path:
    sys.path.insert(0, str(BUILDER_DIR))

from forge_molecule_parser import Molecule, Chain, Residue, Atom, PDBAtomRecord  # noqa: E402

from app.services.structure.forge_service import (  # noqa: E402
    ForgeWriterError,
    _pdb_safe_resname,
    build_forge_meta,
    molecule_to_pdb,
)


def _residue(chain_id, resseq, ff_resname, original_resname, group, atoms, terminus_reason=None):
    return Residue(
        chain_id=chain_id,
        resseq=resseq,
        icode="",
        ff_resname=ff_resname,
        atoms=atoms,
        index_in_chain=0,
        original_resname=original_resname,
        group=group,
        terminus_reason=terminus_reason,
    )


def _atom(name, coord=(1.0, 2.0, 3.0), element=None, occupancy=1.0, bfactor=0.0):
    return Atom(name=name, element=element, coord=coord, occupancy=occupancy, bfactor=bfactor)


class TestPdbSafeResname:
    def test_short_ff_resname_used_directly(self):
        residue = _residue("A", 1, "RU3", "U", "R", {})
        assert _pdb_safe_resname(residue) == "RU3"

    def test_long_ff_resname_falls_back_to_original(self):
        """CGLU (4 znaky, proteinový terminál) -> GLU (original_resname, builder ho nepřepisuje)."""
        residue = _residue("A", 83, "CGLU", "GLU", "P", {})
        assert _pdb_safe_resname(residue) == "GLU"

    def test_raises_when_no_safe_representation_exists(self):
        residue = _residue("A", 1, "UNREPRESENTABLE", "ALSOTOOLONG", "P", {})
        with pytest.raises(ForgeWriterError):
            _pdb_safe_resname(residue)


class TestMoleculeToPdb:
    def test_writes_valid_fixed_width_lines_for_polymer_residue(self):
        residue = _residue("A", 1, "RU5", "U", "R", {"P": _atom("P"), "O5'": _atom("O5'")})
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[residue])})

        pdb_text = molecule_to_pdb(molecule)
        atom_lines = [l for l in pdb_text.splitlines() if l.startswith(("ATOM", "HETATM"))]

        assert len(atom_lines) == 2
        for line in atom_lines:
            assert len(line) >= 66  # dost dlouhé pro occupancy/bfactor sloupce
            assert line[17:20].strip() == "RU5"
            assert line[21] == "A"
            assert int(line[22:26]) == 1

    def test_protein_terminal_residue_writes_original_resname_not_truncated(self):
        residue = _residue("A", 83, "CGLU", "GLU", "P", {"CA": _atom("CA")})
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[residue])})

        pdb_text = molecule_to_pdb(molecule)
        atom_line = next(l for l in pdb_text.splitlines() if l.startswith("ATOM"))
        assert atom_line[17:20].strip() == "GLU"
        # Kdyby se sem zapsalo CGLU beze zkrácení, sloupce by se posunuly a
        # tenhle sloupec by přestal odpovídat chain_id "A".
        assert atom_line[21] == "A"

    def test_ion_uses_element_symbol_for_both_atom_name_and_resname(self):
        residue = _residue(
            "I", 500, "Mg2+", "Mg2+", "Im",
            {"MG": _atom("MG", element="MG")},
        )
        molecule = Molecule(chains={"I": Chain(chain_id="I", residues=[residue])})

        pdb_text = molecule_to_pdb(molecule)
        atom_line = next(l for l in pdb_text.splitlines() if l.startswith("HETATM"))
        assert atom_line[17:20].strip() == "MG"
        assert atom_line[12:16].strip() == "MG"

    def test_water_residue_marked_as_hetatm_with_wat_resname(self):
        residue = _residue(
            "W", 1, "WAT", "WAT", "W3",
            {"O": _atom("O"), "H1": _atom("H1"), "H2": _atom("H2")},
        )
        molecule = Molecule(chains={"W": Chain(chain_id="W", residues=[residue])})

        pdb_text = molecule_to_pdb(molecule)
        atom_lines = [l for l in pdb_text.splitlines() if l.startswith("HETATM") and l[17:20].strip() == "WAT"]
        assert len(atom_lines) == 3

    def test_atoms_without_coordinates_are_skipped(self):
        """Nepostavené atomy (coord=None) se do výstupu nesmí dostat."""
        residue = _residue("A", 1, "ALA", "ALA", "P", {
            "CA": _atom("CA"),
            "CB": Atom(name="CB", element="C", coord=None),
        })
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[residue])})

        pdb_text = molecule_to_pdb(molecule)
        atom_lines = [l for l in pdb_text.splitlines() if l.startswith("ATOM")]
        assert len(atom_lines) == 1

    def test_passthrough_atoms_are_written(self):
        record = PDBAtomRecord(
            record_name="HETATM", serial=1, atom_name="C1", altloc="", resname="LIG",
            chain_id="L", resseq=1, icode="", coord=(0.0, 0.0, 0.0), group=None, ff_resname=None,
        )
        molecule = Molecule(chains={}, passthrough_atoms=[record])

        pdb_text = molecule_to_pdb(molecule)
        assert any(l.startswith("HETATM") and l[17:20].strip() == "LIG" for l in pdb_text.splitlines())

    def test_output_ends_with_end_record(self):
        molecule = Molecule(chains={})
        pdb_text = molecule_to_pdb(molecule)
        assert pdb_text.rstrip("\n").splitlines()[-1] == "END"


class TestGapBoundaryTer:
    """
    1JJ2 chain K: GLU83 (terminus_reason="gap", C-cap/OXT) je bezprostředně
    následováno PHE89 (terminus_reason="gap", N-cap/H1-H2-H3) - mezi nimi
    chybí residua 84-88, takže tam patří explicitní TER, i když číslo
    řetězce (chain_id) zůstává pro obě strany stejné "K".
    """

    def test_ter_inserted_between_artificial_gap_termini(self):
        before = _residue("K", 83, "CGLU", "GLU", "P", {"CA": _atom("CA"), "OXT": _atom("OXT")}, terminus_reason="gap")
        after = _residue("K", 89, "NPHE", "PHE", "P", {"N": _atom("N"), "H1": _atom("H1")}, terminus_reason="gap")
        # `tail` udrží PHE89 mimo pozici skutečného konce řetězce, ať se v
        # tomhle testu izolovaně ověří jen ten TER na hranici mezery (ten na
        # konci celého řetězce testuje samostatně
        # test_gap_boundary_ter_and_end_of_chain_ter_both_present níže).
        tail = _residue("K", 90, "ARG", "ARG", "P", {"CA": _atom("CA")})
        molecule = Molecule(chains={"K": Chain(chain_id="K", residues=[before, after, tail])})

        pdb_text = molecule_to_pdb(molecule)
        lines = pdb_text.splitlines()

        # Dva TERy celkem: jeden na hranici mezery (GLU83/PHE89), jeden na
        # skutečném konci řetězce po ARG90 - tenhle test se soustředí na ten první.
        ter_indices = [i for i, l in enumerate(lines) if l.startswith("TER")]
        assert len(ter_indices) == 2
        ter_line = lines[ter_indices[0]]
        assert ter_line[17:20].strip() == "GLU"
        assert int(ter_line[22:26]) == 83
        # TER musí sedět přesně mezi poslední atom GLU83 a první atom PHE89.
        assert lines[ter_indices[0] - 1].startswith("ATOM") and "GLU" in lines[ter_indices[0] - 1]
        assert lines[ter_indices[0] + 1].startswith("ATOM") and "PHE" in lines[ter_indices[0] + 1]

    def test_no_ter_between_ordinary_consecutive_residues(self):
        r1 = _residue("A", 1, "ALA", "ALA", "P", {"CA": _atom("CA")})
        r2 = _residue("A", 2, "ALA", "ALA", "P", {"CA": _atom("CA")})
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[r1, r2])})

        pdb_text = molecule_to_pdb(molecule)
        # Jediný TER je ten na konci celého řetězce (end-of-chain), žádný navíc uprostřed.
        assert pdb_text.count("TER") == 1

    def test_real_chain_end_still_gets_exactly_one_ter(self):
        """
        terminus_reason="chain_end" (skutečný začátek/konec řetězce) nesmí
        spustit tu samou "oba sousedé jsou umělý terminus" podmínku - jen
        jedno z dvojice sousedů ho kdy má, takže se žádný TER navíc
        uprostřed řetězce neobjeví, jen ten obvyklý na konci.
        """
        first = _residue("A", 1, "ALA", "ALA", "P", {"CA": _atom("CA")}, terminus_reason="chain_end")
        last = _residue("A", 2, "ALA", "ALA", "P", {"CA": _atom("CA")}, terminus_reason="chain_end")
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[first, last])})

        pdb_text = molecule_to_pdb(molecule)
        assert pdb_text.count("TER") == 1

    def test_gap_boundary_ter_and_end_of_chain_ter_both_present(self):
        before = _residue("K", 83, "CGLU", "GLU", "P", {"CA": _atom("CA")}, terminus_reason="gap")
        after = _residue("K", 89, "NPHE", "PHE", "P", {"CA": _atom("CA")}, terminus_reason="gap")
        tail = _residue("K", 92, "ASP", "ASP", "P", {"CA": _atom("CA")})
        molecule = Molecule(chains={"K": Chain(chain_id="K", residues=[before, after, tail])})

        pdb_text = molecule_to_pdb(molecule)
        assert pdb_text.count("TER") == 2


class TestBuildForgeMeta:
    def test_keys_by_chain_resseq_icode(self):
        residue = _residue("A", 83, "CGLU", "GLU", "P", {})
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[residue])})

        meta = build_forge_meta(molecule)
        assert meta == {"A:83:": {"ff_resname": "CGLU", "group": "P"}}
