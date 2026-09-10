"""
Unit testy pro app/services/analysis_service.py - AltLocs, modely, symetrie.

Pokrývá analyze_pdb_altlocs (detekce + doporučení trasy), clean_pdb_altlocs
(fyzická aplikace výběru) a process_structure (výběr modelu + AltLocs najednou).
"""

import pytest

from app.services.analysis_service import (
    analyze_pdb_altlocs,
    clean_pdb_altlocs,
    process_structure,
)

pytestmark = pytest.mark.unit


class TestAnalyzePdbAltlocs:
    def test_detects_altloc_residue(self, pdb_altloc_sample):
        result = analyze_pdb_altlocs(pdb_altloc_sample)
        assert result["hasAltLocs"] is True
        residues = result["residues"]
        assert len(residues) == 1
        assert residues[0]["chain"] == "A"
        assert residues[0]["resseq"] == 42
        assert set(residues[0]["altLocs"].keys()) == {"A", "B"}

    def test_recommends_higher_occupancy_variant(self, pdb_altloc_sample):
        result = analyze_pdb_altlocs(pdb_altloc_sample)
        residue = result["residues"][0]
        # Fixture ma altloc A s occupancy 0.6, B s 0.4 -> A ma vyhrat
        assert residue["recommended_alt"] == "A"

    def test_no_altlocs_reports_empty(self):
        pdb = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        result = analyze_pdb_altlocs(pdb)
        assert result["hasAltLocs"] is False
        assert result["residues"] == []

    def test_detects_multiple_models(self, pdb_two_model_nmr):
        result = analyze_pdb_altlocs(pdb_two_model_nmr)
        assert result["models"] == [1, 2]

    def test_detects_biomt_symmetry(self):
        pdb = (
            "REMARK 350 BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        )
        result = analyze_pdb_altlocs(pdb)
        assert result["hasSymmetry"] is True


class TestCleanPdbAltlocs:
    def test_removes_unselected_variant_and_resets_occupancy(self, pdb_altloc_sample):
        cleaned = clean_pdb_altlocs(pdb_altloc_sample, {"A_42_SER": "A"})

        og_lines = [l for l in cleaned.splitlines() if l[12:16].strip() == "OG"]
        assert len(og_lines) == 1
        kept = og_lines[0]
        assert kept[16] == " "  # altloc indikátor smazán
        assert kept[54:60].strip() == "1.00"  # occupancy vrácena na 1.00

    def test_unselected_residue_left_untouched(self, pdb_altloc_sample):
        """Pokud pro reziduum není v selection žádná volba, řádky se nesmí smazat."""
        cleaned = clean_pdb_altlocs(pdb_altloc_sample, {})
        og_lines = [l for l in cleaned.splitlines() if l[12:16].strip() == "OG"]
        assert len(og_lines) == 2  # obě varianty zůstaly


class TestProcessStructure:
    def test_selects_only_requested_model(self, pdb_two_model_nmr):
        cleaned = process_structure(pdb_two_model_nmr, target_model=2, apply_symmetry=False, selection={})
        atom_lines = [l for l in cleaned.splitlines() if l.startswith("ATOM")]
        assert len(atom_lines) == 5  # jen atomy z MODEL 2
        # Model 2 ma x posunuty o +5.0 (viz fixture generator)
        first_x = float(atom_lines[0][30:38])
        assert first_x == pytest.approx(6.0, abs=0.01)

    def test_applies_altloc_selection_during_model_extraction(self, pdb_altloc_sample):
        cleaned = process_structure(pdb_altloc_sample, target_model=1, apply_symmetry=False,
                                     selection={"A_42_SER": "B"})
        og_lines = [l for l in cleaned.splitlines() if l[12:16].strip() == "OG"]
        assert len(og_lines) == 1
        assert og_lines[0][54:60].strip() == "1.00"

    def test_identity_biomt_when_no_symmetry_present(self, pdb_alanine_single):
        """Bez REMARK 350 se pouzije identita - vystup by mel mit stejny pocet atomu."""
        cleaned = process_structure(pdb_alanine_single, target_model=1, apply_symmetry=True, selection={})
        input_atoms = sum(1 for l in pdb_alanine_single.splitlines() if l.startswith("ATOM"))
        output_atoms = sum(1 for l in cleaned.splitlines() if l.startswith("ATOM"))
        assert output_atoms == input_atoms
