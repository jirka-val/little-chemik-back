"""
Unit testy pro app/services/pdb_service.py::parse_pdb_to_topology_dict.

Klíčové pokrytí: chování BEZE ZMĚNY, když forge_meta chybí (regrese proti
původnímu, čistě heuristickému chování), a korektní přednost forge_meta,
když je k dispozici - to je jádro opravy z Fáze 6 (proteinové terminály
CGLU/NPHE na okraji sekvenční díry, které se do 3sloupcového PDB nevejdou
a bez sidecaru by AMBER_topology mohl přes díru tiše vytvořit vazbu).
"""

import pytest

from app.services.pdb_service import parse_pdb_to_topology_dict

pytestmark = pytest.mark.unit


GAP_PDB = (
    "ATOM      1  CA  GLU A  83       1.000   1.000   1.000  1.00  0.00           C\n"
    "ATOM      2  CA  PHE A  89       2.000   2.000   2.000  1.00  0.00           C\n"
)


class TestWithoutSidecar:
    """Dokumentuje dnešní (pre-existující, mimo rozsah Fáze 6) chování beze změny."""

    def test_protein_defaults_to_mol_type_r(self):
        """
        Známý, zdokumentovaný nedostatek: proteiny bez sidecaru dostanou
        mol_type 'R' (výchozí), ne 'P' - terminální/gap logika pro protein
        v tomhle heuristickém režimu vůbec neexistuje.
        """
        result = parse_pdb_to_topology_dict(GAP_PDB)
        assert all(r["mol_type"] == "R" for r in result["residues"])

    def test_resnames_unchanged_across_gap(self):
        result = parse_pdb_to_topology_dict(GAP_PDB)
        names = [r["resn"] for r in result["residues"]]
        assert names == ["GLU", "PHE"]

    def test_water_and_ion_detection_still_works(self):
        pdb = (
            "ATOM      1  O   HOH A 200       0.000   0.000   0.000  1.00  0.00           O\n"
            "ATOM      2 NA   NA  A 201       1.000   1.000   1.000  1.00  0.00          NA\n"
        )
        result = parse_pdb_to_topology_dict(pdb)
        by_name = {r["resn"]: r["mol_type"] for r in result["residues"]}
        assert by_name["HOH"] == "W"
        assert by_name["NA"] == "I"


class TestWithSidecar:
    """Fáze 6: forge_meta má přednost před heuristikou."""

    def test_uses_authoritative_terminal_resname(self):
        forge_meta = {
            "A:83:": {"ff_resname": "CGLU", "group": "P"},
            "A:89:": {"ff_resname": "NPHE", "group": "P"},
        }
        result = parse_pdb_to_topology_dict(GAP_PDB, forge_meta=forge_meta)
        residues = {(r["chain"], r["resn"]): r["mol_type"] for r in result["residues"]}
        assert residues[("A", "CGLU")] == "P"
        assert residues[("A", "NPHE")] == "P"

    def test_partial_sidecar_falls_back_per_residue(self):
        """
        Sidecar nemusí pokrývat úplně všechno (např. workspace prošel starším
        /prepare) - rezidua bez záznamu spadnou zpět na heuristiku.
        """
        forge_meta = {"A:83:": {"ff_resname": "CGLU", "group": "P"}}
        result = parse_pdb_to_topology_dict(GAP_PDB, forge_meta=forge_meta)
        by_resseq = {}
        for line, r in zip(GAP_PDB.splitlines(), result["residues"]):
            by_resseq[int(line[22:26])] = r
        assert by_resseq[83]["resn"] == "CGLU"
        assert by_resseq[83]["mol_type"] == "P"
        assert by_resseq[89]["resn"] == "PHE"  # heuristika, beze změny
        assert by_resseq[89]["mol_type"] == "R"  # stejný pre-existující default jako bez sidecaru

    def test_forge_group_to_mol_type_mapping(self):
        pdb = (
            "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
            "ATOM      2 MG   MG  A   2       1.000   1.000   1.000  1.00  0.00          MG\n"
        )
        forge_meta = {
            "A:1:": {"ff_resname": "WAT", "group": "W3"},
            "A:2:": {"ff_resname": "MG", "group": "Im"},
        }
        result = parse_pdb_to_topology_dict(pdb, forge_meta=forge_meta)
        by_resn = {r["resn"]: r["mol_type"] for r in result["residues"]}
        assert by_resn["WAT"] == "W"
        assert by_resn["MG"] == "I"
