"""
Unit testy pro app/services/analysis_service.py::required_ff_groups.

Zjišťuje, jaké FORGE mol_type skupiny (P/R/D/W/I1/I1+/Im/Im+) daná struktura
reálně potřebuje, ať se dá zkontrolovat proti ff_selections ještě PŘED
spuštěním buildu (viz forge_service.ForgeStructureService._check_ff_coverage).

Regresní sada pro konverzaci: uživatel na reálném 1JJ2 vybral jen "I1+" a
narazil na "KeyError: Ion parameters missing for Im:Mg2+" po několika
minutách běhu - required_ff_groups musí přesně tenhle případ (Mg2+ krystalové
ionty vyžadují "Im", ne "I1+") odhalit staticky a rychle.
"""

import pytest

from app.services.analysis_service import required_ff_groups

pytestmark = pytest.mark.unit

_PROTEIN_PDB = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
)

_PROTEIN_WITH_MG_PDB = _PROTEIN_PDB + (
    "HETATM    3 MG    MG A 100      10.000  10.000  10.000  1.00  0.00          MG\n"
)


class TestRequiredFfGroups:
    def test_protein_only_no_solvent_requires_only_p(self):
        required = required_ff_groups(_PROTEIN_PDB, add_solvent_and_ions=False)
        assert set(required.keys()) == {"P"}

    def test_solvent_requested_adds_water_and_default_neutralization(self):
        required = required_ff_groups(_PROTEIN_PDB, add_solvent_and_ions=True)
        assert set(required.keys()) == {"P", "W", "I1"}
        assert "neutralization" in required["I1"]["reason"]

    def test_crystal_mg_ion_requires_im_not_i1_plus(self):
        """
        Přesně ten případ, co narazil na reálném 1JJ2: Mg2+ patří pod "Im",
        ne "I1+" (navzdory matoucímu názvu "Im", to není anion).
        """
        required = required_ff_groups(_PROTEIN_WITH_MG_PDB, add_solvent_and_ions=True)
        assert "Im" in required
        assert "Mg2+" in required["Im"]["reason"]
        assert "I1+" not in required

    def test_explicit_salt_mol_type_is_required_instead_of_default(self):
        salts = [{
            "cation": {"mol_type": "I1+", "resname": "Cs+"},
            "anion": {"mol_type": "I1+", "resname": "Br-"},
            "concentration": 0.15,
        }]
        required = required_ff_groups(_PROTEIN_PDB, add_solvent_and_ions=True, salts=salts)
        assert "I1+" in required
        assert "I1" not in required

    def test_no_solvent_no_ion_requirements_even_with_crystal_ions(self):
        required = required_ff_groups(_PROTEIN_WITH_MG_PDB, add_solvent_and_ions=False)
        assert set(required.keys()) == {"P"}
