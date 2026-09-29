"""
Unit testy pro app/services/hmr.py - Hydrogen Mass Repartitioning nad
AMBER topologickým slovníkem. Topologie se staví ručně (jen pole, která
apply_hmr čte), bez silových polí a builderu.
"""

import pytest

from app.services.hmr import DEFAULT_HMR_H_MASS, apply_hmr

pytestmark = pytest.mark.unit


def _methane_and_water():
    # Atomy: 0 C, 1-4 H (methan, residue "MET"), 5 O, 6-7 H (voda "WAT")
    return {
        "ATOM_NAME": ["C", "H1", "H2", "H3", "H4", "O", "H1", "H2"],
        "ATOMIC_NUMBER": [6, 1, 1, 1, 1, 8, 1, 1],
        "MASS": [12.01, 1.008, 1.008, 1.008, 1.008, 16.0, 1.008, 1.008],
        "RESIDUE_LABEL": ["MET", "WAT"],
        "RESIDUE_POINTER": [1, 6],
        # AMBER indexy = 3 * atom; třetí číslo je index parametru vazby.
        "BONDS_INC_HYDROGEN": [0, 3, 1, 0, 6, 1, 0, 9, 1, 0, 12, 1,
                               15, 18, 2, 15, 21, 2, 18, 21, 3],
    }


def test_repartitions_solute_hydrogens_and_keeps_total_mass():
    top = _methane_and_water()
    total_before = sum(top["MASS"])

    count = apply_hmr(top)

    assert count == 4
    assert top["MASS"][1:5] == [DEFAULT_HMR_H_MASS] * 4
    assert top["MASS"][0] == pytest.approx(12.01 - 4 * (DEFAULT_HMR_H_MASS - 1.008))
    assert sum(top["MASS"]) == pytest.approx(total_before)


def test_water_is_left_untouched():
    top = _methane_and_water()
    apply_hmr(top)
    assert top["MASS"][5:] == [16.0, 1.008, 1.008]


def test_bond_order_does_not_matter():
    top = _methane_and_water()
    # H první, C druhý
    top["BONDS_INC_HYDROGEN"][0:3] = [3, 0, 1]
    apply_hmr(top)
    assert top["MASS"][1] == DEFAULT_HMR_H_MASS


def test_non_positive_heavy_mass_is_rejected():
    top = _methane_and_water()
    top["MASS"][0] = 5.0  # nesmyslně lehký těžký atom
    with pytest.raises(ValueError):
        apply_hmr(top)
