"""
Unit testy pro odstranění vod, které se po oříznutí vodní šablony do cílového
boxu překrývají přes periodickou hranici (viz
forge_molecule_solvation._remove_periodic_image_clashes). Dřív tam zůstávaly
dvojice vod s O-O pod 0.5 A - "chrchel vod na hraně boxu" u 2TRA.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "app" / "builder"))

from forge_molecule_solvation import (  # noqa: E402
    WaterGeometry,
    _remove_periodic_image_clashes,
)

pytestmark = pytest.mark.unit

_H = (np.array([0.9572, 0.0, 0.0]), np.array([-0.24, 0.9266, 0.0]))


def _water(x, y=0.0, z=0.0):
    return WaterGeometry(np.array([x, y, z], dtype=float), _H)


def test_waters_overlapping_through_boundary_are_deduplicated():
    vectors = np.diag([20.0, 20.0, 20.0])
    # x = -9.8 a +9.6 jsou přes periodickou hranici jen 0.6 A od sebe.
    waters = [_water(-9.8), _water(0.0), _water(9.6)]
    kept, removed = _remove_periodic_image_clashes(waters, vectors)
    assert removed == 1
    assert [w.oxygen[0] for w in kept] == [-9.8, 0.0]


def test_well_separated_boundary_waters_are_kept():
    vectors = np.diag([20.0, 20.0, 20.0])
    # 3.0 A přes hranici - normální kontakt, nic se nemaže.
    waters = [_water(-8.5), _water(8.5)]
    kept, removed = _remove_periodic_image_clashes(waters, vectors)
    assert removed == 0
    assert len(kept) == 2


def test_hydrogen_contact_through_boundary_counts_as_clash():
    vectors = np.diag([20.0, 20.0, 20.0])
    # O-O 2.6 A (nad limitem), ale vodíky obou vod míří proti sobě (H-H ~0.7 A).
    towards_minus_x = (np.array([-0.9572, 0.0, 0.0]), np.array([0.24, 0.9266, 0.0]))
    waters = [_water(8.7), WaterGeometry(np.array([-8.7, 0.0, 0.0]), towards_minus_x)]
    kept, removed = _remove_periodic_image_clashes(waters, vectors)
    assert removed == 1
    assert len(kept) == 1
