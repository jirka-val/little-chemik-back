"""
Rekonstrukce periodického boxu z vody (app/services/structure/box_reconstruction.py):
syntetické periodické vodní boxy v náhodné orientaci, CRYST1 se "zahodí".
"""

import itertools
import math

import numpy as np
import pytest

from app.services.structure.box_reconstruction import reconstruct_box, standard_vectors

pytestmark = pytest.mark.unit


def _rotation(rng):
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def periodic_water(lattice, wigner_seitz, seed, cavity=10.0, spacing=3.1, jitter=0.45):
    """Kyslíky periodické vody (jitterovaná frakční mřížka - perioda sedí přesně),
    dutina pro solut, náhodná rotace a posun. Vrací i rotaci a posun."""
    rng = np.random.default_rng(seed)
    n = [max(2, round(np.linalg.norm(a) / spacing)) for a in lattice]
    frac = np.array(list(itertools.product(*[(np.arange(k) + 0.5) / k - 0.5 for k in n])))
    X = frac @ lattice + rng.normal(scale=jitter, size=(len(frac), 3))
    if wigner_seitz:   # tvar tělesa (tleap), ne rovnoběžnostěn
        shifts = np.array(list(itertools.product(range(-2, 3), repeat=3))) @ lattice
        X = X - shifts[np.linalg.norm(X[:, None, :] - shifts[None, :, :], axis=2).argmin(axis=1)]
    X = X[np.linalg.norm(X, axis=1) > cavity]
    R, t = _rotation(rng), rng.uniform(-30, 30, size=3)
    return X @ R.T + t


def _params(shape, size):
    angles = {"truncated octahedron": (109.4712,) * 3, "rhombic dodecahedron": (60.0, 60.0, 90.0)}[shape]
    return dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), (size, size, size, *angles)))


@pytest.mark.parametrize("name, lattice, wigner_seitz, shape", [
    ("cube", np.diag([36.0, 36.0, 36.0]), False, "cubic"),
    ("cuboid", np.diag([32.0, 38.0, 44.0]), False, "rectangular"),
    ("truncated octahedron (tleap shape)", standard_vectors(_params("truncated octahedron", 46.0)), True,
     "truncated octahedron"),
    ("truncated octahedron as a triclinic cell", standard_vectors(_params("truncated octahedron", 46.0)), False,
     "truncated octahedron"),
    ("rhombic dodecahedron", standard_vectors(_params("rhombic dodecahedron", 42.0)), True, "rhombic dodecahedron"),
])
def test_shape_size_and_orientation_are_recovered(name, lattice, wigner_seitz, shape):
    oxygens = periodic_water(lattice, wigner_seitz, seed=3)
    fit = reconstruct_box(oxygens)
    assert fit.ok, fit.reason
    assert fit.shape == shape
    true_volume = abs(np.linalg.det(lattice))
    # ~0.4 A na hranu; chyba jde spíš nahoru (díra se v NPT zavře, překryv je horší)
    assert fit.volume == pytest.approx(true_volume, rel=0.05)
    # v nalezené orientaci je voda periodická: přes hranu nejsou překryvy ani díry navíc
    assert fit.clash_fraction <= fit.interior_clash_fraction + max(0.05, 0.25 * fit.interior_clash_fraction)
    assert fit.gap_fraction <= 0.05


def test_cuboid_edges_are_recovered_individually():
    fit = reconstruct_box(periodic_water(np.diag([32.0, 38.0, 44.0]), False, seed=5))
    assert fit.ok and fit.shape == "rectangular"
    assert sorted(fit.params[k] for k in "abc") == pytest.approx([32.0, 38.0, 44.0], abs=0.6)


def test_transform_puts_water_into_the_standard_box():
    lattice = np.diag([36.0, 36.0, 36.0])
    oxygens = periodic_water(lattice, False, seed=7)
    fit = reconstruct_box(oxygens)
    frac = fit.transform(oxygens) @ np.linalg.inv(standard_vectors(fit.params))
    assert frac.min() > -0.03 and frac.max() < 1.03


def test_water_droplet_is_not_a_box():
    rng = np.random.default_rng(1)
    points = rng.uniform(-30, 30, size=(30000, 3))
    droplet = points[np.linalg.norm(points, axis=1) < 26][:2400]
    fit = reconstruct_box(droplet)
    assert not fit.ok and "periodic cell" in fit.reason


def test_smeared_trajectory_frame_is_rejected():
    rng = np.random.default_rng(2)
    oxygens = periodic_water(standard_vectors(_params("truncated octahedron", 46.0)), True, seed=2)
    moved = rng.random(len(oxygens)) < 0.3
    oxygens[moved] += rng.normal(scale=5.0, size=(moved.sum(), 3))
    assert not reconstruct_box(oxygens).ok


def test_too_few_waters():
    fit = reconstruct_box(np.zeros((50, 3)))
    assert not fit.ok and "too few" in fit.reason
    assert math.isfinite(fit.volume)
