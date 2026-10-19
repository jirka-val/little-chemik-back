"""
Rozpoznání už solvatovaného systému a kontroly pro jeho ponechání
(app/services/structure/solvated_system.py).
"""

import numpy as np
import pytest

from app.services.structure.solvated_system import (
    analyze_solvated_system, box_vectors, min_image_distance, restore_cryst1,
)

pytestmark = pytest.mark.unit


def _line(serial, name, resname, chain, resseq, xyz, element, record="ATOM"):
    return (f"{record:<6}{serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}")


def _system(box=30.0, spacing=3.1, solute_at=(15.0, 15.0, 15.0), drop_h=False, clash=False, ions=("NA", "CL")):
    """Kubický box: malý "solut" (2 rezidua), mřížka vod a ionty."""
    lines = [f"CRYST1{box:9.3f}{box:9.3f}{box:9.3f}  90.00  90.00  90.00 P 1           1"]
    serial = 1
    sx, sy, sz = solute_at
    for resseq, dx in ((1, 0.0), (2, 3.8)):
        for name, el, off in (("N", "N", (0, 0, 0)), ("CA", "C", (1.45, 0, 0)), ("C", "C", (2.4, 1.0, 0))):
            lines.append(_line(serial, name, "ALA", "A", resseq, (sx + dx + off[0], sy + off[1], sz + off[2]), el))
            serial += 1
    resseq = 1
    # mezi uzly mřížky vod (ta začíná na spacing/2)
    for ion, pos in zip(ions, ((spacing, spacing, spacing), (box - spacing, box - spacing, spacing))):
        lines.append(_line(serial, ion, ion, "I", resseq, pos, ion[:2] if ion != "K" else "K", "HETATM"))
        serial += 1
        resseq += 1
    grid = np.arange(spacing / 2, box, spacing)
    w = 1
    for x in grid:
        for y in grid:
            for z in grid:
                if np.linalg.norm(np.array([x, y, z]) - np.array(solute_at)) < 6.0:
                    continue
                atoms = [("O", "O", (x, y, z))]
                if not drop_h:
                    atoms += [("H1", "H", (x + 0.96, y, z)), ("H2", "H", (x - 0.24, y + 0.93, z))]
                for name, el, pos in atoms:
                    lines.append(_line(serial % 100000, name, "WAT", "W", w % 10000, pos, el, "HETATM"))
                    serial += 1
                w += 1
    if clash:
        lines.append(_line(serial, "O", "WAT", "X", 1, (sx + 0.3, sy, sz), "O", "HETATM"))
    return "\n".join(lines) + "\nEND\n"


def _statuses(result):
    return {c["id"]: c["status"] for c in result["checks"]}


def test_good_box_is_detected_and_passes():
    result = analyze_solvated_system(_system())
    assert result["solvated"] is True
    assert result["box"]["shape"] == "cubic"
    assert result["ions"] == {"CL": 1, "NA": 1}
    assert result["solute_residues"] == 2
    statuses = _statuses(result)
    assert statuses["clashes"] == "ok" and statuses["water_complete"] == "ok" and statuses["ions"] == "ok"
    assert statuses["image_distance"] == "ok"


def test_crystal_structure_is_not_solvated():
    # pár krystalových vod v buňce = hustota hluboko pod solvatovaným boxem
    pdb = _system(spacing=12.0)
    result = analyze_solvated_system(pdb)
    assert result["solvated"] is False and "checks" not in result


def test_small_box_and_problems_are_reported():
    result = analyze_solvated_system(_system(box=16.0, spacing=2.6, clash=True, drop_h=True, ions=("NA", "XX")))
    assert result["solvated"] is True
    statuses = _statuses(result)
    assert statuses["image_distance"] in ("warning", "error")
    assert statuses["clashes"] == "error"
    assert statuses["water_complete"] == "warning"   # bez vodíků - doplní se
    assert statuses["ions"] == "error"               # XX nemá FF skupinu


def test_min_image_distance_octahedron():
    box = {"a": 50.0, "b": 50.0, "c": 50.0, "alpha": 109.4712, "beta": 109.4712, "gamma": 109.4712}
    one_atom = np.array([[0.0, 0.0, 0.0]])
    # nejbližší obraz v oktaedru je ve vzdálenosti a (délka mřížkového vektoru)
    assert min_image_distance(one_atom, box_vectors(box)) == pytest.approx(50.0, abs=0.01)


def test_restore_cryst1():
    line = "CRYST1   51.551   51.551   51.551  70.53 109.47  70.53 P 1           1"
    assert restore_cryst1("ATOM ...\n", line).startswith("CRYST1")
    assert restore_cryst1(line + "\nATOM ...\n", line).count("CRYST1") == 1
