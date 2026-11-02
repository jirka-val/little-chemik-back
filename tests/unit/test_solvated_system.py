"""
Rozpoznání už solvatovaného systému a kontroly pro jeho ponechání
(app/services/structure/solvated_system.py).
"""

import numpy as np
import pytest

from app.services.structure.solvated_system import (
    analyze_solvated_system, box_vectors, min_image_distance, restore_cryst1, standardize_reconstructed_box,
    strip_solvent,
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


# ---------------------------------------------------------------------------
# Bez CRYST1: box se rekonstruuje z vody (box_reconstruction.py)
# ---------------------------------------------------------------------------

def _without_cryst1(pdb, rotation=None, shift=(0.0, 0.0, 0.0)):
    """Smaže CRYST1 a souřadnice otočí/posune (jako export z jiného programu)."""
    out = []
    for line in pdb.splitlines():
        if line.startswith("CRYST1"):
            continue
        if line.startswith(("ATOM", "HETATM")) and rotation is not None:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ rotation.T + shift
            line = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def _rotation(deg_z, deg_x):
    cz, sz = np.cos(np.radians(deg_z)), np.sin(np.radians(deg_z))
    cx, sx = np.cos(np.radians(deg_x)), np.sin(np.radians(deg_x))
    return np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]]) @ np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])


def test_missing_cryst1_box_is_reconstructed():
    pdb = _without_cryst1(_system(box=34.0, solute_at=(17.0, 17.0, 17.0)), _rotation(27, 41), (5.0, -12.0, 30.0))
    result = analyze_solvated_system(pdb)
    assert result["solvated"] is True
    assert result["box_source"] == "reconstructed"
    assert result["box"]["shape"] == "cubic"
    assert result["box"]["a"] == pytest.approx(34.0, abs=0.6)
    assert result["checks"][0]["id"] == "box_source" and result["checks"][0]["status"] == "warning"
    statuses = _statuses(result)
    # kontroly běží ve standardní orientaci rekonstruovaného boxu (stejně jako s CRYST1)
    with_box = _statuses(analyze_solvated_system(_system(box=34.0, solute_at=(17.0, 17.0, 17.0))))
    assert {k: v for k, v in statuses.items() if k != "box_source"} == with_box
    assert result["box_reconstruction"]["shape"] == "cubic"


def test_quick_detection_skips_reconstruction():
    pdb = _without_cryst1(_system(box=34.0, solute_at=(17.0, 17.0, 17.0)))
    result = analyze_solvated_system(pdb, with_checks=False, reconstruct=False)
    assert result["solvated"] is True and result["box_source"] == "missing" and result["box"] is None


def test_standardize_reconstructed_box():
    pdb = _without_cryst1(_system(box=34.0, solute_at=(17.0, 17.0, 17.0)), _rotation(27, 41), (5.0, -12.0, 30.0))
    text, cryst1 = standardize_reconstructed_box(pdb)
    assert text.startswith(cryst1) and cryst1.startswith("CRYST1")
    again = analyze_solvated_system(text)
    assert again["box_source"] == "cryst1" and again["box"]["shape"] == "cubic"
    waters = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in text.splitlines() if l.startswith("HETATM") and l[17:20] == "WAT" and l[12:16].strip() == "O"])
    a = again["box"]["a"]
    assert waters.min() > -1.0 and waters.max() < a + 1.0
    assert standardize_reconstructed_box(_system(box=34.0)) is None      # CRYST1 už je


def test_liquid_water_without_a_box_shape():
    # kulová kapka vody: kapalina, ale žádný periodický box -> "keep" nejde
    lines, serial = [], 1
    grid = np.arange(-20.0, 20.1, 3.1)
    for x in grid:
        for y in grid:
            for z in grid:
                if np.linalg.norm([x, y, z]) < 19.0:
                    lines.append(_line(serial, "O", "WAT", "W", serial % 10000, (x, y, z), "O", "HETATM"))
                    serial += 1
    result = analyze_solvated_system("\n".join(lines) + "\nEND\n")
    assert result["solvated"] is True and result["box_source"] == "missing" and result["box"] is None
    assert [c["status"] for c in result["checks"]] == ["error"]


def test_crystal_waters_without_cryst1_are_not_solvated():
    pdb = _without_cryst1(_system(spacing=12.0))
    result = analyze_solvated_system(pdb)
    assert result["solvated"] is False and result["box_source"] is None


def test_strip_solvent_removes_water_monovalent_ions_and_box():
    mg = _line(99998, "MG", "MG", "M", 1, (4.0, 26.0, 4.0), "MG", "HETATM")
    pdb = _system().replace("END\n", mg + "\nTER\nCONECT    1    2\nCONECT 9999 9998\nEND\n")
    stripped, removed = strip_solvent(pdb)
    assert removed["ions"] == {"CL": 1, "NA": 1}
    assert removed["waters"] == analyze_solvated_system(pdb)["waters"]
    assert "CRYST1" not in stripped and " WAT " not in stripped
    assert " MG " in stripped                               # multivalentní iont zůstává
    assert sum(1 for l in stripped.splitlines() if l.startswith("ATOM")) == 6
    result = analyze_solvated_system(stripped)
    assert result["solvated"] is False and result["waters"] == 0


def test_strip_solvent_keeps_a_structure_without_solvent():
    lines = [_line(i + 1, name, "ALA", "A", 1, (i * 1.5, 0.0, 0.0), el) for i, (name, el) in
             enumerate((("N", "N"), ("CA", "C"), ("C", "C")))]
    pdb = "\n".join(lines) + "\nTER\nEND\n"
    stripped, removed = strip_solvent(pdb)
    assert removed == {"waters": 0, "ions": {}}
    assert stripped == pdb


def test_restore_cryst1():
    line = "CRYST1   51.551   51.551   51.551  70.53 109.47  70.53 P 1           1"
    assert restore_cryst1("ATOM ...\n", line).startswith("CRYST1")
    assert restore_cryst1(line + "\nATOM ...\n", line).count("CRYST1") == 1
