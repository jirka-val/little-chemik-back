"""
Rozpozná, že nahrané PDB je už hotový solvatovaný systém (box v CRYST1,
voda kolem solutu, ionty), a zkontroluje, jestli se dá použít tak, jak je.

Uživatel pak volí (PreparationRequest.existing_solvent):
- "keep"    nechat vodu, ionty i box, nic nepřidávat (jen dostavět, co chybí)
- "rebuild" vodu a monovalentní ionty odstranit a solvatovat znovu

Krystalová struktura z RCSB má CRYST1 taky (jednotková buňka) a pár desítek
krystalových vod - od solvatovaného boxu se liší hustotou vody v buňce.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from app.services.pdb_service import ion_group_for_resname

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP", "DOD"}

# Čistá voda při 298 K: 0.997 g/cm3 = 0.0334 molekul na A^3.
WATER_NUMBER_DENSITY = 0.0334
# Parciální specifický objem proteinu ~0.73 cm3/g = 1.21 A^3 na Da.
SOLUTE_VOLUME_PER_DALTON = 1.21

# Solvatovaný box má v buňce řádově 0.02-0.03 vod/A^3, krystal ~0.001.
SOLVATED_MIN_WATERS = 50
SOLVATED_MIN_WATER_DENSITY = 0.01

MIN_IMAGE_DISTANCE_OK = 10.0     # A - aspoň nonbonded cutoff
MIN_IMAGE_DISTANCE_SPLIT = 3.0   # pod tím solut nejspíš přesahuje hranu boxu
CLASH_DISTANCE = 1.2             # A - těžké atomy různých reziduí blíž = překryv
_IMAGE_SAMPLE_ATOMS = 3000

_MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "S": 32.06,
           "MG": 24.305, "ZN": 65.38, "CA": 40.078, "NA": 22.99, "K": 39.098, "CL": 35.45,
           "F": 18.998, "BR": 79.904, "I": 126.9, "FE": 55.845, "MN": 54.938, "CU": 63.546,
           "SE": 78.971, "LI": 6.94, "RB": 85.468, "CS": 132.91}


@dataclass
class _Atom:
    name: str
    resname: str
    res_key: tuple
    element: str
    xyz: tuple


@dataclass
class _Parsed:
    box: Optional[Dict[str, float]] = None
    atoms: List[_Atom] = field(default_factory=list)


def _element(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    name = line[12:16].strip().upper()
    return name.lstrip("0123456789")[:1]


def _parse(pdb_text: str) -> _Parsed:
    parsed = _Parsed()
    for line in pdb_text.splitlines():
        if line.startswith("CRYST1"):
            try:
                a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
                alpha, beta, gamma = float(line[33:40]), float(line[40:47]), float(line[47:54])
            except ValueError:
                continue
            # 1x1x1 A je zástupná buňka (NMR, modely) - žádný box
            if min(a, b, c) > 2.0:
                parsed.box = {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}
        elif line.startswith(("ATOM", "HETATM")):
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            parsed.atoms.append(_Atom(
                name=line[12:16].strip(),
                resname=line[17:21].strip(),
                res_key=(line[21], line[22:27]),
                element=_element(line),
                xyz=xyz,
            ))
    return parsed


def box_vectors(box: Dict[str, float]) -> np.ndarray:
    """Mřížkové vektory z CRYST1 (a podél x, b v rovině xy)."""
    a, b, c = box["a"], box["b"], box["c"]
    alpha, beta, gamma = (math.radians(box[k]) for k in ("alpha", "beta", "gamma"))
    bx, by = b * math.cos(gamma), b * math.sin(gamma)
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]])


def box_shape(box: Dict[str, float]) -> str:
    angles = (box["alpha"], box["beta"], box["gamma"])
    if all(abs(x - 90.0) < 0.1 for x in angles):
        return "cubic" if abs(box["a"] - box["b"]) < 0.01 and abs(box["a"] - box["c"]) < 0.01 else "rectangular"
    if all(abs(x - 109.4712) < 0.1 for x in angles) or (
            abs(angles[0] - 70.5288) < 0.1 and abs(angles[1] - 109.4712) < 0.1 and abs(angles[2] - 70.5288) < 0.1):
        return "truncated octahedron"
    return "other"


def _residues(atoms: List[_Atom]) -> Dict[tuple, List[_Atom]]:
    residues: Dict[tuple, List[_Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.res_key, []).append(atom)
    return residues


def _ion_charge(resname: str) -> Optional[int]:
    """Náboj iontu ze jména (Na+, Cl-, Mg2+, MG, ZN ...), None když ho nepoznáme."""
    name = resname.strip()
    if name.endswith(("+", "-")):
        sign = 1 if name.endswith("+") else -1
        digits = "".join(ch for ch in name if ch.isdigit())
        return sign * (int(digits) if digits else 1)
    return {"NA": 1, "K": 1, "LI": 1, "RB": 1, "CS": 1, "CL": -1, "BR": -1, "F": -1, "I": -1,
            "MG": 2, "CA": 2, "ZN": 2, "MN": 2, "CO": 2, "NI": 2, "CU": 2, "FE": 2}.get(name.upper())


# ---------------------------------------------------------------------------
# Geometrie
# ---------------------------------------------------------------------------

def _close_pairs(coords: np.ndarray, cutoff: float) -> np.ndarray:
    """Dvojice indexů (i < j) atomů blíž než cutoff - mřížka buněk o hraně cutoff."""
    if len(coords) < 2:
        return np.empty((0, 2), dtype=np.int64)
    cells = np.floor((coords - coords.min(axis=0)) / cutoff).astype(np.int64)
    dims = cells.max(axis=0) + 3
    keys = ((cells[:, 0] + 1) * dims[1] + (cells[:, 1] + 1)) * dims[2] + (cells[:, 2] + 1)
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    pairs = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                shift = (dx * dims[1] + dy) * dims[2] + dz
                lo = np.searchsorted(sorted_keys, keys + shift, side="left")
                hi = np.searchsorted(sorted_keys, keys + shift, side="right")
                counts = hi - lo
                if not counts.any():
                    continue
                i = np.repeat(np.arange(len(coords)), counts)
                starts = np.repeat(lo, counts)
                offsets = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
                j = order[starts + offsets]
                keep = i < j
                i, j = i[keep], j[keep]
                d = np.linalg.norm(coords[i] - coords[j], axis=1)
                hit = d < cutoff
                pairs.append(np.stack([i[hit], j[hit]], axis=1))
    return np.concatenate(pairs) if pairs else np.empty((0, 2), dtype=np.int64)


def min_image_distance(solute: np.ndarray, vectors: np.ndarray) -> float:
    """Nejmenší vzdálenost solutu od jeho periodického obrazu (vzorek těžkých atomů)."""
    if len(solute) == 0:
        return math.inf
    if len(solute) > _IMAGE_SAMPLE_ATOMS:
        solute = solute[np.linspace(0, len(solute) - 1, _IMAGE_SAMPLE_ATOMS).astype(int)]
    points = solute.astype(np.float32)
    best = math.inf
    seen = set()
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if (i, j, k) == (0, 0, 0) or (-i, -j, -k) in seen:
                    continue
                seen.add((i, j, k))
                shift = (i * vectors[0] + j * vectors[1] + k * vectors[2]).astype(np.float32)
                shifted = points + shift
                for start in range(0, len(points), 500):
                    chunk = points[start:start + 500]
                    d2 = ((chunk[:, None, :] - shifted[None, :, :]) ** 2).sum(axis=2)
                    best = min(best, float(np.sqrt(d2.min())))
    return best


# ---------------------------------------------------------------------------
# Detekce + kontroly
# ---------------------------------------------------------------------------

def _check(check_id: str, label: str, status: str, detail: str) -> Dict[str, str]:
    return {"id": check_id, "label": label, "status": status, "detail": detail}


def analyze_solvated_system(pdb_text: str, with_checks: bool = True) -> Dict[str, Any]:
    """
    {"solvated": bool, ...}. Když je systém solvatovaný, obsahuje i box,
    počty a seznam kontrol (status ok/warning/error) pro volbu "keep".
    """
    parsed = _parse(pdb_text)
    residues = _residues(parsed.atoms)
    waters = {k: v for k, v in residues.items() if v[0].resname in WATER_RESNAMES}
    ions = {k: v for k, v in residues.items() if k not in waters and len(v) == 1 and v[0].element != "H"}
    solute_atoms = [a for k, v in residues.items() if k not in waters and k not in ions for a in v]

    box = parsed.box
    volume = float(abs(np.linalg.det(box_vectors(box)))) if box else None
    water_density = len(waters) / volume if volume else 0.0
    solvated = bool(box) and len(waters) >= SOLVATED_MIN_WATERS and water_density >= SOLVATED_MIN_WATER_DENSITY

    ion_counts: Dict[str, int] = {}
    for atoms in ions.values():
        ion_counts[atoms[0].resname] = ion_counts.get(atoms[0].resname, 0) + 1

    result: Dict[str, Any] = {
        "solvated": solvated,
        "box": {**box, "shape": box_shape(box), "volume": round(volume, 1)} if box else None,
        "waters": len(waters),
        "ions": dict(sorted(ion_counts.items())),
        "solute_atoms": len(solute_atoms),
        "solute_residues": len({a.res_key for a in solute_atoms}),
    }
    if solvated and with_checks:
        result["checks"] = _checks(box, volume, waters, ion_counts, solute_atoms, parsed.atoms)
    return result


def _checks(box, volume, waters, ion_counts, solute_atoms, all_atoms) -> List[Dict[str, str]]:
    checks = []
    shape = box_shape(box)

    # 1. tvar boxu
    if shape == "other":
        checks.append(_check("box_shape", "Box shape", "warning",
                             f"Unusual cell angles {box['alpha']:g}/{box['beta']:g}/{box['gamma']:g} deg."))
    else:
        checks.append(_check("box_shape", "Box shape", "ok",
                             f"{shape}, {box['a']:.1f} x {box['b']:.1f} x {box['c']:.1f} A"))

    # 2. kompletní molekuly vody
    incomplete = sum(1 for atoms in waters.values()
                     if sum(1 for a in atoms if a.element == "O") != 1
                     or sum(1 for a in atoms if a.element == "H") not in (0, 2))
    no_h = sum(1 for atoms in waters.values() if not any(a.element == "H" for a in atoms))
    if incomplete:
        checks.append(_check("water_complete", "Water molecules", "error",
                             f"{incomplete} of {len(waters)} water molecules have missing or extra atoms."))
    elif no_h:
        checks.append(_check("water_complete", "Water molecules", "warning",
                             f"{no_h} water molecules have no hydrogens - they will be added."))
    else:
        checks.append(_check("water_complete", "Water molecules", "ok", f"{len(waters)} complete molecules"))

    # 3. hustota vody ve volném objemu boxu
    solute_mass = sum(_MASSES.get(a.element, 12.0) for a in solute_atoms)
    free_volume = max(volume - solute_mass * SOLUTE_VOLUME_PER_DALTON, 1.0)
    ratio = len(waters) / free_volume / WATER_NUMBER_DENSITY
    text = f"{ratio * 100:.0f} % of liquid water density in the free volume of the box"
    if ratio < 0.6 or ratio > 1.2:
        checks.append(_check("water_density", "Water density", "error", text + " - the box is not filled properly."))
    elif ratio < 0.75 or ratio > 1.1:
        checks.append(_check("water_density", "Water density", "warning",
                             text + " - expect the box to change noticeably during NPT equilibration."))
    else:
        checks.append(_check("water_density", "Water density", "ok", text))

    # 4. vzdálenost solutu od jeho periodického obrazu
    solute_heavy = np.array([a.xyz for a in solute_atoms if a.element != "H"], dtype=float)
    if len(solute_heavy):
        distance = min_image_distance(solute_heavy, box_vectors(box))
        if distance < MIN_IMAGE_DISTANCE_SPLIT:
            checks.append(_check("image_distance", "Distance to periodic image", "error",
                                 f"{distance:.1f} A - the solute crosses the box boundary or touches its own "
                                 "image (it may need re-imaging/centering)."))
        elif distance < MIN_IMAGE_DISTANCE_OK:
            checks.append(_check("image_distance", "Distance to periodic image", "warning",
                                 f"{distance:.1f} A - shorter than a typical 10 A cutoff; the solute can "
                                 "interact with its own image."))
        else:
            checks.append(_check("image_distance", "Distance to periodic image", "ok", f"{distance:.1f} A"))

    # 5. překryvy těžkých atomů různých reziduí
    heavy = [a for a in all_atoms if a.element != "H" and not a.name.upper().startswith("EP")]
    coords = np.array([a.xyz for a in heavy], dtype=float)
    pairs = _close_pairs(coords, CLASH_DISTANCE)
    clashes = [(heavy[i], heavy[j]) for i, j in pairs if heavy[i].res_key != heavy[j].res_key]
    if clashes:
        a, b = clashes[0]
        checks.append(_check("clashes", "Overlapping atoms", "error",
                             f"{len(clashes)} pair(s) closer than {CLASH_DISTANCE} A, e.g. "
                             f"{a.resname} {a.res_key[0]}{a.res_key[1].strip()} {a.name} - "
                             f"{b.resname} {b.res_key[0]}{b.res_key[1].strip()} {b.name}."))
    else:
        checks.append(_check("clashes", "Overlapping atoms", "ok", "none"))

    # 6. ionty - známé FF a jejich náboj
    if ion_counts:
        unknown = [name for name in ion_counts if ion_group_for_resname(name) is None]
        charges = [(_ion_charge(name), n) for name, n in ion_counts.items()]
        ion_charge = sum(q * n for q, n in charges if q is not None)
        listing = ", ".join(f"{name} {n}" for name, n in sorted(ion_counts.items()))
        if unknown:
            checks.append(_check("ions", "Ions", "error",
                                 f"No force field group for: {', '.join(sorted(unknown))} ({listing})."))
        else:
            checks.append(_check("ions", "Ions", "ok",
                                 f"{listing}; total ion charge {ion_charge:+d} e. The net charge of the whole "
                                 "system is checked when the topology is generated."))
    else:
        checks.append(_check("ions", "Ions", "warning",
                             "No ions - the system is neutral only if the solute itself is neutral."))
    return checks


def cryst1_line(pdb_text: str) -> Optional[str]:
    return next((line for line in pdb_text.splitlines() if line.startswith("CRYST1")), None)


def restore_cryst1(pdb_text: str, line: Optional[str]) -> str:
    """Vrátí CRYST1 do připraveného PDB - builder ho bez solvatace nezapíše."""
    if not line or cryst1_line(pdb_text):
        return pdb_text
    return line.rstrip("\n") + "\n" + pdb_text
