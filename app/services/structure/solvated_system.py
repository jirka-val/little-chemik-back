"""
Rozpozná, že nahrané PDB je už hotový solvatovaný systém (box v CRYST1,
voda kolem solutu, ionty), a zkontroluje, jestli se dá použít tak, jak je.

Uživatel pak volí (PreparationRequest.existing_solvent):
- "keep"    nechat vodu, ionty i box, nic nepřidávat (jen dostavět, co chybí)
- "rebuild" vodu a monovalentní ionty odstranit a solvatovat znovu

Krystalová struktura z RCSB má CRYST1 taky (jednotková buňka) a pár desítek
krystalových vod - od solvatovaného boxu se liší hustotou vody v buňce.

PDB bez CRYST1 (export z VMD/PyMOL, převod z .gro ...): když je v něm voda
s hustotou kapaliny, je to taky solvatovaný systém. Box se pak rekonstruuje
z rozložení vody (box_reconstruction.py). Kontroly i volba "keep" pracují
se strukturou otočenou a posunutou do standardní orientace boxu
(standardize_reconstructed_box). Když tvar boxu poznat nejde, "keep" není
možné - zbývá "rebuild".
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from app.services.pdb_service import ion_group_for_resname
from app.services.structure.box_reconstruction import MIN_WATERS as LIQUID_MIN_WATERS
from app.services.structure.box_reconstruction import BoxReconstruction, reconstruct_box, standard_vectors

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP", "DOD"}

# Čistá voda při 298 K: 0.997 g/cm3 = 0.0334 molekul na A^3.
WATER_NUMBER_DENSITY = 0.0334
# Parciální specifický objem proteinu ~0.73 cm3/g = 1.21 A^3 na Da.
SOLUTE_VOLUME_PER_DALTON = 1.21

# Solvatovaný box má v buňce řádově 0.02-0.03 vod/A^3, krystal ~0.001.
SOLVATED_MIN_WATERS = 50
SOLVATED_MIN_WATER_DENSITY = 0.01

# Bez CRYST1: voda je kapalina, když má medián nejbližší vzdálenosti O-O pod
# touto mezí (kapalná voda ~2.8 A, krystalové vody jsou od sebe řádově dál).
LIQUID_MAX_NN = 3.3

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
    # rhombický dodekaedr: dvakrát 60 a jednou 90 stupňů (GROMACS 60/60/90,
    # cpptraj/AMBER 60/90/60)
    if sorted(round(x) for x in angles) == [60, 60, 90] and all(
            min(abs(x - 60), abs(x - 90)) < 0.1 for x in angles):
        return "rhombic dodecahedron"
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


# O-O vody blíž = dvě molekuly přes sebe (v kapalině nejblíž ~2,5 A).
WATER_OVERLAP_DISTANCE = 2.2


def periodic_water_overlaps(oxygens: np.ndarray, vectors: np.ndarray,
                            cutoff: float = WATER_OVERLAP_DISTANCE) -> int:
    """Kolik molekul vody se překrývá s periodickým obrazem jiné molekuly.

    Ve správném boxu voda přes hranu plynule navazuje - žádné překryvy. Box,
    který k souřadnicím nepatří (jiná mřížka, box z jiného snímku simulace,
    prmtop s jedním úhlem místo celé buňky), dá na hranách stovky překryvů,
    přitom ostatní kontroly (hustota, překryvy uvnitř) projdou.
    Kyslíky se zabalí do rovnoběžnostěnu mřížky (platí pro každý tvar buňky)
    a k těm u stěn se přidají obrazy.
    """
    if len(oxygens) < 2:
        return 0
    vectors = np.asarray(vectors, dtype=float)
    frac = np.asarray(oxygens, dtype=float) @ np.linalg.inv(vectors)
    frac -= np.floor(frac)
    wrapped = frac @ vectors
    volume = abs(np.linalg.det(vectors))
    heights = [volume / np.linalg.norm(np.cross(vectors[(i + 1) % 3], vectors[(i + 2) % 3])) for i in range(3)]
    margin = np.array([cutoff / h for h in heights])
    near_low, near_high = frac < margin, frac > 1.0 - margin
    images, owners = [], []
    for shift in itertools.product((-1, 0, 1), repeat=3):
        if shift == (0, 0, 0):
            continue
        # obraz posunutý o +1 je blízko spodní stěny, jen když atom leží u horní
        mask = np.ones(len(frac), dtype=bool)
        for axis, step in enumerate(shift):
            if step == 1:
                mask &= near_high[:, axis]
            elif step == -1:
                mask &= near_low[:, axis]
        if mask.any():
            images.append(wrapped[mask] - np.asarray(shift, dtype=float) @ vectors)
            owners.append(np.flatnonzero(mask))
    if not images:
        return 0
    n = len(wrapped)
    owner = np.concatenate([np.arange(n)] + owners)
    pairs = _close_pairs(np.vstack([wrapped] + images), cutoff)
    across = (pairs[:, 0] >= n) ^ (pairs[:, 1] >= n)          # původní atom vs. obraz
    pairs = pairs[across]
    a, b = owner[pairs[:, 0]], owner[pairs[:, 1]]
    keep = a != b
    return int(len(np.unique(np.concatenate([a[keep], b[keep]]))))


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


def _water_oxygens(atoms: List[_Atom]) -> np.ndarray:
    return np.array([a.xyz for a in atoms if a.resname in WATER_RESNAMES and a.element == "O"], dtype=float)


def _is_liquid_water(oxygens: np.ndarray) -> bool:
    """Hodně vody s hustotou kapaliny (bez boxu nejde spočítat hustotu v buňce)."""
    if len(oxygens) < LIQUID_MIN_WATERS:
        return False
    sample = oxygens[np.linspace(0, len(oxygens) - 1, min(400, len(oxygens))).astype(int)]
    nearest = []
    for start in range(0, len(sample), 100):
        chunk = sample[start:start + 100]
        d2 = ((chunk[:, None, :] - oxygens[None, :, :]) ** 2).sum(axis=2)
        d2[d2 < 1e-6] = np.inf
        nearest.append(np.sqrt(d2.min(axis=1)))
    return float(np.median(np.concatenate(nearest))) < LIQUID_MAX_NN


def _cryst1_box(params: Dict[str, float]) -> Dict[str, float]:
    """Parametry v přesnosti CRYST1, aby kontroly i zapsaný řádek seděly."""
    return {k: round(float(params[k]), 3 if k in ("a", "b", "c") else 2)
            for k in ("a", "b", "c", "alpha", "beta", "gamma")}


def format_cryst1(box: Dict[str, float]) -> str:
    return (f"CRYST1{box['a']:9.3f}{box['b']:9.3f}{box['c']:9.3f}"
            f"{box['alpha']:7.2f}{box['beta']:7.2f}{box['gamma']:7.2f} P 1           1")


def _override_params(box_override: Optional[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    return _cryst1_box(box_override) if box_override else None


def _adjusted(reconstruction: BoxReconstruction, params: Dict[str, float]) -> BoxReconstruction:
    """Rekonstrukce s rozměry/tvarem od uživatele: natočení a střed zůstanou,
    struktura se vycentruje do nového boxu."""
    return replace(
        reconstruction, params=dict(params),
        std_center=standard_vectors(params).sum(axis=0) / 2,
        volume=float(abs(np.linalg.det(standard_vectors(params)))),
    )


def _reconstruct(parsed: _Parsed) -> Optional[BoxReconstruction]:
    """Rekonstrukce boxu pro PDB bez CRYST1, None když v něm není kapalná voda."""
    if parsed.box:
        return None
    oxygens = _water_oxygens(parsed.atoms)
    if not _is_liquid_water(oxygens):
        return None
    return reconstruct_box(oxygens)


def analyze_solvated_system(
    pdb_text: str, with_checks: bool = True, reconstruct: bool = True,
    box_override: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    {"solvated": bool, ...}. Když je systém solvatovaný, obsahuje i box,
    počty a seznam kontrol (status ok/warning/error) pro volbu "keep".

    box_source: "cryst1" | "reconstructed" (bez CRYST1, box z rozložení
    vody, kontroly jsou ve standardní orientaci boxu) | "missing" (kapalná
    voda bez CRYST1, tvar boxu poznat nejde) | None (žádný box).
    reconstruct=False jen pozná solvataci, box nerekonstruuje (rychlé).

    box_override (Expert, jen PDB bez CRYST1): {"mode": "reconstructed" |
    "manual", a, b, c, alpha, beta, gamma}. "reconstructed" = rekonstruovaný
    box s upraveným tvarem/rozměry (natočení z rekonstrukce), "manual" = box
    zadaný ručně pro souřadnice tak, jak jsou (jako dopsaný CRYST1) ->
    box_source "manual".
    """
    parsed = _parse(pdb_text)
    box_source = "cryst1" if parsed.box else None
    reconstruction: Optional[BoxReconstruction] = None
    liquid_without_box = False
    if not parsed.box:
        if reconstruct:
            reconstruction = _reconstruct(parsed)
            liquid_without_box = reconstruction is not None
        else:
            liquid_without_box = _is_liquid_water(_water_oxygens(parsed.atoms))
    override = _override_params(box_override) if box_source is None else None
    if override is not None and box_override.get("mode") == "manual":
        box_source = "manual"
        parsed.box = override
    elif reconstruction is not None and reconstruction.ok:
        if override is not None:
            reconstruction = _adjusted(reconstruction, override)
        box_source = "reconstructed"
        parsed.box = _cryst1_box(reconstruction.params)
        moved = reconstruction.transform(np.array([a.xyz for a in parsed.atoms], dtype=float))
        parsed.atoms = [replace(a, xyz=tuple(float(v) for v in xyz)) for a, xyz in zip(parsed.atoms, moved)]
    elif liquid_without_box:
        box_source = "missing"

    residues = _residues(parsed.atoms)
    waters = {k: v for k, v in residues.items() if v[0].resname in WATER_RESNAMES}
    ions = {k: v for k, v in residues.items() if k not in waters and len(v) == 1 and v[0].element != "H"}
    solute_atoms = [a for k, v in residues.items() if k not in waters and k not in ions for a in v]

    box = parsed.box
    volume = float(abs(np.linalg.det(box_vectors(box)))) if box else None
    water_density = len(waters) / volume if volume else 0.0
    solvated = box_source == "missing" or (
        bool(box) and len(waters) >= SOLVATED_MIN_WATERS and water_density >= SOLVATED_MIN_WATER_DENSITY)

    ion_counts: Dict[str, int] = {}
    for atoms in ions.values():
        ion_counts[atoms[0].resname] = ion_counts.get(atoms[0].resname, 0) + 1

    result: Dict[str, Any] = {
        "solvated": solvated,
        "box": {**box, "shape": box_shape(box), "volume": round(volume, 1)} if box else None,
        "box_source": box_source if solvated else None,
        "waters": len(waters),
        "ions": dict(sorted(ion_counts.items())),
        "solute_atoms": len(solute_atoms),
        "solute_residues": len({a.res_key for a in solute_atoms}),
    }
    if box_source == "reconstructed":
        result["box_reconstruction"] = reconstruction.summary()
    if override is not None and box_source in ("reconstructed", "manual"):
        result["box_override"] = box_override.get("mode")
    if solvated and with_checks:
        if box_source == "missing":
            reason = reconstruction.reason if reconstruction is not None else "not attempted"
            result["checks"] = [_check(
                "box_source", "Box", "error",
                f"No CRYST1 record and the box could not be reconstructed from the water ({reason}). "
                "The existing solvent cannot be kept - remove it and solvate again.")]
        else:
            checks = _checks(box, volume, waters, ion_counts, solute_atoms, parsed.atoms)
            size = f"{box_shape(box)}, {box['a']:.2f} x {box['b']:.2f} x {box['c']:.2f} A"
            if box_source == "reconstructed":
                how = "reconstructed from the water and adjusted by hand" if override else "reconstructed from the water"
                checks.insert(0, _check(
                    "box_source", "Box", "warning",
                    f"No CRYST1 record - the box was {how}: {size}. Keeping it moves the structure "
                    "into the standard orientation of this box; NPT equilibration will settle its size."))
            elif box_source == "manual":
                checks.insert(0, _check(
                    "box_source", "Box", "warning",
                    f"No CRYST1 record - box entered by hand: {size}, angles {box['alpha']:g}/{box['beta']:g}/"
                    f"{box['gamma']:g}. The coordinates are used as they are."))
            result["checks"] = checks
    return result


def standardize_reconstructed_box(
    pdb_text: str, box_override: Optional[Dict[str, Any]] = None,
) -> Optional[Tuple[str, str]]:
    """
    Pro "keep" u PDB bez CRYST1: (PDB otočené a posunuté do standardní
    orientace rekonstruovaného boxu s CRYST1 na začátku, řádek CRYST1).
    None, když box rekonstruovat nejde (nebo CRYST1 v souboru je).
    ANISOU se po rotaci zahodí. box_override viz analyze_solvated_system
    ("manual" = jen dopsaný CRYST1, souřadnice beze změny).
    """
    parsed = _parse(pdb_text)
    if parsed.box:
        return None
    override = _override_params(box_override)
    if override is not None and box_override.get("mode") == "manual":
        cryst1 = format_cryst1(override)
        lines = [line for line in pdb_text.splitlines() if not line.startswith("CRYST1")]
        return cryst1 + "\n" + "\n".join(lines) + "\n", cryst1
    reconstruction = _reconstruct(parsed)
    if reconstruction is None or not reconstruction.ok:
        return None
    if override is not None:
        reconstruction = _adjusted(reconstruction, override)
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ANISOU"):
            continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]])
            except ValueError:
                out.append(line)
                continue
            x, y, z = reconstruction.transform(xyz)[0]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        out.append(line)
    cryst1 = format_cryst1(_cryst1_box(reconstruction.params))
    return cryst1 + "\n" + "\n".join(out) + "\n", cryst1


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

    # 6. voda navazuje přes hranu boxu (box patří k souřadnicím)
    oxygens = np.array([a.xyz for atoms in waters.values() for a in atoms if a.element == "O"], dtype=float)
    if len(oxygens):
        overlapping = periodic_water_overlaps(oxygens, box_vectors(box))
        if overlapping > max(3, 0.002 * len(oxygens)):
            checks.append(_check("periodic_continuity", "Water across the box boundary", "error",
                                 f"{overlapping} water molecules overlap with a periodic image of another one - "
                                 "the box (CRYST1) does not belong to these coordinates (different lattice, or "
                                 "the box comes from another frame of the simulation)."))
        elif overlapping:
            checks.append(_check("periodic_continuity", "Water across the box boundary", "warning",
                                 f"{overlapping} water molecules are closer than {WATER_OVERLAP_DISTANCE} A to a "
                                 "periodic image of another one - minimization should resolve it."))
        else:
            checks.append(_check("periodic_continuity", "Water across the box boundary", "ok",
                                 "continuous - no overlaps with periodic images"))

    # 7. ionty - známé FF a jejich náboj
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


def strip_solvent(pdb_text: str) -> Tuple[str, Dict[str, Any]]:
    """
    Odstraní vodu, monovalentní ionty a CRYST1 nahraného solvatovaného
    systému ("solvatovat znovu" hned po načtení). Multivalentní ionty
    (Mg, Zn ... často strukturální) zůstanou - o nich rozhodne příprava.
    Spolu s atomy zmizí i jejich ANISOU, CONECT a osiřelé TER.
    Vrací (PDB, {"waters": n, "ions": {jméno: počet}}).
    """
    lines = pdb_text.splitlines()
    # Rezidua jako souvislé běhy řádků se stejným (řetězec, číslo, jméno) -
    # čísla vod se u velkých boxů opakují (přetečení přes 9999).
    runs: List[List[int]] = []
    last_key = None
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21], line[22:27], line[17:21].strip())
        if key != last_key:
            runs.append([])
            last_key = key
        runs[-1].append(i)

    removed_lines = set()
    removed_serials = set()
    waters = 0
    ions: Dict[str, int] = {}
    for run in runs:
        first = lines[run[0]]
        resname = first[17:21].strip()
        if resname in WATER_RESNAMES:
            waters += 1
        elif len(run) == 1 and _element(first) != "H" and _ion_charge(resname) in (1, -1):
            ions[resname] = ions.get(resname, 0) + 1
        else:
            continue
        removed_lines.update(run)
        removed_serials.update(lines[i][6:11].strip() for i in run)

    out = []
    previous_atom_removed = False
    for i, line in enumerate(lines):
        if line.startswith("CRYST1") or i in removed_lines:
            previous_atom_removed = previous_atom_removed or i in removed_lines
            continue
        if line.startswith("ANISOU") and line[6:11].strip() in removed_serials:
            continue
        if line.startswith("CONECT"):
            serials = {line[k:k + 5].strip() for k in range(6, len(line), 5)} - {""}
            if serials & removed_serials:
                continue
        if line.startswith("TER") and previous_atom_removed:
            continue
        if line.startswith(("ATOM", "HETATM", "TER")):
            previous_atom_removed = False
        out.append(line)
    return "\n".join(out) + "\n", {"waters": waters, "ions": dict(sorted(ions.items()))}


def cryst1_line(pdb_text: str) -> Optional[str]:
    return next((line for line in pdb_text.splitlines() if line.startswith("CRYST1")), None)


def restore_cryst1(pdb_text: str, line: Optional[str]) -> str:
    """Vrátí CRYST1 do připraveného PDB - builder ho bez solvatace nezapíše."""
    if not line or cryst1_line(pdb_text):
        return pdb_text
    return line.rstrip("\n") + "\n" + pdb_text
