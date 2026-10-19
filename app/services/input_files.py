"""
Topologie a souřadnice nahrané spolu s PDB (upload dialog). Zatím se jen
uloží do workspace a zobrazí se z nich souhrn - do přípravy ani exportu
nevstupují (navázání na funkce typu get_box přijde později).

Ukládají se pod vlastními jmény, aby je nepřepsala topologie, kterou
generuje TopologyService (structure.prmtop / structure.crd).
"""

from __future__ import annotations

import json
import struct
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.utils.adams4sims_processing_library.utils import AMBER_topology

TOPOLOGY_EXTENSIONS = (".prmtop", ".parm7", ".top")
COORDINATE_EXTENSIONS = (".crd", ".inpcrd", ".rst7", ".rst", ".ncrst")

UPLOADED_TOPOLOGY = "uploaded.prmtop"
UPLOADED_COORDINATES = "uploaded_coordinates"  # + původní přípona
MANIFEST = "uploaded_inputs.json"

AMBER_CHARGE_UNIT = 18.2223  # prmtop CHARGE je v e * 18.2223
_WATER = {"WAT", "HOH", "SOL"}
_GROUP_LABELS = {"P": "protein", "R": "RNA", "D": "DNA"}


def load_manifest(workspace_dir: Path) -> Dict[str, Any]:
    path = workspace_dir / MANIFEST
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_manifest(workspace_dir: Path, manifest: Dict[str, Any]) -> None:
    (workspace_dir / MANIFEST).write_text(json.dumps(manifest, indent=1), encoding="utf-8")


def coordinates_filename(original: str) -> str:
    return UPLOADED_COORDINATES + Path(original).suffix.lower()


# ---------------------------------------------------------------------------
# Topologie
# ---------------------------------------------------------------------------

def summarize_topology(path: Path) -> Dict[str, Any]:
    top = AMBER_topology.read_AMBER_topology(str(path))
    if not top or "POINTERS" not in top:
        return {"ok": False, "error": "The file could not be read as an AMBER topology (prmtop)."}
    pointers = top["POINTERS"]
    natom, nres = pointers[0], pointers[11]
    labels = top.get("RESIDUE_LABEL", [])
    pointer = top.get("RESIDUE_POINTER", [])
    sizes = [(pointer[i + 1] if i + 1 < nres else natom + 1) - pointer[i] for i in range(nres)]

    # Typ molekuly po reziduích stejným slovníkem, jaký používá generování topologie.
    try:
        AMBER_topology.read_mol_type(top)
        per_res = top["MOLECULE_TYPE_PER_RES"]
    except Exception:
        per_res = ["unk"] * nres

    water = sum(1 for name in labels if name.strip() in _WATER)
    ions: Counter = Counter()
    groups: Counter = Counter()
    other: Counter = Counter()
    for name, size, mol_type in zip(labels, sizes, per_res):
        name = name.strip()
        if name in _WATER:
            continue
        if size == 1:
            ions[name] += 1
        elif mol_type in _GROUP_LABELS:
            groups[_GROUP_LABELS[mol_type]] += 1
        else:
            other[name] += 1

    atom_names = top.get("ATOM_NAME", [])
    atomic_numbers = top.get("ATOMIC_NUMBER", [])
    masses = top.get("MASS", [])
    extra_points = sum(1 for z in atomic_numbers if z == 0) if atomic_numbers else \
        sum(1 for n in atom_names if n.startswith("EP"))
    hydrogen_masses = [m for z, m in zip(atomic_numbers, masses) if z == 1]

    box = None
    ifbox = pointers[27] if len(pointers) > 27 else 0
    if ifbox:
        dims = top.get("BOX_DIMENSIONS", [])
        if len(dims) == 4:
            beta, a, b, c = dims
            box = {"a": a, "b": b, "c": c, "beta": beta,
                   "shape": "truncated octahedron" if abs(beta - 109.4712190) < 0.01 else
                            "rectangular" if abs(beta - 90.0) < 0.01 else "other"}

    molecules = len(top["ATOMS_PER_MOLECULE"]) if top.get("ATOMS_PER_MOLECULE") else None
    # TITLE je 20a4 - čtečka ho vrací po čtyřech znacích
    title = "".join(top.get("TITLE", [])).strip() or None
    return {
        "ok": True,
        "title": title,
        "atoms": natom,
        "residues": nres,
        "molecules": molecules,
        "composition": {
            **{k: v for k, v in sorted(groups.items())},
            "water": water,
            "ions": dict(sorted(ions.items())),
            "other": dict(sorted(other.items())),
        },
        "water_model_sites": _water_sites(labels, sizes),
        "extra_points": extra_points,
        "total_charge": round(sum(top.get("CHARGE", [])) / AMBER_CHARGE_UNIT, 3) + 0.0,  # bez -0.0
        "hmr": bool(hydrogen_masses) and max(hydrogen_masses) > 2.0,
        "box": box,
    }


def _water_sites(labels: List[str], sizes: List[int]) -> Optional[int]:
    for name, size in zip(labels, sizes):
        if name.strip() in _WATER:
            return size
    return None


# ---------------------------------------------------------------------------
# Souřadnice
# ---------------------------------------------------------------------------

def _netcdf_atom_count(data: bytes) -> Optional[int]:
    """Rozměr "atom" z hlavičky NetCDF3 (classic / 64-bit offset) - restart z pmemd."""
    if len(data) < 16 or data[:3] != b"CDF":
        return None
    pos = 8  # magic + verze + numrecs
    tag, count = struct.unpack(">ii", data[pos:pos + 8])
    pos += 8
    if tag != 0x0A:  # NC_DIMENSION
        return None
    for _ in range(count):
        (length,) = struct.unpack(">i", data[pos:pos + 4])
        name = data[pos + 4:pos + 4 + length].decode("ascii", errors="replace")
        pos += 4 + (length + 3) // 4 * 4
        (size,) = struct.unpack(">i", data[pos:pos + 4])
        pos += 4
        if name == "atom":
            return size
    return None


def summarize_coordinates(path: Path) -> Dict[str, Any]:
    data = path.read_bytes()
    if data[:3] == b"CDF":
        return {"ok": True, "format": "NetCDF restart", "atoms": _netcdf_atom_count(data), "box": None}
    if data[:3] == b"\x89HD":
        return {"ok": True, "format": "NetCDF4 (not inspected)", "atoms": None, "box": None}

    lines = data.decode("utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        return {"ok": False, "error": "The coordinate file is empty or truncated."}
    fields = lines[1].split()
    try:
        natom = int(fields[0])
    except (IndexError, ValueError):
        return {"ok": False, "error": "Not an AMBER ASCII coordinate file (atom count on line 2 is missing)."}
    time = float(fields[1]) if len(fields) > 1 else None

    # Box je poslední řádek se 6 hodnotami, pokud soubor má víc než jen
    # souřadnice (a případně rychlosti).
    values_per_line = 6
    coord_lines = (natom * 3 + values_per_line - 1) // values_per_line
    body = [line for line in lines[2:] if line.strip()]
    has_velocities = len(body) >= 2 * coord_lines
    box = None
    if len(body) in (coord_lines + 1, 2 * coord_lines + 1):
        last = body[-1]
        try:
            vals = [float(last[i:i + 12]) for i in range(0, 72, 12) if last[i:i + 12].strip()]
            if len(vals) == 6:
                box = {"a": vals[0], "b": vals[1], "c": vals[2], "alpha": vals[3], "beta": vals[4], "gamma": vals[5]}
        except ValueError:
            pass
    return {
        "ok": True,
        "format": "ASCII restart" if has_velocities else "ASCII coordinates",
        "title": lines[0].strip() or None,
        "atoms": natom,
        "time_ps": time,
        "velocities": has_velocities,
        "box": box,
    }


# ---------------------------------------------------------------------------
# Celkový souhrn
# ---------------------------------------------------------------------------

def _pdb_atoms(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    with open(path, encoding="utf-8", errors="replace") as fh:
        return sum(1 for line in fh if line.startswith(("ATOM", "HETATM")))


def summarize_inputs(workspace_dir: Path) -> Dict[str, Any]:
    manifest = load_manifest(workspace_dir)
    summary: Dict[str, Any] = {
        "pdb": {"atoms": _pdb_atoms(workspace_dir / "structure.pdb")},
        "topology": None,
        "coordinates": None,
        "warnings": [],
    }
    if "topology" in manifest:
        summary["topology"] = {"filename": manifest["topology"],
                               **summarize_topology(workspace_dir / UPLOADED_TOPOLOGY)}
    if "coordinates" in manifest:
        stored = workspace_dir / coordinates_filename(manifest["coordinates"])
        summary["coordinates"] = {"filename": manifest["coordinates"], **summarize_coordinates(stored)}

    top, crd, pdb_atoms = summary["topology"], summary["coordinates"], summary["pdb"]["atoms"]
    top_atoms = top.get("atoms") if top and top.get("ok") else None
    crd_atoms = crd.get("atoms") if crd and crd.get("ok") else None
    if top_atoms and crd_atoms and top_atoms != crd_atoms:
        summary["warnings"].append(f"Topology has {top_atoms} atoms but the coordinates have {crd_atoms}.")
    if top_atoms and pdb_atoms and top_atoms != pdb_atoms:
        summary["warnings"].append(f"Topology has {top_atoms} atoms but the PDB has {pdb_atoms}.")
    top_box = top.get("box") if top and top.get("ok") else None
    crd_box = crd.get("box") if crd and crd.get("ok") else None
    if top_box and crd_box and abs(top_box["beta"] - crd_box["beta"]) > 0.01:
        summary["warnings"].append(
            f"Box angle differs: topology {top_box['beta']:g} deg, coordinates {crd_box['beta']:g} deg."
        )
    if top and top.get("ok") and abs(top["total_charge"]) > 0.01:
        summary["warnings"].append(f"The system is not neutral (total charge {top['total_charge']:+g} e).")
    return summary
