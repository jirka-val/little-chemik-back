"""
Topologie a souřadnice nahrané v upload dialogu.

- Spolu s PDB: jen se uloží do workspace a zobrazí se z nich souhrn.
- Bez PDB (uživatel má jen simulaci z jiného FF): structure.pdb se zapíše
  z topologie a souřadnic (write_pdb_from_inputs -> AMBER_topology.savePdb)
  a dál jde běžná příprava. Box se bere ze souřadnic (rst7/inpcrd/ncrst mají
  aktuální box po NPT, se všemi třemi úhly); prmtop BOX_DIMENSIONS je jen
  záloha - je to box ze začátku a umí jen jeden úhel.

Ukládají se pod vlastními jmény, aby je nepřepsala topologie, kterou
generuje TopologyService (structure.prmtop / structure.crd).
"""

from __future__ import annotations

import json
import math
import struct
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from app.services.structure.solvated_system import box_shape, format_cryst1
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

    # Box po cpptraj (AMBER_topology.getBox): BOX_DIMENSIONS = [beta, a, b, c],
    # 60 = rhombický dodekaedr (alpha = gamma = 60, beta = 90).
    box = None
    try:
        cell = AMBER_topology.getBox(top)
    except ValueError:
        cell = []
    if cell:
        params = dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), cell))
        box = {**params, "shape": box_shape(params)}

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

class CoordinateError(ValueError):
    """Souřadnice nejdou přečíst - text je pro uživatele."""


_BOX_KEYS = ("a", "b", "c", "alpha", "beta", "gamma")
# NetCDF3 nc_type -> (struct formát, velikost)
_NC_TYPES = {1: ("b", 1), 2: ("c", 1), 3: ("h", 2), 4: ("i", 4), 5: ("f", 4), 6: ("d", 8)}


def _box_dict(values) -> Optional[Dict[str, float]]:
    box = {k: float(v) for k, v in zip(_BOX_KEYS, values)}
    if len(box) != 6 or not all(math.isfinite(v) for v in box.values()):
        return None
    if min(box["a"], box["b"], box["c"]) <= 0 or not all(0 < box[k] < 180 for k in ("alpha", "beta", "gamma")):
        return None
    return box


def _read_netcdf3(data: bytes) -> Tuple[Dict[str, int], Dict[str, np.ndarray]]:
    """Rozměry a proměnné NetCDF3 (classic / 64-bit offset) - restart z pmemd/cpptraj.
    U proměnných přes neomezený rozměr (trajektorie) jen první záznam."""
    if len(data) < 8 or data[:3] != b"CDF" or data[3] not in (1, 2):
        raise CoordinateError("Not a NetCDF3 coordinate file.")
    offset_size = 4 if data[3] == 1 else 8
    pos = 8  # magic + verze + numrecs

    def take(fmt: str):
        nonlocal pos
        size = struct.calcsize(fmt)
        if pos + size > len(data):
            raise CoordinateError("The NetCDF file is truncated.")
        values = struct.unpack(fmt, data[pos:pos + size])
        pos += size
        return values

    def name() -> str:
        nonlocal pos
        (length,) = take(">i")
        text = data[pos:pos + length].decode("ascii", errors="replace")
        pos += (length + 3) // 4 * 4
        return text

    def skip_attributes():
        nonlocal pos
        tag, count = take(">ii")
        for _ in range(count if tag == 0x0C else 0):
            name()
            nc_type, nelems = take(">ii")
            pos += (nelems * _NC_TYPES.get(nc_type, ("b", 1))[1] + 3) // 4 * 4

    dims: List[Tuple[str, int]] = []
    tag, count = take(">ii")
    for _ in range(count if tag == 0x0A else 0):
        dims.append((name(), take(">i")[0]))
    skip_attributes()  # globální atributy
    variables: Dict[str, np.ndarray] = {}
    tag, count = take(">ii")
    for _ in range(count if tag == 0x0B else 0):
        var_name = name()
        (ndims,) = take(">i")
        dim_ids = take(f">{ndims}i") if ndims else ()
        skip_attributes()
        nc_type, _vsize = take(">ii")
        (begin,) = take(">i" if offset_size == 4 else ">q")
        if nc_type not in _NC_TYPES or nc_type == 2:
            continue
        shape = [dims[i][1] for i in dim_ids if dims[i][1] != 0]  # 0 = neomezený rozměr
        fmt, size = _NC_TYPES[nc_type]
        n = int(np.prod(shape)) if shape else 1
        if begin + n * size > len(data):
            continue
        variables[var_name] = np.frombuffer(data, dtype=f">{fmt}", count=n, offset=begin).reshape(shape)
    return dict(dims), variables


def read_coordinates(path: Path) -> Dict[str, Any]:
    """Souřadnice (n x 3, A) a box z AMBER restartu/inpcrd (ASCII nebo NetCDF3).
    Vrací {"format", "title", "atoms", "time_ps", "velocities", "xyz", "box"}."""
    data = path.read_bytes()
    if data[:4] == b"\x89HDF":
        raise CoordinateError("NetCDF4 (HDF5) coordinates are not supported - convert them to .rst7 with cpptraj.")
    if data[:3] == b"CDF":
        dims, variables = _read_netcdf3(data)
        xyz = variables.get("coordinates")
        if xyz is None or xyz.ndim != 2 or xyz.shape[1] != 3:
            raise CoordinateError("The NetCDF file has no coordinates (n x 3).")
        box = None
        if "cell_lengths" in variables and "cell_angles" in variables:
            box = _box_dict(list(variables["cell_lengths"].ravel()) + list(variables["cell_angles"].ravel()))
        time = variables.get("time")
        return {
            "format": "NetCDF restart",
            "title": None,
            "atoms": int(xyz.shape[0]),
            "time_ps": float(time.ravel()[0]) if time is not None and time.size else None,
            "velocities": "velocities" in variables,
            "xyz": np.asarray(xyz, dtype=float),
            "box": box,
        }

    lines = data.decode("utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise CoordinateError("The coordinate file is empty or truncated.")
    fields = lines[1].split()
    try:
        natom = int(fields[0])
    except (IndexError, ValueError):
        raise CoordinateError("Not an AMBER ASCII coordinate file (atom count on line 2 is missing).")
    try:
        time = float(fields[1]) if len(fields) > 1 else None
    except ValueError:
        time = None

    # Formát 6F12.7: pevná pole po 12 znacích. Souřadnice, volitelně rychlosti
    # a volitelně box (a b c alpha beta gamma) na posledním řádku.
    values: List[float] = []
    try:
        for line in lines[2:]:
            text = line.rstrip()
            values.extend(float(text[i:i + 12]) for i in range(0, len(text), 12) if text[i:i + 12].strip())
    except ValueError:
        raise CoordinateError("The coordinate file contains a value that is not a number (expected 6F12.7).")
    n3 = 3 * natom
    if len(values) not in (n3, n3 + 6, 2 * n3, 2 * n3 + 6):
        raise CoordinateError(f"The coordinate file should hold {natom} atoms, but it has {len(values)} values "
                              f"(expected {n3} or {n3 + 6}, with velocities {2 * n3} or {2 * n3 + 6}).")
    box = _box_dict(values[-6:]) if len(values) in (n3 + 6, 2 * n3 + 6) else None
    return {
        "format": "ASCII restart" if len(values) >= 2 * n3 else "ASCII coordinates",
        "title": lines[0].strip() or None,
        "atoms": natom,
        "time_ps": time,
        "velocities": len(values) >= 2 * n3,
        "xyz": np.array(values[:n3], dtype=float).reshape(natom, 3),
        "box": box,
    }


def summarize_coordinates(path: Path) -> Dict[str, Any]:
    try:
        coords = read_coordinates(path)
    except CoordinateError as exc:
        return {"ok": False, "error": str(exc)}
    coords.pop("xyz")
    if coords["box"]:
        coords["box"] = {**coords["box"], "shape": box_shape(coords["box"])}
    return {"ok": True, **coords}


# ---------------------------------------------------------------------------
# PDB z topologie a souřadnic (uživatel nemá PDB)
# ---------------------------------------------------------------------------

class InputsError(ValueError):
    """Z nahrané topologie a souřadnic nejde zapsat PDB - text je pro uživatele."""


_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN", "GLU", "GLH", "GLY", "HIS", "HID", "HIE",
    "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# 4znakové vody (CHARMM, OPC ...) -> název, který zná savePdb i zbytek přípravy
_WATER_ALIASES = {"TIP3", "TIP4", "TIP5", "SPCE", "T3P", "T4P", "T4E", "OPC3"}


def _pdb_residue_name(label: str) -> str:
    """PDB má na název rezidua 3 znaky. AMBER knihovny mají terminály po 4
    (NTHR, CASN) - jako ambpdb/pdb4amber se zapíše jen reziduum."""
    if len(label) == 4 and label[0] in "NC" and label[1:] in _AMINO_ACIDS:
        return label[1:]
    if label.upper() in _WATER_ALIASES:
        return "WAT"
    return label


_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
_MAX_RESIDUES_PER_CHAIN = 9999


def _split_long_chains(pdb_text: str) -> str:
    """savePdb dá všechnu vodu (ionty) do jednoho řetězce a čísla reziduí za
    9999 začnou znovu od 1 - různé vody pak mají stejné (řetězec, číslo) a
    analýza je slije dohromady. Jako náš builder (forge_molecule_solvation):
    nad 9999 reziduí se řetězec rozdělí na další s čísly 1-9999."""
    sections: List[List[List[str]]] = [[]]   # sekce (do TER) -> rezidua -> řádky ATOM
    previous = None
    tail: List[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            key = (line[17:20], line[21], line[22:27])
            if key != previous:
                sections[-1].append([])
                previous = key
            sections[-1][-1].append(line)
        elif line.startswith("TER"):
            sections.append([])
            previous = None
        elif line.strip() and not line.startswith("END"):
            tail.append(line)
    sections = [s for s in sections if s]
    if all(len(s) <= _MAX_RESIDUES_PER_CHAIN for s in sections):
        return pdb_text

    used = {s[0][0][21] for s in sections}
    spare = iter(c for c in _CHAIN_IDS if c not in used)
    out: List[str] = list(tail)
    serial = 0

    def next_serial() -> int:
        nonlocal serial
        serial = serial % 99999 + 1
        return serial

    for section in sections:
        chain = section[0][0][21]
        for start in range(0, len(section), _MAX_RESIDUES_PER_CHAIN):
            chunk = section[start:start + _MAX_RESIDUES_PER_CHAIN]
            if start:
                chain = next(spare, None)
                if chain is None:
                    raise InputsError("The system has too many residues to be written to a PDB (out of chain IDs).")
            renumber = len(section) > _MAX_RESIDUES_PER_CHAIN
            for number, residue in enumerate(chunk, start=1):
                for line in residue:
                    resseq = f"{number:4d} " if renumber else line[22:27]
                    out.append(f"{line[:6]}{next_serial():5d}{line[11:21]}{chain}{resseq}{line[27:]}")
            last = chunk[-1][0]
            ter = f"TER   {next_serial():5d}      {last[17:20]} {chain}{(len(chunk) if renumber else int(last[22:26])):4d}"
            out.append(ter.ljust(80))
    out.append("END".ljust(80))
    return "\n".join(out) + "\n"


def write_pdb_from_inputs(topology_path: Path, coordinates_path: Path, pdb_path: Path) -> Dict[str, Any]:
    """structure.pdb z AMBER topologie a souřadnic (AMBER_topology.savePdb).

    Box: ze souřadnic, jinak z prmtopu (getBox), jinak žádný CRYST1. Vrací
    {"atoms", "box", "box_from": "coordinates" | "topology" | None, "renamed"}.
    """
    top = AMBER_topology.read_AMBER_topology(str(topology_path))
    if not top or "POINTERS" not in top:
        raise InputsError("The topology could not be read as an AMBER topology (prmtop). "
                          "GROMACS .top files are not supported.")
    if not top.get("ATOMIC_NUMBER"):
        raise InputsError("The topology has no ATOMIC_NUMBER section (written by an old AMBER version) - "
                          "elements cannot be assigned. Re-save it with ParmEd or tleap.")
    try:
        coords = read_coordinates(coordinates_path)
    except CoordinateError as exc:
        raise InputsError(str(exc))
    natom = top["POINTERS"][0]
    if coords["atoms"] != natom:
        raise InputsError(f"The topology has {natom} atoms but the coordinates have {coords['atoms']} - "
                          "they do not belong to the same system.")

    labels = [str(value).strip() for value in top["RESIDUE_LABEL"]]
    names = [_pdb_residue_name(label) for label in labels]
    renamed = dict(sorted(Counter(f"{a} -> {b}" for a, b in zip(labels, names) if a != b).items()))
    too_long = sorted({name for name in names if len(name) > 3})
    if too_long:
        raise InputsError(f"Residue names longer than 3 characters cannot be written to a PDB: {', '.join(too_long)}.")
    long_atoms = sorted({str(a).strip() for a in top["ATOM_NAME"] if len(str(a).strip()) > 4})
    if long_atoms:
        raise InputsError(f"Atom names longer than 4 characters cannot be written to a PDB: {', '.join(long_atoms)}.")

    box, box_from = coords["box"], "coordinates" if coords["box"] else None
    if box is None:
        try:
            cell = AMBER_topology.getBox(top)
        except ValueError:
            cell = []
        if cell:
            box, box_from = dict(zip(_BOX_KEYS, cell)), "topology"

    # CRYST1 se dopíše zvlášť - savePdb by ho vzal z prmtopu (getBox).
    pdb_top = {k: v for k, v in top.items() if k != "BOX_DIMENSIONS"}
    pdb_top["RESIDUE_LABEL"] = names
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "structure.pdb"
        try:
            AMBER_topology.savePdb(str(out), pdb_top, coords["xyz"].tolist())
        except ValueError as exc:
            raise InputsError(f"The PDB could not be written: {exc}.")
        body = _split_long_chains(out.read_text(encoding="ascii"))
    header = format_cryst1(box).ljust(80) + "\n" if box else ""
    pdb_path.write_text(header + body, encoding="utf-8")
    return {"atoms": natom, "box": box, "box_from": box_from, "renamed": renamed}


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
