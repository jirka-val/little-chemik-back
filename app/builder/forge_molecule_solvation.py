"""Periodic explicit-water solvation for built FORGE Molecule objects.

This post-builder layer supports orthorhombic, cubic and two coordinate
representations of a BCC/truncated-octahedral periodic cell.  It deliberately
adds only W3/WAT waters; downstream force-field processing may convert these
to four- or five-site water models.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from itertools import product
import json
import math
from pathlib import Path
from time import perf_counter
from typing import Any, Dict, Iterable, Iterator, Mapping, Optional, Sequence, Tuple

import numpy as np

from forge_molecule_builder import MMParameterProvider
from forge_molecule_parser import (
    Atom,
    Chain,
    Coord,
    Molecule,
    PDBAtomRecord,
    PeriodicBox,
    Residue,
    infer_element,
)


SUPPORTED_BOX_SHAPES = (
    "orthorhombic",
    "cubic",
    "truncated_octahedron",
    "truncated_octahedron_primitive",
)
WATER_GROUPS = frozenset(("W3", "W4", "W5"))


@dataclass(frozen=True)
class WaterGeometry:
    oxygen: np.ndarray
    hydrogen_offsets: Tuple[np.ndarray, np.ndarray]


@dataclass(frozen=True)
class WaterSolvationTemplate:
    mol_type: str
    residue_name: str
    oxygen_atom: str
    hydrogen_atoms: Tuple[str, str]
    cell_vectors: np.ndarray
    waters: Tuple[WaterGeometry, ...]


@dataclass
class SolvationSettings:
    box_shape: str = "truncated_octahedron"
    padding_angstrom: float = 10.0
    rotate_principal_axes: bool = True
    keep_crystal_waters: bool = False
    solvent_chain_id: Optional[str] = None

    def validate(self) -> None:
        if self.box_shape not in SUPPORTED_BOX_SHAPES:
            raise ValueError(
                f"Unsupported box shape {self.box_shape!r}; expected one of "
                f"{', '.join(SUPPORTED_BOX_SHAPES)}"
            )
        if not math.isfinite(self.padding_angstrom) or self.padding_angstrom < 0.0:
            raise ValueError("padding_angstrom must be a finite non-negative value")
        if self.solvent_chain_id is not None and len(self.solvent_chain_id) != 1:
            raise ValueError("solvent_chain_id must be one PDB-compatible character")


@dataclass
class SolvationReport:
    box_shape: str
    padding_angstrom: float
    box_vectors: Tuple[Coord, Coord, Coord]
    input_crystal_waters: int
    retained_crystal_waters: int
    template_candidates_in_box: int
    waters_removed_by_solute: int
    generated_waters: int
    total_waters: int
    solvent_chain_id: str
    solvent_chain_ids: Tuple[str, ...] = ()
    environment_atoms: int = 0
    spatial_pair_checks: int = 0
    timings_seconds: Dict[str, float] = field(default_factory=dict)


@dataclass
class SolvationVdwParameters:
    """mol_type-indexed force-field providers used only for LJ sigma lookup."""

    providers: Dict[str, MMParameterProvider]

    @classmethod
    def from_force_field_directories(
        cls, directories: Iterable[str | Path]
    ) -> "SolvationVdwParameters":
        providers: Dict[str, MMParameterProvider] = {}
        for raw_directory in directories:
            directory = Path(raw_directory)
            mol_type = directory.name.rsplit("_", 1)[-1]
            if mol_type in providers:
                raise ValueError(f"Duplicate force field for mol_type {mol_type!r}")
            providers[mol_type] = _load_mm_provider_from_directory(directory)
        return cls(providers)

    @classmethod
    def from_force_field_root(
        cls, root: str | Path
    ) -> "SolvationVdwParameters":
        directories = sorted(
            path for path in Path(root).iterdir() if path.is_dir() and "_" in path.name
        )
        return cls.from_force_field_directories(directories)

    def atom_params(self, mol_type: str, residue_name: str, atom_name: str) -> Any:
        try:
            residue = self.providers[mol_type].residues[residue_name]
        except KeyError as exc:
            raise KeyError(
                f"LJ sigma missing for {mol_type}:{residue_name}:{atom_name}"
            ) from exc
        if atom_name in residue:
            return residue[atom_name]
        folded = [value for name, value in residue.items() if name.casefold() == atom_name.casefold()]
        if len(folded) == 1:
            return folded[0]
        # Monatomic ions are unambiguous even when PDB/UMFFF capitalization
        # differs (for example CL versus Cl-).
        if len(residue) == 1:
            return next(iter(residue.values()))
        raise KeyError(f"LJ sigma missing for {mol_type}:{residue_name}:{atom_name}")

    def sigma(self, mol_type: str, residue_name: str, atom_name: str) -> float:
        return float(self.atom_params(mol_type, residue_name, atom_name).sigma)

    def merged_builder_provider(
        self, mol_types: Iterable[str] = ("R", "D", "P")
    ) -> MMParameterProvider:
        """Combine disjoint polymer parameter sets for mixed-system building."""
        selected = [self.providers[mol_type] for mol_type in mol_types if mol_type in self.providers]
        if not selected:
            raise ValueError("No requested polymer force-field providers are available")
        fudge_lj = selected[0].fudge_lj
        fudge_qq = selected[0].fudge_qq
        residues: Dict[str, Any] = {}
        bonds: Dict[str, Any] = {}
        dihedrals: Dict[Any, Any] = {}
        for provider in selected:
            if not math.isclose(provider.fudge_lj, fudge_lj) or not math.isclose(
                provider.fudge_qq, fudge_qq
            ):
                raise ValueError("Cannot merge force fields with different 1-4 scaling")
            for name, value in provider.residues.items():
                if name in residues and residues[name] != value:
                    raise ValueError(f"Conflicting residue parameters for {name!r}")
                residues[name] = value
            for name, value in provider.residue_bonds.items():
                if name in bonds and bonds[name] != value:
                    raise ValueError(f"Conflicting residue bonds for {name!r}")
                bonds[name] = value
            for key, value in provider.dihedraltypes.items():
                if key in dihedrals and dihedrals[key] != value:
                    raise ValueError(f"Conflicting torsion parameters for {key!r}")
                dihedrals[key] = value
        return MMParameterProvider(
            residues=residues,
            residue_bonds=bonds,
            dihedraltypes=dihedrals,
            fudge_lj=fudge_lj,
            fudge_qq=fudge_qq,
        )


@dataclass
class _RetainedWater:
    original_resname: str
    atoms: Dict[str, np.ndarray]
    original_atom_names: frozenset[str] = frozenset()


class _SpatialHash:
    """Small fixed-coordinate cell list; instances are local to solvation."""

    def __init__(self, coordinates: Sequence[np.ndarray], cell_size: float):
        if cell_size <= 0.0:
            raise ValueError("Spatial-hash cell size must be positive")
        self.coordinates = coordinates
        self.cell_size = float(cell_size)
        self.cells: Dict[Tuple[int, int, int], list[int]] = {}
        for index, coord in enumerate(coordinates):
            self.cells.setdefault(self._key(coord), []).append(index)

    def _key(self, coord: np.ndarray) -> Tuple[int, int, int]:
        return tuple(int(math.floor(float(value) / self.cell_size)) for value in coord)  # type: ignore[return-value]

    def candidates(self, coord: np.ndarray, radius: float) -> Iterator[int]:
        center = self._key(coord)
        reach = max(1, int(math.ceil(radius / self.cell_size)))
        for dx, dy, dz in product(range(-reach, reach + 1), repeat=3):
            yield from self.cells.get(
                (center[0] + dx, center[1] + dy, center[2] + dz), ()
            )


def _one_matching_file(
    directory: Path,
    labels: Sequence[str],
    exclude: Sequence[str] = (),
) -> Path:
    matches = [
        path
        for path in directory.iterdir()
        if path.is_file()
        and any(label in path.name.lower() for label in labels)
        and not any(label in path.name.lower() for label in exclude)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {'/'.join(labels)} file in {directory}, found "
            f"{[path.name for path in matches]}"
        )
    return matches[0]


def _load_mm_provider_from_directory(directory: Path) -> MMParameterProvider:
    residue_lib = _one_matching_file(directory, ("residue-lib", "residue_lib"))
    nonbonded = _one_matching_file(directory, ("nonbonded",))
    bonded = _one_matching_file(directory, ("bonded",), exclude=("nonbonded",))
    force_field = _one_matching_file(directory, ("forcefield", "force_field"))
    return MMParameterProvider.from_files(
        residue_lib, nonbonded, bonded, force_field
    )


def load_solvation_template(
    metadata_path: str | Path, mol_type: str = "W3"
) -> WaterSolvationTemplate:
    metadata_path = Path(metadata_path)
    data = json.loads(metadata_path.read_text(encoding="utf-8"))
    try:
        entry = data[mol_type]
    except KeyError as exc:
        raise KeyError(f"No solvation template for mol_type {mol_type!r}") from exc

    cell = np.asarray(entry["cell_vectors"], dtype=float)
    if cell.shape != (3, 3) or abs(float(np.linalg.det(cell))) < 1.0e-9:
        raise ValueError("Solvation-template cell_vectors must be a nonsingular 3x3 matrix")
    oxygen_name = str(entry["oxygen_atom"])
    hydrogen_names = tuple(str(name) for name in entry["hydrogen_atoms"])
    if len(hydrogen_names) != 2:
        raise ValueError("W3 solvation template must define exactly two hydrogens")

    try:
        raw_waters = entry["waters"]
    except KeyError as exc:
        raise ValueError("Solvation template must define inline waters") from exc

    waters: list[WaterGeometry] = []
    for index, raw_water in enumerate(raw_waters):
        try:
            oxygen = np.asarray(raw_water["oxygen"], dtype=float)
            offsets = np.asarray(raw_water["hydrogen_offsets"], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Malformed inline solvation-template water at index {index}"
            ) from exc
        if oxygen.shape != (3,) or offsets.shape != (2, 3):
            raise ValueError(
                f"Inline water {index} must contain oxygen[3] and "
                "hydrogen_offsets[2][3]"
            )
        if not np.all(np.isfinite(oxygen)) or not np.all(np.isfinite(offsets)):
            raise ValueError(f"Inline water {index} contains a non-finite coordinate")
        waters.append(WaterGeometry(oxygen, (offsets[0], offsets[1])))

    if not waters:
        raise ValueError("Solvation template contains no waters")
    rounded_oxygens = {tuple(np.round(water.oxygen, 6)) for water in waters}
    if len(rounded_oxygens) != len(waters):
        raise ValueError("Solvation template contains periodic duplicate oxygen positions")
    return WaterSolvationTemplate(
        mol_type=mol_type,
        residue_name=str(entry["residue_name"]),
        oxygen_atom=oxygen_name,
        hydrogen_atoms=(hydrogen_names[0], hydrogen_names[1]),
        cell_vectors=cell,
        waters=tuple(waters),
    )


def _is_water_identity(group: Optional[str], ff_resname: Optional[str], pdb_resname: str) -> bool:
    return bool(
        group in WATER_GROUPS
        or ff_resname == "WAT"
        or pdb_resname in {"WAT", "HOH"}
    )


def _extract_crystal_waters(molecule: Molecule) -> list[_RetainedWater]:
    grouped: Dict[Tuple[str, int, str], list[PDBAtomRecord]] = {}
    kept_passthrough: list[PDBAtomRecord] = []
    water_record_ids: set[int] = set()
    for record in molecule.passthrough_atoms:
        if _is_water_identity(record.group, record.ff_resname, record.resname):
            grouped.setdefault((record.chain_id, record.resseq, record.icode), []).append(record)
            water_record_ids.add(id(record))
        else:
            kept_passthrough.append(record)
    molecule.passthrough_atoms = kept_passthrough
    molecule.unassigned_records = [
        record for record in molecule.unassigned_records if id(record) not in water_record_ids
    ]

    # Also support callers that parsed W3/W4/W5 as canonical Molecule residues.
    for chain_id in list(molecule.chains):
        chain = molecule.chains[chain_id]
        kept_residues: list[Residue] = []
        for residue in chain.residues:
            if residue.group in WATER_GROUPS or residue.ff_resname == "WAT":
                records: list[PDBAtomRecord] = []
                for atom in residue.atoms.values():
                    if atom.coord is not None:
                        records.append(
                            PDBAtomRecord(
                                record_name="HETATM",
                                serial=atom.serial,
                                atom_name=atom.name,
                                altloc=atom.altloc,
                                resname=residue.original_resname or "WAT",
                                chain_id=residue.chain_id,
                                resseq=residue.resseq,
                                icode=residue.icode,
                                coord=atom.coord,
                                occupancy=atom.occupancy,
                                bfactor=atom.bfactor,
                                element=atom.element,
                                group=residue.group,
                                ff_resname=residue.ff_resname,
                            )
                        )
                grouped.setdefault((chain_id, residue.resseq, residue.icode), []).extend(records)
            else:
                kept_residues.append(residue)
        chain.residues = kept_residues
        for index, residue in enumerate(chain.residues):
            residue.index_in_chain = index
        if not chain.residues:
            del molecule.chains[chain_id]

    waters: list[_RetainedWater] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1], item[2])):
        records = grouped[key]
        atoms = {record.atom_name: np.asarray(record.coord, dtype=float) for record in records}
        oxygen_name = "O" if "O" in atoms else "OW" if "OW" in atoms else None
        if oxygen_name is None:
            continue
        canonical: Dict[str, np.ndarray] = {"O": atoms[oxygen_name]}
        hydrogen_coords = [
            coord for name, coord in atoms.items() if infer_element(name) == "H"
        ]
        for index, coord in enumerate(hydrogen_coords[:2], 1):
            canonical[f"H{index}"] = coord
        waters.append(
            _RetainedWater(
                records[0].resname if records else "HOH",
                canonical,
                frozenset(canonical),
            )
        )
    return waters


def _iter_polymer_atoms(molecule: Molecule) -> Iterator[Tuple[Residue, Atom]]:
    for chain in molecule.chains.values():
        for residue in chain.residues:
            for atom in residue.atoms.values():
                yield residue, atom


def _require_complete_and_collect_solute_coords(molecule: Molecule) -> list[np.ndarray]:
    coordinates: list[np.ndarray] = []
    missing: list[str] = []
    for residue, atom in _iter_polymer_atoms(molecule):
        if atom.coord is None:
            missing.append(f"{residue.chain_id}:{residue.resseq}{residue.icode} {residue.ff_resname}:{atom.name}")
        else:
            coordinates.append(np.asarray(atom.coord, dtype=float))
    for record in molecule.passthrough_atoms:
        coordinates.append(np.asarray(record.coord, dtype=float))
    if missing:
        preview = ", ".join(missing[:8])
        raise ValueError(
            f"Solvation requires a fully built solute; {len(missing)} atoms lack coordinates: {preview}"
        )
    if not coordinates:
        raise ValueError("Cannot solvate an empty molecule")
    return coordinates


def _deterministic_pca_rotation(coordinates: Sequence[np.ndarray]) -> np.ndarray:
    array = np.asarray(coordinates, dtype=float)
    if len(array) < 3:
        return np.eye(3)
    centered = array - array.mean(axis=0)
    covariance = centered.T @ centered / float(len(array))
    values, vectors = np.linalg.eigh(covariance)
    rotation = vectors[:, np.argsort(values)[::-1]]
    for column in range(3):
        vector = rotation[:, column]
        pivot = int(np.argmax(np.abs(vector)))
        if vector[pivot] < 0.0:
            rotation[:, column] *= -1.0
    if np.linalg.det(rotation) < 0.0:
        rotation[:, 2] *= -1.0
    return rotation


def _transform_coord(
    coord: np.ndarray, center: np.ndarray, rotation: np.ndarray, recenter: np.ndarray
) -> np.ndarray:
    return (coord - center) @ rotation - recenter


def _transform_molecule(
    molecule: Molecule, center: np.ndarray, rotation: np.ndarray, recenter: np.ndarray
) -> None:
    for _residue, atom in _iter_polymer_atoms(molecule):
        if atom.coord is not None:
            transformed = _transform_coord(np.asarray(atom.coord), center, rotation, recenter)
            atom.coord = tuple(float(value) for value in transformed)
    seen: set[int] = set()
    for collection in (molecule.passthrough_atoms, molecule.unassigned_records):
        for record in collection:
            if id(record) in seen:
                continue
            seen.add(id(record))
            transformed = _transform_coord(np.asarray(record.coord), center, rotation, recenter)
            record.coord = tuple(float(value) for value in transformed)


def _unit_truncated_octahedron_vectors() -> np.ndarray:
    # GROMACS reduced form; all vector lengths equal the periodic image distance d.
    return np.asarray(
        [
            (1.0, 0.0, 0.0),
            (1.0 / 3.0, 2.0 * math.sqrt(2.0) / 3.0, 0.0),
            (-1.0 / 3.0, math.sqrt(2.0) / 3.0, math.sqrt(6.0) / 3.0),
        ],
        dtype=float,
    )


def _ws_neighbor_vectors(unit_vectors: np.ndarray) -> list[np.ndarray]:
    translations: list[np.ndarray] = []
    for coeffs in product(range(-2, 3), repeat=3):
        if coeffs == (0, 0, 0):
            continue
        vector = np.asarray(coeffs, dtype=float) @ unit_vectors
        # The BCC Wigner-Seitz cell has 8 nearest and 6 next-nearest facets.
        if np.linalg.norm(vector) <= 2.0 / math.sqrt(3.0) + 1.0e-8:
            translations.append(vector)
    return translations


def _build_periodic_vectors(
    centered_points: np.ndarray, shape: str, padding: float
) -> np.ndarray:
    if shape in {"orthorhombic", "cubic"}:
        extent = np.ptp(centered_points, axis=0)
        lengths = extent + 2.0 * padding
        if shape == "cubic":
            lengths[:] = float(np.max(lengths))
        lengths = np.maximum(lengths, 1.0e-6)
        return np.diag(lengths)

    unit = _unit_truncated_octahedron_vectors()
    if shape == "truncated_octahedron":
        required = 0.0
        for translation in _ws_neighbor_vectors(unit):
            norm = float(np.linalg.norm(translation))
            unit_normal = translation / norm
            projection = float(np.max(centered_points @ unit_normal))
            required = max(required, 2.0 * (projection + padding) / norm)
        return unit * max(required, 1.0e-6)

    if shape == "truncated_octahedron_primitive":
        inverse = np.linalg.inv(unit)
        required = 0.0
        for axis in range(3):
            gradient = inverse[:, axis]
            coordinate = centered_points @ gradient
            required = max(
                required,
                2.0 * (
                    float(np.max(np.abs(coordinate)))
                    + padding * float(np.linalg.norm(gradient))
                ),
            )
        return unit * max(required, 1.0e-6)
    raise ValueError(f"Unsupported box shape {shape!r}")


def _inside_target_region(coord: np.ndarray, vectors: np.ndarray, shape: str) -> bool:
    tolerance = 1.0e-9
    if shape in {"orthorhombic", "cubic"}:
        half = np.diag(vectors) / 2.0
        return bool(np.all(coord >= -half - tolerance) and np.all(coord < half - tolerance))
    if shape == "truncated_octahedron_primitive":
        fractional = coord @ np.linalg.inv(vectors)
        return bool(np.all(fractional >= -0.5 - tolerance) and np.all(fractional < 0.5 - tolerance))
    if shape == "truncated_octahedron":
        unit = vectors / float(np.linalg.norm(vectors[0]))
        scale = float(np.linalg.norm(vectors[0]))
        for translation_unit in _ws_neighbor_vectors(unit):
            translation = translation_unit * scale
            if float(coord @ translation) > 0.5 * float(translation @ translation) + tolerance:
                return False
        return True
    raise ValueError(f"Unsupported box shape {shape!r}")


def _tile_template(
    template: WaterSolvationTemplate, vectors: np.ndarray, shape: str
) -> list[WaterGeometry]:
    # Vectorized replication in a safe Cartesian envelope.  The source box is
    # periodic, so each source oxygen is represented exactly once per tile.
    bound = 0.5 * np.sum(np.abs(vectors), axis=0)
    source_lengths = np.linalg.norm(template.cell_vectors, axis=1)
    if not np.allclose(template.cell_vectors, np.diag(np.diag(template.cell_vectors)), atol=1.0e-8):
        raise ValueError("Current W3 tiler expects an orthorhombic source water box")
    source_oxygens = np.asarray([water.oxygen for water in template.waters])
    source_min = np.min(source_oxygens, axis=0)
    source_max = np.max(source_oxygens, axis=0)
    minima = np.ceil((-bound - source_max) / source_lengths).astype(int)
    maxima = np.floor((bound - source_min) / source_lengths).astype(int)
    integer_tiles = np.asarray(
        list(
            product(
                range(minima[0], maxima[0] + 1),
                range(minima[1], maxima[1] + 1),
                range(minima[2], maxima[2] + 1),
            )
        ),
        dtype=float,
    )
    if integer_tiles.size == 0:
        return []
    translations = integer_tiles @ template.cell_vectors
    oxygens = (
        translations[:, None, :] + source_oxygens[None, :, :]
    ).reshape(-1, 3)

    tolerance = 1.0e-9
    if shape in {"orthorhombic", "cubic"}:
        half = np.diag(vectors) / 2.0
        mask = np.all(oxygens >= -half - tolerance, axis=1) & np.all(
            oxygens < half - tolerance, axis=1
        )
    elif shape == "truncated_octahedron_primitive":
        fractional = oxygens @ np.linalg.inv(vectors)
        mask = np.all(fractional >= -0.5 - tolerance, axis=1) & np.all(
            fractional < 0.5 - tolerance, axis=1
        )
    elif shape == "truncated_octahedron":
        scale = float(np.linalg.norm(vectors[0]))
        unit = vectors / scale
        translations_ws = np.asarray(_ws_neighbor_vectors(unit)) * scale
        projections = oxygens @ translations_ws.T
        limits = 0.5 * np.sum(translations_ws * translations_ws, axis=1)
        mask = np.all(projections <= limits[None, :] + tolerance, axis=1)
    else:
        raise ValueError(f"Unsupported box shape {shape!r}")

    source_indices = np.tile(np.arange(len(template.waters)), len(translations))[mask]
    selected_oxygens = oxygens[mask]
    order = np.lexsort(
        (selected_oxygens[:, 2], selected_oxygens[:, 1], selected_oxygens[:, 0])
    )
    return [
        WaterGeometry(
            selected_oxygens[index],
            template.waters[int(source_indices[index])].hydrogen_offsets,
        )
        for index in order
    ]


def _complete_crystal_water_hydrogens(
    waters: list[_RetainedWater], candidates: Sequence[WaterGeometry]
) -> None:
    if not waters:
        return
    if not candidates:
        raise ValueError("Cannot orient retained crystal waters: target box contains no template waters")
    oxygen_coords = [water.oxygen for water in candidates]
    index = _SpatialHash(oxygen_coords, 4.0)
    for water in waters:
        if "H1" in water.atoms and "H2" in water.atoms:
            continue
        oxygen = water.atoms["O"]
        nearest_index: Optional[int] = None
        nearest_d2 = float("inf")
        for radius in (4.0, 8.0, 16.0):
            for candidate_index in index.candidates(oxygen, radius):
                delta = oxygen_coords[candidate_index] - oxygen
                distance2 = float(delta @ delta)
                if distance2 < nearest_d2:
                    nearest_d2 = distance2
                    nearest_index = candidate_index
            if nearest_index is not None:
                break
        if nearest_index is None:
            nearest_index = min(
                range(len(candidates)),
                key=lambda item: float(np.sum((oxygen_coords[item] - oxygen) ** 2)),
            )
        offsets = candidates[nearest_index].hydrogen_offsets
        water.atoms.setdefault("H1", oxygen + offsets[0])
        water.atoms.setdefault("H2", oxygen + offsets[1])


def _solute_environment(
    molecule: Molecule,
    retained_waters: Sequence[_RetainedWater],
    parameters: SolvationVdwParameters,
    water_sigma: float,
) -> tuple[list[np.ndarray], list[float], float]:
    coordinates: list[np.ndarray] = []
    cutoffs: list[float] = []
    for residue, atom in _iter_polymer_atoms(molecule):
        if atom.coord is None:
            continue
        if residue.group is None:
            raise KeyError(f"Residue {residue.ff_resname} has no mol_type for sigma lookup")
        sigma = parameters.sigma(residue.group, residue.ff_resname, atom.name)
        coordinates.append(np.asarray(atom.coord, dtype=float))
        cutoffs.append(5.0 * (water_sigma + sigma))  # nm sigma -> Angstrom radius sum
    for record in molecule.passthrough_atoms:
        if record.group is None or record.ff_resname is None:
            raise KeyError(
                f"Passthrough atom {record.resname}:{record.atom_name} at "
                f"{record.chain_id}:{record.resseq}{record.icode} lacks converting identity"
            )
        sigma = parameters.sigma(record.group, record.ff_resname, record.atom_name)
        coordinates.append(np.asarray(record.coord, dtype=float))
        cutoffs.append(5.0 * (water_sigma + sigma))
    water_cutoff = 10.0 * water_sigma
    for water in retained_waters:
        coordinates.append(water.atoms["O"])
        cutoffs.append(water_cutoff)
    maximum = max(cutoffs, default=water_cutoff)
    return coordinates, cutoffs, maximum


def _remove_solute_clashes(
    candidates: Sequence[WaterGeometry],
    environment_coords: Sequence[np.ndarray],
    environment_cutoffs: Sequence[float],
    maximum_cutoff: float,
) -> tuple[list[WaterGeometry], int, int]:
    if not environment_coords:
        return list(candidates), 0, 0
    cell_size = max(maximum_cutoff, 1.0)
    index = _SpatialHash(environment_coords, cell_size)
    kept: list[WaterGeometry] = []
    removed = 0
    pair_checks = 0
    for water in candidates:
        collides = False
        for partner_index in index.candidates(water.oxygen, maximum_cutoff):
            pair_checks += 1
            delta = water.oxygen - environment_coords[partner_index]
            cutoff = environment_cutoffs[partner_index]
            if float(delta @ delta) < cutoff * cutoff:
                collides = True
                break
        if collides:
            removed += 1
        else:
            kept.append(water)
    return kept, removed, pair_checks


def _choose_solvent_chain_ids(
    molecule: Molecule, requested: Optional[str], water_count: int
) -> Tuple[str, ...]:
    required = max(1, int(math.ceil(water_count / 9999.0)))
    candidates = list("WXYZSVUTABCDEFGHIJKLMNOPQRabcdefghijklmnopqrstuvwxyz0123456789")
    if requested is not None:
        if requested in molecule.chains:
            raise ValueError(f"Requested solvent chain {requested!r} already exists")
        if requested in candidates:
            candidates.remove(requested)
        candidates.insert(0, requested)
    available = [candidate for candidate in candidates if candidate not in molecule.chains]
    if len(available) < required:
        raise ValueError(
            f"Need {required} unused PDB chain identifiers for {water_count} waters, "
            f"but only {len(available)} remain"
        )
    return tuple(available[:required])


def _append_waters(
    molecule: Molecule,
    retained: Sequence[_RetainedWater],
    generated: Sequence[WaterGeometry],
    template: WaterSolvationTemplate,
    chain_ids: Sequence[str],
    final_shift: np.ndarray,
) -> None:
    residues: list[Residue] = []
    first_chain_id = chain_ids[0]

    def make_atom(name: str, coord: np.ndarray, source: str) -> Atom:
        shifted = coord + final_shift
        return Atom(
            name=name,
            element=infer_element(name),
            coord=tuple(float(value) for value in shifted),
            built=source != "input",
            build_source=source,
        )

    for water in retained:
        atoms = {
            "O": make_atom("O", water.atoms["O"], "input"),
            "H1": make_atom(
                "H1", water.atoms["H1"],
                "input" if "H1" in water.original_atom_names else "crystal_water_orientation",
            ),
            "H2": make_atom(
                "H2", water.atoms["H2"],
                "input" if "H2" in water.original_atom_names else "crystal_water_orientation",
            ),
        }
        residues.append(
            Residue(
                chain_id=first_chain_id,
                resseq=len(residues) + 1,
                icode="",
                ff_resname=template.residue_name,
                atoms=atoms,
                index_in_chain=len(residues),
                original_resname=water.original_resname,
                group="W3",
            )
        )
    for water in generated:
        atoms = {
            template.oxygen_atom: make_atom(
                template.oxygen_atom, water.oxygen, "solvation_template"
            ),
            template.hydrogen_atoms[0]: make_atom(
                template.hydrogen_atoms[0],
                water.oxygen + water.hydrogen_offsets[0],
                "solvation_template",
            ),
            template.hydrogen_atoms[1]: make_atom(
                template.hydrogen_atoms[1],
                water.oxygen + water.hydrogen_offsets[1],
                "solvation_template",
            ),
        }
        residues.append(
            Residue(
                chain_id=first_chain_id,
                resseq=len(residues) + 1,
                icode="",
                ff_resname=template.residue_name,
                atoms=atoms,
                index_in_chain=len(residues),
                original_resname="WAT",
                group="W3",
            )
        )
    for chunk_index, chain_id in enumerate(chain_ids):
        chunk = residues[chunk_index * 9999 : (chunk_index + 1) * 9999]
        for residue_index, residue in enumerate(chunk):
            residue.chain_id = chain_id
            residue.resseq = residue_index + 1
            residue.index_in_chain = residue_index
        molecule.chains[chain_id] = Chain(chain_id=chain_id, residues=chunk)


def _translate_molecule(molecule: Molecule, shift: np.ndarray) -> None:
    for _residue, atom in _iter_polymer_atoms(molecule):
        if atom.coord is not None:
            coord = np.asarray(atom.coord, dtype=float) + shift
            atom.coord = tuple(float(value) for value in coord)
    seen: set[int] = set()
    for collection in (molecule.passthrough_atoms, molecule.unassigned_records):
        for record in collection:
            if id(record) in seen:
                continue
            seen.add(id(record))
            coord = np.asarray(record.coord, dtype=float) + shift
            record.coord = tuple(float(value) for value in coord)


def solvate_molecule(
    molecule: Molecule,
    template: WaterSolvationTemplate,
    parameters: SolvationVdwParameters,
    settings: Optional[SolvationSettings] = None,
    *,
    modify_myself: bool = False,
) -> tuple[Molecule, SolvationReport]:
    """Solvate a fully built molecule and return `(molecule, report)`.

    By default both input coordinates and existing crystal waters remain
    untouched because the operation is performed on a deep copy.
    """
    started = perf_counter()
    settings = settings or SolvationSettings()
    settings.validate()
    target = molecule if modify_myself else copy.deepcopy(molecule)
    target.periodic_box = None

    t0 = perf_counter()
    crystal_waters = _extract_crystal_waters(target)
    input_crystal_count = len(crystal_waters)
    retained = crystal_waters if settings.keep_crystal_waters else []
    solute_coords = _require_complete_and_collect_solute_coords(target)
    timings = {"prepare": perf_counter() - t0}

    t0 = perf_counter()
    center = np.mean(np.asarray(solute_coords), axis=0)
    rotation = (
        _deterministic_pca_rotation(solute_coords)
        if settings.rotate_principal_axes
        else np.eye(3)
    )
    rotated_solute = np.asarray([(coord - center) @ rotation for coord in solute_coords])
    rotated_retained_oxygen = [
        (water.atoms["O"] - center) @ rotation for water in retained
    ]
    box_points = (
        np.vstack([rotated_solute, np.asarray(rotated_retained_oxygen)])
        if rotated_retained_oxygen
        else rotated_solute
    )
    recenter = 0.5 * (np.min(box_points, axis=0) + np.max(box_points, axis=0))
    centered_points = box_points - recenter
    _transform_molecule(target, center, rotation, recenter)
    for water in retained:
        for atom_name in list(water.atoms):
            water.atoms[atom_name] = _transform_coord(
                water.atoms[atom_name], center, rotation, recenter
            )
    timings["pca_and_center"] = perf_counter() - t0

    t0 = perf_counter()
    vectors = _build_periodic_vectors(
        centered_points, settings.box_shape, settings.padding_angstrom
    )
    final_shift = 0.5 * np.sum(vectors, axis=0)
    target.periodic_box = PeriodicBox(
        vectors=tuple(
            tuple(float(value) for value in vector) for vector in vectors
        ),  # type: ignore[arg-type]
        origin=(0.0, 0.0, 0.0),
        shape=settings.box_shape,
    )
    timings["box_geometry"] = perf_counter() - t0

    t0 = perf_counter()
    candidates = _tile_template(template, vectors, settings.box_shape)
    timings["water_tiling_and_clipping"] = perf_counter() - t0

    t0 = perf_counter()
    if retained:
        _complete_crystal_water_hydrogens(retained, candidates)
    timings["crystal_water_completion"] = perf_counter() - t0

    t0 = perf_counter()
    water_sigma = parameters.sigma("W3", template.residue_name, template.oxygen_atom)
    environment_coords, environment_cutoffs, maximum_cutoff = _solute_environment(
        target, retained, parameters, water_sigma
    )
    generated, removed_by_solute, spatial_pair_checks = _remove_solute_clashes(
        candidates, environment_coords, environment_cutoffs, maximum_cutoff
    )
    timings["spatial_index_and_clash_filter"] = perf_counter() - t0

    t0 = perf_counter()
    _translate_molecule(target, final_shift)
    solvent_chains = _choose_solvent_chain_ids(
        target, settings.solvent_chain_id, len(retained) + len(generated)
    )
    _append_waters(target, retained, generated, template, solvent_chains, final_shift)
    timings["append_and_finalize"] = perf_counter() - t0
    timings["total"] = perf_counter() - started

    report = SolvationReport(
        box_shape=settings.box_shape,
        padding_angstrom=settings.padding_angstrom,
        box_vectors=target.periodic_box.vectors,
        input_crystal_waters=input_crystal_count,
        retained_crystal_waters=len(retained),
        template_candidates_in_box=len(candidates),
        waters_removed_by_solute=removed_by_solute,
        generated_waters=len(generated),
        total_waters=len(retained) + len(generated),
        solvent_chain_id=solvent_chains[0],
        solvent_chain_ids=solvent_chains,
        environment_atoms=len(environment_coords),
        spatial_pair_checks=spatial_pair_checks,
        timings_seconds=timings,
    )
    return target, report


def solvation_report_text(report: SolvationReport) -> str:
    lines = [
        "FORGE solvation report",
        "=======================",
        f"Box shape: {report.box_shape}",
        f"Padding from atom centers: {report.padding_angstrom:.3f} A",
        "Box vectors (A):",
    ]
    lines.extend(
        "  " + " ".join(f"{value:12.6f}" for value in vector)
        for vector in report.box_vectors
    )
    lines.extend(
        [
            f"Input crystal waters: {report.input_crystal_waters}",
            f"Retained crystal waters: {report.retained_crystal_waters}",
            f"Template candidates inside box: {report.template_candidates_in_box}",
            f"Waters removed by solute sigma clash: {report.waters_removed_by_solute}",
            f"Generated waters: {report.generated_waters}",
            f"Total W3/WAT waters: {report.total_waters}",
            f"Solvent chain(s): {', '.join(report.solvent_chain_ids or (report.solvent_chain_id,))}",
            f"Solute/environment atoms indexed: {report.environment_atoms}",
            f"Candidate-environment pair checks: {report.spatial_pair_checks}",
            "Timings (s):",
        ]
    )
    lines.extend(
        f"  {name}: {seconds:.6f}"
        for name, seconds in report.timings_seconds.items()
    )
    return "\n".join(lines) + "\n"
