#!/usr/bin/env python3
"""Reusable local molecular-mechanics support for FORGE structure building.

This module deliberately implements a small, local MM model rather than a
general simulation engine.  It owns force-field parameter loading, switched
Lennard-Jones/Coulomb evaluation, spatial indexing, proper-torsion evaluation,
and deterministic periodic-DOF optimization.  Molecule-specific coordinate
construction remains in :mod:`forge_molecule_builder` and is supplied through
callbacks, which keeps this module independent of build-plan classes.

Units follow Gromacs/UMFFF conventions: nm, kJ/mol, elementary charge, and
degrees.  Coordinate callbacks used by the optimizer return Angstrom values;
the local-energy helpers explicitly convert them to nm.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Dict,
    Hashable,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

import numpy as np


COULOMB_KJ_MOL_NM = 138.935458


def normalize_periodic_angle(angle: float) -> float:
    """Normalize degrees to ``(-180, 180]`` deterministically."""

    out = ((float(angle) + 180.0) % 360.0) - 180.0
    return 180.0 if out == -180.0 else out


@dataclass(frozen=True)
class MMAtomParams:
    atom_type: str
    charge: float
    sigma: float  # nm, Gromacs convention
    epsilon: float  # kJ/mol


@dataclass(frozen=True)
class MMTorsionTerm:
    k: float  # kJ/mol
    phase: float  # degrees
    periodicity: int


@dataclass
class FreeRotorSearchSettings:
    grid_step_deg: float = 15.0
    preselection_radius_nm: float = 0.9
    switch_radius_nm: float = 0.6
    cutoff_radius_nm: float = 0.8
    include_torsions: bool = True
    include_lj: bool = True
    include_electrostatics: bool = True
    collapse_other_free_rotors: bool = True


@dataclass
class SidechainOptimizationSettings:
    """Search and local-MM settings for residue-local missing-DOF branches."""

    grid_step_deg: float = 15.0
    beam_width: int = 8
    refined_seed_count: int = 4
    cartesian_tolerance_angstrom: float = 0.1
    max_refinement_sweeps: int = 12
    preselection_radius_nm: float = 0.9
    switch_radius_nm: float = 0.6
    cutoff_radius_nm: float = 0.8
    include_torsions: bool = True
    include_lj: bool = True
    include_electrostatics: bool = True
    coupling_energy_threshold_kj_mol: float = 1.0
    coupling_sigma_factor: float = 1.2

    def validate(self) -> None:
        if self.grid_step_deg <= 0.0 or self.grid_step_deg > 180.0:
            raise ValueError("grid_step_deg must lie in (0, 180]")
        if self.beam_width < 1:
            raise ValueError("beam_width must be positive")
        if self.refined_seed_count < 1:
            raise ValueError("refined_seed_count must be positive")
        if self.cartesian_tolerance_angstrom <= 0.0:
            raise ValueError("cartesian_tolerance_angstrom must be positive")
        if self.max_refinement_sweeps < 1:
            raise ValueError("max_refinement_sweeps must be positive")
        if not (
            0.0 < self.switch_radius_nm < self.cutoff_radius_nm
            <= self.preselection_radius_nm
        ):
            raise ValueError(
                "Expected 0 < switch_radius_nm < cutoff_radius_nm <= "
                "preselection_radius_nm"
            )


@dataclass
class MMParameterProvider:
    """Thin force-field parameter layer shared by local FORGE MM tasks."""

    residues: Dict[str, Dict[str, MMAtomParams]]
    residue_bonds: Dict[str, List[Tuple[str, str]]]
    dihedraltypes: Dict[Tuple[str, str, str, str], List[MMTorsionTerm]]
    fudge_lj: float = 0.5
    fudge_qq: float = 0.8333

    @classmethod
    def from_files(
        cls,
        residue_lib_file: Any,
        nonbonded_file: Any,
        bonded_file: Any,
        force_field_file: Optional[Any] = None,
    ) -> "MMParameterProvider":
        residues_raw, bonds = _read_umfff_residue_lib(residue_lib_file)
        atomtypes = _read_umfff_nonbonded(nonbonded_file)
        dihedraltypes = _read_umfff_bonded_dihedrals(bonded_file)
        fudge_lj, fudge_qq = (
            _read_umfff_defaults(force_field_file)
            if force_field_file is not None
            else (0.5, 0.8333)
        )

        residues: Dict[str, Dict[str, MMAtomParams]] = {}
        for resname, atoms in residues_raw.items():
            residues[resname] = {}
            for atom_name, atom_type, charge in atoms:
                if atom_type not in atomtypes:
                    raise ValueError(
                        f"Atom type {atom_type!r} for {resname}:{atom_name} "
                        "missing in nonbonded parameters"
                    )
                sigma, epsilon = atomtypes[atom_type]
                residues[resname][atom_name] = MMAtomParams(
                    atom_type, charge, sigma, epsilon
                )
        return cls(
            residues=residues,
            residue_bonds=bonds,
            dihedraltypes=dihedraltypes,
            fudge_lj=fudge_lj,
            fudge_qq=fudge_qq,
        )

    @classmethod
    def from_ff_ida_object(
        cls,
        ff_obj: Any,
        force_field_file: Optional[Any] = None,
    ) -> "MMParameterProvider":
        residues: Dict[str, Dict[str, MMAtomParams]] = {}
        bonds: Dict[str, List[Tuple[str, str]]] = {}
        for resname, unit in ff_obj.units.items():
            residues[resname] = {}
            names = unit["atoms"]["name"]
            types = unit["atoms"]["type"]
            charges = unit["atoms"]["charge"]
            sigmas = unit["atoms"].get("R", unit["atoms"].get("sigma"))
            epsilons = unit["atoms"].get("eps", unit["atoms"].get("epsilon"))
            if sigmas is None or epsilons is None:
                raise ValueError(
                    "FF_IDA object does not contain nonbonded sigma/R and eps arrays"
                )
            for name, atype, charge, sigma, epsilon in zip(
                names, types, charges, sigmas, epsilons
            ):
                residues[resname][name] = MMAtomParams(
                    str(atype),
                    float(charge),
                    float(sigma),
                    float(epsilon),
                )
            bonds[resname] = [tuple(b[:2]) for b in unit.get("bonds", [])]

        dihedraltypes: Dict[
            Tuple[str, str, str, str], List[MMTorsionTerm]
        ] = {}
        for key, terms in ff_obj.b.get("dihedraltypes", {}).items():
            dihedraltypes[tuple(key)] = [
                MMTorsionTerm(float(term[0]), float(term[1]), int(float(term[2])))
                for term in terms
            ]
        fudge_lj, fudge_qq = (
            _read_umfff_defaults(force_field_file)
            if force_field_file is not None
            else (0.5, 0.8333)
        )
        return cls(
            residues=residues,
            residue_bonds=bonds,
            dihedraltypes=dihedraltypes,
            fudge_lj=fudge_lj,
            fudge_qq=fudge_qq,
        )

    def atom_params(self, residue: Any, atom_name: str) -> MMAtomParams:
        try:
            return self.residues[residue.ff_resname][atom_name]
        except KeyError as exc:
            raise KeyError(
                f"MM parameters missing for {residue.ff_resname}:{atom_name}"
            ) from exc

    def proper_dihedral_terms(
        self,
        atom_types: Tuple[str, str, str, str],
    ) -> List[MMTorsionTerm]:
        return self.dihedraltypes.get(canonical_proper_dihedral_key(atom_types), [])


def _iter_text_lines(source: Any) -> Iterable[str]:
    if hasattr(source, "read"):
        try:
            source.seek(0)
        except Exception:
            pass
        for line in source:
            yield line.decode("utf-8") if isinstance(line, bytes) else line
        return
    with open(source, "r", encoding="utf-8") as handle:
        yield from handle


def _strip_umfff_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def _read_umfff_residue_lib(
    path: Any,
) -> Tuple[
    Dict[str, List[Tuple[str, str, float]]],
    Dict[str, List[Tuple[str, str]]],
]:
    residues: Dict[str, List[Tuple[str, str, float]]] = {}
    bonds: Dict[str, List[Tuple[str, str]]] = {}
    current_residue: Optional[str] = None
    section: Optional[str] = None
    sub_sections = {"atoms", "bonds", "impropers"}
    for raw in _iter_text_lines(path):
        line = _strip_umfff_comment(raw)
        if not line:
            continue
        if line.startswith("["):
            fields = line.replace("[", " [ ").replace("]", " ] ").split()
            if len(fields) >= 3 and fields[0] == "[" and fields[2] == "]":
                name = fields[1]
                if name in sub_sections:
                    section = name
                else:
                    current_residue = name
                    section = None
                    residues.setdefault(current_residue, [])
                    bonds.setdefault(current_residue, [])
            continue
        if current_residue is None or section is None:
            continue
        fields = line.split()
        if section == "atoms":
            if len(fields) < 3:
                raise ValueError(f"Bad atoms line in residue lib: {raw.rstrip()}")
            residues[current_residue].append(
                (fields[0], fields[1], float(fields[2]))
            )
        elif section == "bonds":
            if len(fields) < 2:
                raise ValueError(f"Bad bonds line in residue lib: {raw.rstrip()}")
            bonds[current_residue].append((fields[0], fields[1]))
    return residues, bonds


def _read_umfff_nonbonded(path: Any) -> Dict[str, Tuple[float, float]]:
    atomtypes: Dict[str, Tuple[float, float]] = {}
    in_atomtypes = False
    for raw in _iter_text_lines(path):
        line = _strip_umfff_comment(raw)
        if not line:
            continue
        if line.startswith("["):
            fields = line.replace("[", " [ ").replace("]", " ] ").split()
            in_atomtypes = len(fields) >= 3 and fields[1] == "atomtypes"
            continue
        if not in_atomtypes:
            continue
        fields = line.split()
        if len(fields) < 7:
            raise ValueError(f"Bad atomtypes line in nonbonded file: {raw.rstrip()}")
        atomtypes[fields[0]] = (float(fields[5]), float(fields[6]))
    return atomtypes


def _read_umfff_bonded_dihedrals(
    path: Any,
) -> Dict[Tuple[str, str, str, str], List[MMTorsionTerm]]:
    dihedrals: Dict[Tuple[str, str, str, str], List[MMTorsionTerm]] = {}
    section: Optional[str] = None
    for raw in _iter_text_lines(path):
        line = _strip_umfff_comment(raw)
        if not line:
            continue
        if line.startswith("["):
            fields = line.replace("[", " [ ").replace("]", " ] ").split()
            section = fields[1] if len(fields) >= 3 else None
            continue
        if section != "dihedraltypes":
            continue
        fields = line.split()
        if len(fields) < 8:
            raise ValueError(f"Bad dihedraltypes line in bonded file: {raw.rstrip()}")
        if fields[4] != "9":
            continue
        atom_types = (fields[0], fields[1], fields[2], fields[3])
        key = canonical_proper_dihedral_key(atom_types)
        dihedrals.setdefault(key, []).append(
            MMTorsionTerm(float(fields[6]), float(fields[5]), int(float(fields[7])))
        )
    return dihedrals


def _read_umfff_defaults(path: Any) -> Tuple[float, float]:
    if path is None:
        return 0.5, 0.8333
    section: Optional[str] = None
    for raw in _iter_text_lines(path):
        line = _strip_umfff_comment(raw)
        if not line:
            continue
        if line.startswith("["):
            fields = line.replace("[", " [ ").replace("]", " ] ").split()
            section = fields[1] if len(fields) >= 3 else None
            continue
        if section == "defaults":
            fields = line.split()
            if len(fields) >= 5:
                return float(fields[3]), float(fields[4])
    return 0.5, 0.8333


def canonical_proper_dihedral_key(
    atom_types: Tuple[str, str, str, str],
) -> Tuple[str, str, str, str]:
    reverse = (atom_types[3], atom_types[2], atom_types[1], atom_types[0])
    return min(tuple(atom_types), reverse)


@dataclass(frozen=True)
class SpatialAtom:
    atom_id: Hashable
    coord_nm: np.ndarray
    params: MMAtomParams


class MMSpatialIndex:
    """Small immutable cell list for repeated local MM queries."""

    def __init__(
        self,
        atoms: Iterable[SpatialAtom],
        cell_size_nm: float,
    ) -> None:
        if cell_size_nm <= 0.0:
            raise ValueError("cell_size_nm must be positive")
        self.cell_size_nm = float(cell_size_nm)
        cells: Dict[Tuple[int, int, int], List[SpatialAtom]] = defaultdict(list)
        for atom in atoms:
            cells[self.cell_key(atom.coord_nm)].append(atom)
        self.cells = dict(cells)

    def cell_key(self, coord_nm: np.ndarray) -> Tuple[int, int, int]:
        return tuple(
            int(math.floor(float(value) / self.cell_size_nm))
            for value in coord_nm
        )  # type: ignore[return-value]

    def query(self, coord_nm: np.ndarray, radius_nm: float) -> List[SpatialAtom]:
        if radius_nm < 0.0:
            raise ValueError("radius_nm cannot be negative")
        center = self.cell_key(coord_nm)
        shell = int(math.ceil(radius_nm / self.cell_size_nm))
        radius_squared = radius_nm * radius_nm
        result: List[SpatialAtom] = []
        for dx in range(-shell, shell + 1):
            for dy in range(-shell, shell + 1):
                for dz in range(-shell, shell + 1):
                    for atom in self.cells.get(
                        (center[0] + dx, center[1] + dy, center[2] + dz),
                        (),
                    ):
                        delta = atom.coord_nm - coord_nm
                        if float(np.dot(delta, delta)) <= radius_squared:
                            result.append(atom)
        return result


def switch_weights(
    distances_nm: np.ndarray,
    switch_radius_nm: float,
    cutoff_radius_nm: float,
) -> np.ndarray:
    """Return cubic smoothstep weights, one below switch and zero at cutoff."""

    distances = np.asarray(distances_nm, dtype=float)
    weights = np.ones_like(distances)
    weights[distances >= cutoff_radius_nm] = 0.0
    active = (distances > switch_radius_nm) & (distances < cutoff_radius_nm)
    if np.any(active):
        x = (distances[active] - switch_radius_nm) / (
            cutoff_radius_nm - switch_radius_nm
        )
        weights[active] = 1.0 - (3.0 * x * x - 2.0 * x * x * x)
    return weights


def vectorized_nonbonded_energy(
    moving_coord_nm: np.ndarray,
    moving_params: MMAtomParams,
    partner_coords_nm: np.ndarray,
    partner_charges: np.ndarray,
    partner_sigmas: np.ndarray,
    partner_epsilons: np.ndarray,
    scale_lj: np.ndarray,
    scale_qq: np.ndarray,
    *,
    switch_radius_nm: float,
    cutoff_radius_nm: float,
    include_lj: bool = True,
    include_electrostatics: bool = True,
) -> float:
    """Evaluate switched nonbonded energy for one moving atom."""

    if partner_coords_nm.shape[0] == 0:
        return 0.0
    delta = partner_coords_nm - np.asarray(moving_coord_nm, dtype=float)
    r_squared = np.einsum("ij,ij->i", delta, delta)
    active = (r_squared < cutoff_radius_nm**2) & (r_squared > 1.0e-24)
    if not np.any(active):
        return 0.0
    r_nm = np.sqrt(r_squared[active])
    weights = switch_weights(r_nm, switch_radius_nm, cutoff_radius_nm)
    pair_energy = np.zeros_like(r_nm)
    if include_lj:
        sigma = 0.5 * (moving_params.sigma + partner_sigmas[active])
        epsilon = np.sqrt(
            np.maximum(0.0, moving_params.epsilon * partner_epsilons[active])
        )
        valid = (sigma > 0.0) & (epsilon > 0.0)
        if np.any(valid):
            sr6 = (sigma[valid] / r_nm[valid]) ** 6
            pair_energy[valid] += (
                scale_lj[active][valid]
                * 4.0
                * epsilon[valid]
                * (sr6 * sr6 - sr6)
            )
    if include_electrostatics and moving_params.charge != 0.0:
        pair_energy += (
            scale_qq[active]
            * COULOMB_KJ_MOL_NM
            * moving_params.charge
            * partner_charges[active]
            / r_nm
        )
    return float(np.sum(weights * pair_energy))


def pair_nonbonded_energy(
    coord_a_nm: np.ndarray,
    params_a: MMAtomParams,
    coord_b_nm: np.ndarray,
    params_b: MMAtomParams,
    *,
    scale_lj: float,
    scale_qq: float,
    switch_radius_nm: float,
    cutoff_radius_nm: float,
    include_lj: bool = True,
    include_electrostatics: bool = True,
) -> float:
    """Scalar pair wrapper used for moving-moving and coupling energies."""

    return vectorized_nonbonded_energy(
        np.asarray(coord_a_nm, dtype=float),
        params_a,
        np.asarray([coord_b_nm], dtype=float),
        np.asarray([params_b.charge], dtype=float),
        np.asarray([params_b.sigma], dtype=float),
        np.asarray([params_b.epsilon], dtype=float),
        np.asarray([scale_lj], dtype=float),
        np.asarray([scale_qq], dtype=float),
        switch_radius_nm=switch_radius_nm,
        cutoff_radius_nm=cutoff_radius_nm,
        include_lj=include_lj,
        include_electrostatics=include_electrostatics,
    )


def proper_torsion_energy(
    dihedral_degrees: float,
    terms: Sequence[MMTorsionTerm],
) -> float:
    return sum(
        term.k
        * (
            1.0
            + math.cos(
                math.radians(
                    term.periodicity * float(dihedral_degrees) - term.phase
                )
            )
        )
        for term in terms
    )


def topological_shells_upto_three(
    graph: Mapping[Hashable, Set[Hashable]],
    start: Hashable,
) -> Tuple[Set[Hashable], Set[Hashable], Set[Hashable]]:
    visited = {start}
    frontier = {start}
    shells: List[Set[Hashable]] = []
    for _depth in range(3):
        next_frontier: Set[Hashable] = set()
        for atom_id in frontier:
            next_frontier.update(graph.get(atom_id, set()))
        next_frontier.difference_update(visited)
        shells.append(next_frontier)
        visited.update(next_frontier)
        frontier = next_frontier
    return shells[0], shells[1], shells[2]


def topological_distance_upto(
    graph: Mapping[Hashable, Set[Hashable]],
    start: Hashable,
    target: Hashable,
    max_depth: int = 3,
) -> Optional[int]:
    """Return shortest graph distance up to ``max_depth``, otherwise ``None``."""

    if start == target:
        return 0
    visited = {start}
    frontier = {start}
    for depth in range(1, max_depth + 1):
        next_frontier: Set[Hashable] = set()
        for node in frontier:
            for neighbor in graph.get(node, set()):
                if neighbor == target:
                    return depth
                if neighbor not in visited:
                    visited.add(neighbor)
                    next_frontier.add(neighbor)
        frontier = next_frontier
        if not frontier:
            break
    return None


@dataclass(frozen=True)
class PeriodicOptimizationResult:
    values: Mapping[Hashable, float]
    energy_kj_mol: float
    evaluations: int
    coarse_candidates: int
    refinement_sweeps: int


def optimize_symmetric_periodic_1d(
    template_degrees: float,
    symmetry_order: int,
    grid_step_degrees: float,
    score: Callable[[float], float],
) -> float:
    """Grid-scan one symmetric rotor and interpolate its local minimum.

    This is the established free-rotor-H search policy.  It remains separate
    from the multi-DOF beam/pattern optimizer because equivalent H permutations
    reduce its physical search period to ``360 / symmetry_order``.
    """

    if symmetry_order < 1:
        raise ValueError("symmetry_order must be positive")
    if grid_step_degrees <= 0.0:
        raise ValueError("grid_step_degrees must be positive")
    search_period = 360.0 / symmetry_order
    n_grid = max(3, int(round(search_period / grid_step_degrees)))
    step = search_period / n_grid
    phis = [float(template_degrees) + index * step for index in range(n_grid)]
    energies = [float(score(normalize_periodic_angle(phi))) for phi in phis]
    minimum_index = min(
        range(n_grid),
        key=lambda index: (
            energies[index],
            abs(normalize_periodic_angle(phis[index] - template_degrees)),
        ),
    )
    previous_index = (minimum_index - 1) % n_grid
    next_index = (minimum_index + 1) % n_grid
    denominator = (
        energies[previous_index]
        - 2.0 * energies[minimum_index]
        + energies[next_index]
    )
    refined = phis[minimum_index]
    if denominator > 1.0e-12 and math.isfinite(denominator):
        delta = (
            0.5
            * step
            * (energies[previous_index] - energies[next_index])
            / denominator
        )
        if math.isfinite(delta) and abs(delta) <= step:
            refined += delta
    return normalize_periodic_angle(refined)


def _value_tuple(
    dof_order: Sequence[Hashable],
    values: Mapping[Hashable, float],
) -> Tuple[float, ...]:
    return tuple(normalize_periodic_angle(values[key]) for key in dof_order)


def _periodic_distance(a: float, b: float) -> float:
    return abs(normalize_periodic_angle(float(a) - float(b)))


def _refine_periodic_seeds(
    order: Sequence[Hashable],
    initial: Mapping[Hashable, float],
    seeds: Sequence[Mapping[Hashable, float]],
    evaluate: Callable[[Mapping[Hashable, float]], float],
    coordinates: Callable[
        [Mapping[Hashable, float]], Mapping[Hashable, Sequence[float]]
    ],
    affected_atoms_by_dof: Mapping[Hashable, Sequence[Hashable]],
    settings: SidechainOptimizationSettings,
    initial_step_degrees: float,
) -> Tuple[Dict[Hashable, float], float, int]:
    """Pattern-refine one or more periodic seeds with a shared evaluator."""

    def max_displacement(
        state: Mapping[Hashable, float],
        dof: Hashable,
        delta: float,
    ) -> float:
        altered = dict(state)
        altered[dof] = normalize_periodic_angle(altered[dof] + delta)
        before = coordinates(state)
        after = coordinates(altered)
        affected = affected_atoms_by_dof.get(dof, ())
        if not affected:
            raise ValueError(f"Optimized DOF {dof!r} has no affected atoms")
        missing = [
            atom_id
            for atom_id in affected
            if atom_id not in before or atom_id not in after
        ]
        if missing:
            raise ValueError(
                f"Coordinate callback omitted atoms affected by {dof!r}: {missing!r}"
            )
        displacements = [
            float(
                np.linalg.norm(
                    np.asarray(after[atom_id], dtype=float)
                    - np.asarray(before[atom_id], dtype=float)
                )
            )
            for atom_id in affected
        ]
        return max(displacements)

    if initial_step_degrees <= 0.0:
        raise ValueError("initial_step_degrees must be positive")
    if not seeds:
        raise ValueError("At least one refinement seed is required")

    best_state = {
        dof: normalize_periodic_angle(seeds[0][dof]) for dof in order
    }
    best_energy = evaluate(best_state)
    total_sweeps = 0
    for seed in seeds:
        state = {dof: normalize_periodic_angle(seed[dof]) for dof in order}
        energy = evaluate(state)
        steps = {dof: float(initial_step_degrees) for dof in order}
        sweeps = 0
        while sweeps < settings.max_refinement_sweeps:
            sweeps += 1
            changed = False
            all_converged = True
            for dof in order:
                step_size = steps[dof]
                if max_displacement(state, dof, step_size) < (
                    settings.cartesian_tolerance_angstrom
                ):
                    continue
                all_converged = False
                candidates: List[Tuple[float, Dict[Hashable, float]]] = []
                for sign in (-1.0, 1.0):
                    candidate = dict(state)
                    candidate[dof] = normalize_periodic_angle(
                        candidate[dof] + sign * step_size
                    )
                    candidates.append((evaluate(candidate), candidate))
                candidate_energy, candidate_state = min(
                    candidates,
                    key=lambda item: (
                        item[0],
                        _periodic_distance(item[1][dof], initial[dof]),
                        _value_tuple(order, item[1]),
                    ),
                )
                if candidate_energy < energy - 1.0e-10:
                    state = candidate_state
                    energy = candidate_energy
                    changed = True
                else:
                    steps[dof] *= 0.5
            if all_converged:
                break
            if not changed and all(
                max_displacement(state, dof, steps[dof])
                < settings.cartesian_tolerance_angstrom
                for dof in order
            ):
                break
        total_sweeps += sweeps
        if (energy, _value_tuple(order, state)) < (
            best_energy,
            _value_tuple(order, best_state),
        ):
            best_state = state
            best_energy = energy
    return best_state, best_energy, total_sweeps


def refine_periodic_dofs(
    dof_order: Sequence[Hashable],
    initial_values: Mapping[Hashable, float],
    score: Callable[[Mapping[Hashable, float]], float],
    coordinates: Callable[
        [Mapping[Hashable, float]], Mapping[Hashable, Sequence[float]]
    ],
    affected_atoms_by_dof: Mapping[Hashable, Sequence[Hashable]],
    settings: SidechainOptimizationSettings,
    *,
    initial_step_degrees: Optional[float] = None,
) -> PeriodicOptimizationResult:
    """Locally refine periodic DOFs from exactly one supplied configuration.

    Unlike :func:`optimize_periodic_dofs`, this function performs no coarse
    grid or beam search.  It is intended for interactive relaxation of a GUI
    configuration into the local basin containing that configuration.
    """

    settings.validate()
    order = tuple(dof_order)
    if not order:
        return PeriodicOptimizationResult({}, 0.0, 0, 0, 0)
    if set(initial_values) != set(order):
        raise ValueError("initial_values must contain exactly every optimized DOF")
    initial = {
        dof: normalize_periodic_angle(initial_values[dof]) for dof in order
    }
    cache: Dict[Tuple[float, ...], float] = {}

    def evaluate(values: Mapping[Hashable, float]) -> float:
        clean = {
            dof: normalize_periodic_angle(values[dof]) for dof in order
        }
        key = tuple(round(clean[dof], 10) for dof in order)
        if key not in cache:
            energy = float(score(clean))
            cache[key] = energy if math.isfinite(energy) else math.inf
        return cache[key]

    n_grid = max(3, int(round(360.0 / settings.grid_step_deg)))
    grid_step = 360.0 / n_grid
    best_state, best_energy, total_sweeps = _refine_periodic_seeds(
        order,
        initial,
        (initial,),
        evaluate,
        coordinates,
        affected_atoms_by_dof,
        settings,
        0.5 * grid_step
        if initial_step_degrees is None
        else float(initial_step_degrees),
    )
    return PeriodicOptimizationResult(
        values=best_state,
        energy_kj_mol=best_energy,
        evaluations=len(cache),
        coarse_candidates=0,
        refinement_sweeps=total_sweeps,
    )


def optimize_periodic_dofs(
    dof_order: Sequence[Hashable],
    initial_values: Mapping[Hashable, float],
    score: Callable[[Mapping[Hashable, float]], float],
    coordinates: Callable[
        [Mapping[Hashable, float]], Mapping[Hashable, Sequence[float]]
    ],
    affected_atoms_by_dof: Mapping[Hashable, Sequence[Hashable]],
    settings: SidechainOptimizationSettings,
) -> PeriodicOptimizationResult:
    """Deterministically optimize periodic DOFs with beam scan and refinement.

    One-dimensional problems are exhaustively scanned.  Multi-dimensional
    problems use deterministic forward and reverse beam scans, avoiding the
    exponential full Cartesian product.  Local coordinate refinement stops
    when a trial DOF perturbation moves every affected atom by less than the
    requested Cartesian tolerance.
    """

    settings.validate()
    order = tuple(dof_order)
    if not order:
        return PeriodicOptimizationResult({}, 0.0, 0, 0, 0)
    if set(initial_values) != set(order):
        raise ValueError("initial_values must contain exactly every optimized DOF")

    cache: Dict[Tuple[float, ...], float] = {}

    def normalized(values: Mapping[Hashable, float]) -> Dict[Hashable, float]:
        return {key: normalize_periodic_angle(values[key]) for key in order}

    def evaluate(values: Mapping[Hashable, float]) -> float:
        clean = normalized(values)
        key = tuple(round(clean[dof], 10) for dof in order)
        if key not in cache:
            energy = float(score(clean))
            cache[key] = energy if math.isfinite(energy) else math.inf
        return cache[key]

    n_grid = max(3, int(round(360.0 / settings.grid_step_deg)))
    grid_step = 360.0 / n_grid
    initial = normalized(initial_values)

    def beam_pass(scan_order: Sequence[Hashable]) -> List[Dict[Hashable, float]]:
        beam: List[Dict[Hashable, float]] = [dict(initial)]
        for dof in scan_order:
            generated: Dict[Tuple[float, ...], Dict[Hashable, float]] = {}
            base = initial[dof]
            for state in beam:
                for index in range(n_grid):
                    candidate = dict(state)
                    candidate[dof] = normalize_periodic_angle(base + index * grid_step)
                    generated[_value_tuple(order, candidate)] = candidate
            ranked = sorted(
                generated.values(),
                key=lambda state: (
                    evaluate(state),
                    sum(
                        _periodic_distance(state[key], initial[key])
                        for key in order
                    ),
                    _value_tuple(order, state),
                ),
            )
            beam = ranked[: settings.beam_width]
        return beam

    seeds_by_key: Dict[Tuple[float, ...], Dict[Hashable, float]] = {
        _value_tuple(order, initial): dict(initial)
    }
    if len(order) == 1:
        for candidate in beam_pass(order):
            seeds_by_key[_value_tuple(order, candidate)] = candidate
        # beam_width must not truncate distinct one-dimensional minima before
        # ranking; add the complete grid explicitly.
        dof = order[0]
        for index in range(n_grid):
            candidate = dict(initial)
            candidate[dof] = normalize_periodic_angle(initial[dof] + index * grid_step)
            seeds_by_key[_value_tuple(order, candidate)] = candidate
    else:
        for scan_order in (order, tuple(reversed(order))):
            for candidate in beam_pass(scan_order):
                seeds_by_key[_value_tuple(order, candidate)] = candidate

    ranked_seeds = sorted(
        seeds_by_key.values(),
        key=lambda state: (
            evaluate(state),
            sum(_periodic_distance(state[key], initial[key]) for key in order),
            _value_tuple(order, state),
        ),
    )
    seeds = ranked_seeds[: settings.refined_seed_count]
    coarse_candidates = len(seeds_by_key)

    best_state, best_energy, total_sweeps = _refine_periodic_seeds(
        order,
        initial,
        seeds,
        evaluate,
        coordinates,
        affected_atoms_by_dof,
        settings,
        0.5 * grid_step,
    )

    return PeriodicOptimizationResult(
        values=normalized(best_state),
        energy_kj_mol=best_energy,
        evaluations=len(cache),
        coarse_candidates=coarse_candidates,
        refinement_sweeps=total_sweeps,
    )
