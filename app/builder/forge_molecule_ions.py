"""Crystal-ion cleanup and post-solvation ion placement for FORGE molecules."""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from itertools import combinations, product
import math
from time import perf_counter
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np

from forge_molecule_builder import AtomID, build_mm_bond_graph
from forge_molecule_mm import topological_distance_upto
from forge_molecule_parser import Atom, Chain, Molecule, PDBAtomRecord, Residue, infer_element
from forge_molecule_solvation import SolvationVdwParameters


ION_MOL_TYPES = frozenset(("I1", "I1+", "Im", "Im+"))
WATER_MOL_TYPES = frozenset(("W3", "W4", "W5"))
_IMAGE_OFFSETS = np.asarray(list(product((-1.0, 0.0, 1.0), repeat=3)))


@dataclass(frozen=True)
class IonIdentity:
    mol_type: str
    resname: str


@dataclass(frozen=True)
class SaltSpecification:
    cation: IonIdentity
    anion: IonIdentity
    concentration_mol_l: float


@dataclass
class IonPlacementSettings:
    clean_crystal_ions: bool = True
    replace_structural_multivalent_with_mg: bool = False
    concentration_mode: str = "water_ratio"
    inner_shell_sigma_factor: float = 1.2
    monovalent_cation_inner_shell_sigma_factor: float = 1.0
    potential_switch_angstrom: float = 10.0
    potential_cutoff_angstrom: float = 12.0
    minimum_ion_separation_angstrom: float = 3.0
    salt_random_seed: int = 0
    default_cation: IonIdentity = IonIdentity("I1", "K+")
    default_anion: IonIdentity = IonIdentity("I1", "Cl-")
    magnesium_identity: IonIdentity = IonIdentity("Im", "Mg2+")

    def validate(self) -> None:
        if self.concentration_mode not in {"water_ratio", "box_volume"}:
            raise ValueError("concentration_mode must be 'water_ratio' or 'box_volume'")
        if self.inner_shell_sigma_factor <= 0.0:
            raise ValueError("inner_shell_sigma_factor must be positive")
        if self.monovalent_cation_inner_shell_sigma_factor <= 0.0:
            raise ValueError(
                "monovalent_cation_inner_shell_sigma_factor must be positive"
            )
        if not 0.0 <= self.potential_switch_angstrom < self.potential_cutoff_angstrom:
            raise ValueError("Potential radii must satisfy 0 <= switch < cutoff")
        if self.minimum_ion_separation_angstrom < 0.0:
            raise ValueError("minimum_ion_separation_angstrom must be non-negative")
        if not isinstance(self.salt_random_seed, int):
            raise ValueError("salt_random_seed must be an integer")


@dataclass
class CrystalIonCleanupReport:
    input_ions: int = 0
    removed_monovalent: int = 0
    removed_nonstructural_multivalent: int = 0
    retained_structural_monovalent: int = 0
    retained_structural_multivalent: int = 0
    replaced_by_magnesium: int = 0
    structural_contacts: Dict[str, int] = field(default_factory=dict)


@dataclass
class IonAdditionReport:
    initial_fixed_charge: float
    final_system_charge: float
    initial_waters: int
    final_waters: int
    neutralization_ions: Dict[str, int]
    salt_formula_units: list[Dict[str, Any]]
    added_ions: Dict[str, int]
    ion_chain_ids: Tuple[str, ...]
    potential_pair_evaluations: int
    neutralization_candidate_waters: int
    salt_random_seed: int
    minimum_added_ion_distance_angstrom: Optional[float]
    added_ion_pairs_below_3_angstrom: int
    timings_seconds: Dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class _FixedAtom:
    coord: np.ndarray
    charge: float
    sigma: float


@dataclass
class _WaterCandidate:
    chain_id: str
    residue: Residue
    coord: np.ndarray
    active: bool = True


@dataclass(frozen=True)
class _PlacedIon:
    identity: IonIdentity
    coord: np.ndarray
    charge: float
    sigma: float
    source: str
    candidate_index: Optional[int] = None


def load_salt_specifications(data: Mapping[str, Any]) -> list[SaltSpecification]:
    salts: list[SaltSpecification] = []
    for index, raw in enumerate(data.get("salts", [])):
        try:
            cation = IonIdentity(str(raw["cation"]["mol_type"]), str(raw["cation"]["resname"]))
            anion = IonIdentity(str(raw["anion"]["mol_type"]), str(raw["anion"]["resname"]))
            concentration = float(raw["concentration"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Malformed salt specification at index {index}") from exc
        if not math.isfinite(concentration) or concentration < 0.0:
            raise ValueError(f"Salt concentration at index {index} must be finite and non-negative")
        salts.append(SaltSpecification(cation, anion, concentration))
    return salts


def _ion_params(parameters: SolvationVdwParameters, identity: IonIdentity) -> tuple[str, float, float]:
    if identity.mol_type not in ION_MOL_TYPES:
        raise ValueError(f"{identity} is not an ion mol_type")
    try:
        residue = parameters.providers[identity.mol_type].residues[identity.resname]
    except KeyError as exc:
        raise KeyError(f"Ion parameters missing for {identity.mol_type}:{identity.resname}") from exc
    if len(residue) != 1:
        raise ValueError(f"Ion {identity} must contain exactly one force-field atom")
    atom_name, params = next(iter(residue.items()))
    return atom_name, float(params.charge), float(params.sigma)


def _ion_element(identity: IonIdentity) -> str:
    letters = "".join(character for character in identity.resname if character.isalpha())
    if not letters:
        raise ValueError(f"Cannot infer chemical element for ion {identity}")
    return letters[0].upper() if len(letters) == 1 else letters[:2].title()


def _integer_charge(charge: float, label: str, tolerance: float = 0.05) -> int:
    rounded = int(round(charge))
    if abs(charge - rounded) > tolerance:
        raise ValueError(f"{label} charge {charge:.6f} is not sufficiently close to an integer")
    return rounded


def _record_is_ion(record: PDBAtomRecord) -> bool:
    return bool(record.group in ION_MOL_TYPES)


def _acceptor_atoms(
    molecule: Molecule, building_template: Mapping[str, Any]
) -> list[tuple[AtomID, Residue, Atom]]:
    result: list[tuple[AtomID, Residue, Atom]] = []
    for chain_id, chain in molecule.chains.items():
        for residue in chain.residues:
            sites = (
                building_template.get(residue.group or "", {})
                .get(residue.ff_resname, {})
                .get("hbond_sites", {})
            )
            for atom_name, site in sites.items():
                atom = residue.atoms.get(atom_name)
                if site.get("acceptor") is True and atom is not None and atom.coord is not None:
                    result.append(
                        (AtomID(chain_id, residue.index_in_chain, atom_name), residue, atom)
                    )
    return result


def _independent_contact_count(
    contacts: Sequence[AtomID], graph: Mapping[AtomID, set[AtomID]]
) -> int:
    # Maximum mutually independent subset.  This avoids over-collapsing a
    # transitive A--B--C relation when A and C themselves are topologically
    # independent. Coordination shells are tiny, so exact enumeration is cheap.
    unique = tuple(dict.fromkeys(contacts))
    for size in range(len(unique), 0, -1):
        for subset in combinations(unique, size):
            if all(
                topological_distance_upto(graph, left, right, max_depth=2)
                is None
                for left, right in combinations(subset, 2)
            ):
                return size
    return 0


def _add_explicit_covalent_edges(
    molecule: Molecule, graph: Dict[AtomID, set[AtomID]]
) -> None:
    lookup: Dict[tuple[str, int, str, str], AtomID] = {}
    for chain_id, chain in molecule.chains.items():
        for residue in chain.residues:
            for atom_name in residue.atoms:
                lookup[(chain_id, residue.resseq, residue.icode, atom_name)] = AtomID(
                    chain_id, residue.index_in_chain, atom_name
                )
    for left_key, right_key in molecule.covalent_bonds:
        left, right = lookup.get(left_key), lookup.get(right_key)
        if left is None or right is None:
            continue
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)


def clean_crystal_ions(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    parameters: SolvationVdwParameters,
    settings: Optional[IonPlacementSettings] = None,
    *,
    modify_myself: bool = False,
) -> tuple[Molecule, CrystalIonCleanupReport]:
    """Remove diffuse crystal ions and retain only independently multi-coordinated multivalents."""
    settings = settings or IonPlacementSettings()
    settings.validate()
    target = molecule if modify_myself else copy.deepcopy(molecule)
    report = CrystalIonCleanupReport()
    ion_records = [record for record in target.passthrough_atoms if _record_is_ion(record)]
    report.input_ions = len(ion_records)
    if not settings.clean_crystal_ions:
        return target, report

    acceptors = _acceptor_atoms(target, building_template)
    builder_params = parameters.merged_builder_provider()
    graph = build_mm_bond_graph(target, builder_params)
    _add_explicit_covalent_edges(target, graph)
    keep_ids: set[int] = set()
    replacements: Dict[int, tuple[str, float, float]] = {}
    for record in ion_records:
        identity = IonIdentity(str(record.group), str(record.ff_resname))
        _atom_name, charge, ion_sigma = _ion_params(parameters, identity)
        formal_charge = _integer_charge(charge, f"{identity.mol_type}:{identity.resname}")
        # The ligand model below is intentionally cation-specific: it searches
        # atoms marked as hydrogen-bond acceptors.  Applying it to monovalent
        # anions would retain opportunistic contacts to the wrong chemical
        # environment.  Structural-anion recognition needs a separate donor-
        # oriented geometric model and is therefore not attempted here.
        if formal_charge == -1:
            report.removed_monovalent += 1
            continue
        coord = np.asarray(record.coord, dtype=float)
        contacts: list[AtomID] = []
        sigma_factor = (
            settings.monovalent_cation_inner_shell_sigma_factor
            if formal_charge == 1
            else settings.inner_shell_sigma_factor
        )
        for atom_id, residue, atom in acceptors:
            ligand_sigma = parameters.sigma(
                str(residue.group), residue.ff_resname, atom.name
            )
            cutoff = sigma_factor * 5.0 * (
                ion_sigma + ligand_sigma
            )
            if float(np.linalg.norm(coord - np.asarray(atom.coord))) <= cutoff:
                contacts.append(atom_id)
        independent = _independent_contact_count(contacts, graph)
        label = f"{record.chain_id}:{record.resseq}{record.icode}:{identity.resname}"
        report.structural_contacts[label] = independent
        if independent < 2:
            if abs(formal_charge) <= 1:
                report.removed_monovalent += 1
            else:
                report.removed_nonstructural_multivalent += 1
            continue
        keep_ids.add(id(record))
        if abs(formal_charge) <= 1:
            report.retained_structural_monovalent += 1
        else:
            report.retained_structural_multivalent += 1
        if (
            abs(formal_charge) > 1
            and settings.replace_structural_multivalent_with_mg
            and identity != settings.magnesium_identity
        ):
            replacements[id(record)] = _ion_params(parameters, settings.magnesium_identity)
            record.group = settings.magnesium_identity.mol_type
            record.ff_resname = settings.magnesium_identity.resname
            record.resname = settings.magnesium_identity.resname
            record.atom_name = replacements[id(record)][0]
            record.element = _ion_element(settings.magnesium_identity)
            report.replaced_by_magnesium += 1

    target.passthrough_atoms = [
        record
        for record in target.passthrough_atoms
        if not _record_is_ion(record) or id(record) in keep_ids
    ]
    target.unassigned_records = [
        record
        for record in target.unassigned_records
        if not _record_is_ion(record) or id(record) in keep_ids
    ]
    return target, report


class _PeriodicCellList:
    def __init__(self, coordinates: Sequence[np.ndarray], vectors: np.ndarray, cutoff: float):
        self.coords = np.asarray(coordinates, dtype=float)
        self.vectors = np.asarray(vectors, dtype=float)
        self.inverse = np.linalg.inv(self.vectors)
        spacings = np.asarray(
            [1.0 / np.linalg.norm(self.inverse[:, axis]) for axis in range(3)]
        )
        self.spacings = spacings
        self.counts = np.maximum(
            1, np.floor(spacings / max(cutoff / 3.0, 1.0e-6)).astype(int)
        )
        self.cells: Dict[tuple[int, int, int], list[int]] = {}
        if len(self.coords):
            fractional = (self.coords @ self.inverse) % 1.0
            keys = np.floor(fractional * self.counts).astype(int) % self.counts
            for index, key in enumerate(keys):
                self.cells.setdefault(tuple(int(value) for value in key), []).append(index)

    def _minimum_vectors(self, point: np.ndarray, indices: np.ndarray) -> np.ndarray:
        delta_frac = (self.coords[indices] - point) @ self.inverse
        base = np.rint(delta_frac)
        trials = (
            delta_frac[:, None, :] - base[:, None, :] - _IMAGE_OFFSETS[None, :, :]
        ) @ self.vectors
        distances2 = np.einsum("ijk,ijk->ij", trials, trials)
        best = np.argmin(distances2, axis=1)
        return trials[np.arange(len(indices)), best]

    def candidate_indices(self, point: np.ndarray, radius: float) -> np.ndarray:
        fractional = (point @ self.inverse) % 1.0
        center = np.floor(fractional * self.counts).astype(int) % self.counts
        reach = np.maximum(
            1,
            np.ceil(radius * self.counts / self.spacings).astype(int),
        )
        keys = {
            tuple(int(value) for value in ((center + np.asarray(offset)) % self.counts))
            for offset in product(
                range(-int(reach[0]), int(reach[0]) + 1),
                range(-int(reach[1]), int(reach[1]) + 1),
                range(-int(reach[2]), int(reach[2]) + 1),
            )
        }
        indices = np.asarray(
            sorted({index for key in keys for index in self.cells.get(key, ())}),
            dtype=int,
        )
        return indices

    def query(self, point: np.ndarray, cutoff: float) -> tuple[np.ndarray, np.ndarray]:
        indices = self.candidate_indices(point, cutoff)
        if not len(indices):
            return indices, np.empty((0,), dtype=float)
        vectors = self._minimum_vectors(point, indices)
        distances = np.linalg.norm(vectors, axis=1)
        mask = distances <= cutoff
        return indices[mask], distances[mask]


def _switch_weights(distances: np.ndarray, switch: float, cutoff: float) -> np.ndarray:
    weights = np.ones_like(distances)
    weights[distances >= cutoff] = 0.0
    middle = (distances > switch) & (distances < cutoff)
    x = (distances[middle] - switch) / (cutoff - switch)
    weights[middle] = 1.0 - (3.0 * x * x - 2.0 * x * x * x)
    return weights


def _water_candidates(molecule: Molecule) -> list[_WaterCandidate]:
    result: list[_WaterCandidate] = []
    for chain_id, chain in molecule.chains.items():
        for residue in chain.residues:
            if residue.group not in WATER_MOL_TYPES:
                continue
            oxygen = residue.atoms.get("O") or residue.atoms.get("OW")
            if oxygen is None or oxygen.coord is None:
                continue
            result.append(
                _WaterCandidate(chain_id, residue, np.asarray(oxygen.coord, dtype=float))
            )
    return result


def _fixed_atoms(
    molecule: Molecule, parameters: SolvationVdwParameters
) -> tuple[list[_FixedAtom], list[_PlacedIon]]:
    fixed: list[_FixedAtom] = []
    retained: list[_PlacedIon] = []
    for chain in molecule.chains.values():
        for residue in chain.residues:
            if residue.group in WATER_MOL_TYPES:
                continue
            for atom in residue.atoms.values():
                if atom.coord is None or residue.group is None:
                    continue
                params = parameters.atom_params(residue.group, residue.ff_resname, atom.name)
                fixed.append(_FixedAtom(np.asarray(atom.coord), float(params.charge), float(params.sigma)))
                if residue.group in ION_MOL_TYPES:
                    retained.append(
                        _PlacedIon(
                            IonIdentity(residue.group, residue.ff_resname),
                            np.asarray(atom.coord), float(params.charge), float(params.sigma),
                            "retained_structural",
                        )
                    )
    for record in molecule.passthrough_atoms:
        if record.group in WATER_MOL_TYPES:
            continue
        if record.group is None or record.ff_resname is None:
            raise KeyError(f"Passthrough atom {record.resname}:{record.atom_name} lacks force-field identity")
        params = parameters.atom_params(record.group, record.ff_resname, record.atom_name)
        item = _FixedAtom(np.asarray(record.coord), float(params.charge), float(params.sigma))
        fixed.append(item)
        if record.group in ION_MOL_TYPES:
            retained.append(
                _PlacedIon(
                    IonIdentity(record.group, record.ff_resname), item.coord,
                    item.charge, item.sigma, "retained_structural",
                )
            )
    return fixed, retained


def _initial_potential(
    candidates: Sequence[_WaterCandidate],
    fixed: Sequence[_FixedAtom],
    vectors: np.ndarray,
    settings: IonPlacementSettings,
    candidate_index: Optional[_PeriodicCellList] = None,
) -> tuple[np.ndarray, _PeriodicCellList, int, np.ndarray]:
    """Calculate the local solute potential and identify its compact water shell.

    A water belongs to the neutralization shell when it lies within the
    electrostatic cutoff of at least one fixed atom.  The atom-centred cell
    queries avoid a solute-by-water Cartesian product and naturally restrict
    subsequent electrostatic selection to waters near the solute.
    """
    coords = [candidate.coord for candidate in candidates]
    if candidate_index is None:
        candidate_index = _PeriodicCellList(
            coords, vectors, settings.potential_cutoff_angstrom
        )
    potential = np.zeros(len(candidates), dtype=float)
    touched = np.zeros(len(candidates), dtype=bool)
    pair_count = 0
    for atom in fixed:
        indices, distances = candidate_index.query(
            atom.coord, settings.potential_cutoff_angstrom
        )
        nonzero = distances > 1.0e-9
        indices, distances = indices[nonzero], distances[nonzero]
        if len(indices):
            weights = _switch_weights(
                distances,
                settings.potential_switch_angstrom,
                settings.potential_cutoff_angstrom,
            )
            contributing = weights > 0.0
            indices, distances, weights = (
                indices[contributing], distances[contributing], weights[contributing]
            )
        if len(indices):
            potential[indices] += (
                atom.charge
                * weights
                / distances
            )
            touched[indices] = True
            pair_count += len(indices)
    return potential, candidate_index, pair_count, np.flatnonzero(touched)


def _update_potential_for_ion(
    potential: np.ndarray,
    candidate_index: _PeriodicCellList,
    ion: _PlacedIon,
    settings: IonPlacementSettings,
) -> np.ndarray:
    """Apply one ion's local potential and return the modified water indices."""
    indices, distances = candidate_index.query(
        ion.coord, settings.potential_cutoff_angstrom
    )
    nonzero = distances > 1.0e-9
    indices, distances = indices[nonzero], distances[nonzero]
    if not len(indices):
        return indices
    weights = _switch_weights(
        distances,
        settings.potential_switch_angstrom,
        settings.potential_cutoff_angstrom,
    )
    contributing = weights > 0.0
    indices, distances, weights = (
        indices[contributing], distances[contributing], weights[contributing]
    )
    if len(indices):
        potential[indices] += ion.charge * weights / distances
    return indices


def _added_ion_distance_diagnostics(
    ions: Sequence[_PlacedIon], vectors: np.ndarray
) -> tuple[Optional[float], int]:
    """Return the local minimum and number of suspiciously close added-ion pairs."""
    if len(ions) < 2:
        return None, 0
    diagnostic_cutoff = 6.0
    index = _PeriodicCellList(
        [ion.coord for ion in ions], vectors, diagnostic_cutoff
    )
    minimum = math.inf
    below_three = 0
    for left, ion in enumerate(ions):
        neighbours, distances = index.query(ion.coord, diagnostic_cutoff)
        for right, distance in zip(neighbours, distances):
            if int(right) <= left:
                continue
            value = float(distance)
            minimum = min(minimum, value)
            if value < 3.0:
                below_three += 1
    return (minimum if math.isfinite(minimum) else None), below_three


def _sterically_allowed(
    coord: np.ndarray,
    ion_sigma: float,
    fixed_index: _PeriodicCellList,
    fixed: Sequence[_FixedAtom],
    maximum_fixed_sigma: float,
) -> bool:
    maximum_sigma = max(maximum_fixed_sigma, ion_sigma)
    max_cutoff = 5.0 * (ion_sigma + maximum_sigma)
    indices, distances = fixed_index.query(coord, max_cutoff)
    for index, distance in zip(indices, distances):
        if distance < 5.0 * (ion_sigma + fixed[int(index)].sigma):
            return False
    return True


def _ion_position_is_separated(
    candidate_number: int,
    ion_sigma: float,
    candidate_index: _PeriodicCellList,
    occupied_candidates: Mapping[int, _PlacedIon],
    minimum_separation: float,
) -> bool:
    """Check only already occupied water sites, using the persistent water index."""
    if not occupied_candidates:
        return True
    coord = candidate_index.coords[candidate_number]
    neighbours, distances = candidate_index.query(coord, 12.0)
    for neighbour, distance in zip(neighbours, distances):
        placed = occupied_candidates.get(int(neighbour))
        if placed is None:
            continue
        required = max(
            minimum_separation,
            5.0 * (ion_sigma + placed.sigma),
        )
        if float(distance) < required:
            return False
    return True


class _MutableMinTree:
    """Array-backed minimum tree with vectorized batches of point updates."""

    def __init__(self, values: np.ndarray) -> None:
        self.count = len(values)
        self.size = 1 if self.count <= 1 else 1 << (self.count - 1).bit_length()
        self.tree = np.full(2 * self.size, np.inf, dtype=float)
        self.tree[self.size : self.size + self.count] = values
        for node in range(self.size - 1, 0, -1):
            self.tree[node] = min(self.tree[2 * node], self.tree[2 * node + 1])

    def argmin(self) -> Optional[int]:
        if self.count == 0 or not math.isfinite(float(self.tree[1])):
            return None
        node = 1
        while node < self.size:
            left = 2 * node
            node = left if self.tree[left] <= self.tree[left + 1] else left + 1
        index = node - self.size
        return index if index < self.count else None

    def set_many(self, indices: np.ndarray, values: np.ndarray) -> None:
        if not len(indices):
            return
        indices = np.asarray(indices, dtype=int)
        self.tree[self.size + indices] = np.asarray(values, dtype=float)
        nodes = np.unique((self.size + indices) // 2)
        while len(nodes):
            self.tree[nodes] = np.minimum(
                self.tree[2 * nodes], self.tree[2 * nodes + 1]
            )
            nodes = np.unique(nodes[nodes > 1] // 2)

    def disable(self, index: int) -> None:
        self.set_many(np.asarray([index]), np.asarray([np.inf]))


class _NeutralizationSelector:
    """Local-potential selector backed by a vectorized mutable minimum tree."""

    def __init__(
        self,
        identity: IonIdentity,
        candidates: Sequence[_WaterCandidate],
        eligible_indices: np.ndarray,
        potential: np.ndarray,
        candidate_index: _PeriodicCellList,
        fixed_index: _PeriodicCellList,
        fixed: Sequence[_FixedAtom],
        maximum_fixed_sigma: float,
        occupied_candidates: Dict[int, _PlacedIon],
        parameters: SolvationVdwParameters,
        minimum_separation: float,
    ) -> None:
        _name, self.charge, self.sigma = _ion_params(parameters, identity)
        self.identity = identity
        self.candidates = candidates
        self.potential = potential
        self.candidate_index = candidate_index
        self.fixed_index = fixed_index
        self.fixed = fixed
        self.maximum_fixed_sigma = maximum_fixed_sigma
        self.occupied_candidates = occupied_candidates
        self.minimum_separation = minimum_separation
        self.eligible_indices = np.asarray(eligible_indices, dtype=int)
        self.global_to_local = np.full(len(candidates), -1, dtype=int)
        self.global_to_local[self.eligible_indices] = np.arange(
            len(self.eligible_indices), dtype=int
        )
        scores = self.charge * potential[self.eligible_indices]
        inactive = np.asarray(
            [not candidates[int(index)].active for index in self.eligible_indices]
        )
        scores[inactive] = np.inf
        self.minimum_tree = _MutableMinTree(scores)

    def select(self) -> _PlacedIon:
        while True:
            local_index = self.minimum_tree.argmin()
            if local_index is None:
                break
            index = int(self.eligible_indices[local_index])
            candidate = self.candidates[index]
            if not candidate.active:
                self.minimum_tree.disable(local_index)
                continue
            if not (
                _sterically_allowed(
                    candidate.coord,
                    self.sigma,
                    self.fixed_index,
                    self.fixed,
                    self.maximum_fixed_sigma,
                )
                and _ion_position_is_separated(
                    index,
                    self.sigma,
                    self.candidate_index,
                    self.occupied_candidates,
                    self.minimum_separation,
                )
            ):
                self.minimum_tree.disable(local_index)
                continue
            candidate.active = False
            self.minimum_tree.disable(local_index)
            ion = _PlacedIon(
                self.identity,
                candidate.coord.copy(),
                self.charge,
                self.sigma,
                "charge_neutralization",
                index,
            )
            self.occupied_candidates[index] = ion
            return ion
        raise RuntimeError(
            f"No sterically valid neutralization-shell water remains for {self.identity}"
        )

    def update(self, ion: _PlacedIon, settings: IonPlacementSettings) -> int:
        affected = _update_potential_for_ion(
            self.potential, self.candidate_index, ion, settings
        )
        local_indices = self.global_to_local[affected]
        valid = local_indices >= 0
        global_indices = affected[valid]
        local_indices = local_indices[valid]
        if len(local_indices):
            active = np.asarray(
                [self.candidates[int(index)].active for index in global_indices]
            )
            global_indices = global_indices[active]
            local_indices = local_indices[active]
            self.minimum_tree.set_many(
                local_indices,
                self.charge * self.potential[global_indices],
            )
        return len(affected)


class _RandomSaltSelector:
    """Deterministic shuffled traversal of all remaining solvent positions."""

    def __init__(
        self,
        candidates: Sequence[_WaterCandidate],
        candidate_index: _PeriodicCellList,
        fixed_index: _PeriodicCellList,
        fixed: Sequence[_FixedAtom],
        maximum_fixed_sigma: float,
        occupied_candidates: Dict[int, _PlacedIon],
        parameters: SolvationVdwParameters,
        seed: int,
        minimum_separation: float,
    ) -> None:
        self.candidates = candidates
        self.candidate_index = candidate_index
        self.fixed_index = fixed_index
        self.fixed = fixed
        self.maximum_fixed_sigma = maximum_fixed_sigma
        self.occupied_candidates = occupied_candidates
        self.parameters = parameters
        self.minimum_separation = minimum_separation
        self.order = np.random.default_rng(seed).permutation(len(candidates))
        self.cursor = 0

    def select(self, identity: IonIdentity) -> _PlacedIon:
        _name, charge, sigma = _ion_params(self.parameters, identity)
        while self.cursor < len(self.order):
            index = int(self.order[self.cursor])
            self.cursor += 1
            candidate = self.candidates[index]
            if not candidate.active:
                continue
            if not (
                _sterically_allowed(
                    candidate.coord,
                    sigma,
                    self.fixed_index,
                    self.fixed,
                    self.maximum_fixed_sigma,
                )
                and _ion_position_is_separated(
                    index,
                    sigma,
                    self.candidate_index,
                    self.occupied_candidates,
                    self.minimum_separation,
                )
            ):
                continue
            candidate.active = False
            ion = _PlacedIon(
                identity,
                candidate.coord.copy(),
                charge,
                sigma,
                "salt_excess_random",
                index,
            )
            self.occupied_candidates[index] = ion
            return ion
        raise RuntimeError(
            f"No sterically valid randomly distributed water position remains for {identity}"
        )


def _formula_stoichiometry(cation_charge: int, anion_charge: int) -> tuple[int, int]:
    if cation_charge <= 0 or anion_charge >= 0:
        raise ValueError("Salt must contain a positive cation and a negative anion")
    divisor = math.gcd(abs(cation_charge), abs(anion_charge))
    return abs(anion_charge) // divisor, abs(cation_charge) // divisor


def _formula_units(
    concentration: float, water_count: int, box_vectors: np.ndarray, mode: str
) -> int:
    if mode == "water_ratio":
        expected = water_count * concentration / (1000.0 / 18.0)
    else:
        volume_litre = abs(float(np.linalg.det(box_vectors))) * 1.0e-27
        expected = concentration * volume_litre * 6.02214076e23
    return int(math.floor(expected + 0.5))


def _remove_selected_waters(molecule: Molecule, candidates: Sequence[_WaterCandidate]) -> None:
    removed = {id(candidate.residue) for candidate in candidates if not candidate.active}
    for chain_id in list(molecule.chains):
        chain = molecule.chains[chain_id]
        chain.residues = [residue for residue in chain.residues if id(residue) not in removed]
        for index, residue in enumerate(chain.residues):
            residue.index_in_chain = index
            residue.resseq = index + 1
        if not chain.residues:
            del molecule.chains[chain_id]


def _available_chain_ids(molecule: Molecule, preferred: str, count: int) -> Tuple[str, ...]:
    order = preferred + "JKLMNOPQRSTUVWXYZABCDEFGHabcdefghijklmnopqrstuvwxyz0123456789"
    unique = []
    for item in order:
        if item not in unique and item not in molecule.chains:
            unique.append(item)
    if len(unique) < count:
        raise ValueError("Not enough one-character PDB chain identifiers for ions")
    return tuple(unique[:count])


def _consolidate_ions(
    molecule: Molecule,
    new_ions: Sequence[_PlacedIon],
    parameters: SolvationVdwParameters,
) -> Tuple[str, ...]:
    all_ions: list[_PlacedIon] = list(new_ions)
    kept_passthrough: list[PDBAtomRecord] = []
    ion_record_ids: set[int] = set()
    for record in molecule.passthrough_atoms:
        if not _record_is_ion(record):
            kept_passthrough.append(record)
            continue
        identity = IonIdentity(str(record.group), str(record.ff_resname))
        _name, charge, sigma = _ion_params(parameters, identity)
        all_ions.append(
            _PlacedIon(identity, np.asarray(record.coord), charge, sigma, "retained_structural")
        )
        ion_record_ids.add(id(record))
    molecule.passthrough_atoms = kept_passthrough
    molecule.unassigned_records = [
        record for record in molecule.unassigned_records if id(record) not in ion_record_ids
    ]

    required_chains = max(1, int(math.ceil(len(all_ions) / 9999.0))) if all_ions else 0
    chain_ids = _available_chain_ids(molecule, "I", required_chains) if required_chains else ()
    for chunk_index, chain_id in enumerate(chain_ids):
        chunk = all_ions[chunk_index * 9999 : (chunk_index + 1) * 9999]
        residues: list[Residue] = []
        for index, ion in enumerate(chunk):
            atom_name, _charge, _sigma = _ion_params(parameters, ion.identity)
            atom = Atom(
                name=atom_name,
                element=_ion_element(ion.identity),
                coord=tuple(float(value) for value in ion.coord),
                built=ion.source != "retained_structural",
                build_source=ion.source,
                occupancy=1.0,
                bfactor=0.0,
            )
            residues.append(
                Residue(
                    chain_id=chain_id,
                    resseq=index + 1,
                    icode="",
                    ff_resname=ion.identity.resname,
                    atoms={atom_name: atom},
                    index_in_chain=index,
                    original_resname=ion.identity.resname,
                    group=ion.identity.mol_type,
                )
            )
        molecule.chains[chain_id] = Chain(chain_id, residues)

    # Preserve semantic ordering in Molecule as well as in the PDB writer.
    solute = {
        key: chain for key, chain in molecule.chains.items()
        if not chain.residues or chain.residues[0].group not in ION_MOL_TYPES | WATER_MOL_TYPES
    }
    ions = {
        key: chain for key, chain in molecule.chains.items()
        if chain.residues and chain.residues[0].group in ION_MOL_TYPES
    }
    waters = {
        key: chain for key, chain in molecule.chains.items()
        if chain.residues and chain.residues[0].group in WATER_MOL_TYPES
    }
    molecule.chains = {**solute, **ions, **waters}
    return chain_ids


def add_ions_to_solvated_molecule(
    molecule: Molecule,
    salts: Sequence[SaltSpecification],
    parameters: SolvationVdwParameters,
    settings: Optional[IonPlacementSettings] = None,
    *,
    modify_myself: bool = False,
) -> tuple[Molecule, IonAdditionReport]:
    """Neutralize locally, then distribute salt by reproducible water replacement."""
    started = perf_counter()
    settings = settings or IonPlacementSettings()
    settings.validate()
    target = molecule if modify_myself else copy.deepcopy(molecule)
    if target.periodic_box is None:
        raise ValueError("Ion placement requires a solvated molecule with periodic_box")
    vectors = np.asarray(target.periodic_box.vectors, dtype=float)
    candidates = _water_candidates(target)
    if not candidates:
        raise ValueError("Ion placement requires explicit W3/WAT waters")
    initial_water_count = len(candidates)
    fixed, retained_ions = _fixed_atoms(target, parameters)
    initial_charge_raw = sum(atom.charge for atom in fixed)
    initial_charge = _integer_charge(initial_charge_raw, "Fixed system")

    t0 = perf_counter()
    candidate_index = _PeriodicCellList(
        [candidate.coord for candidate in candidates],
        vectors,
        settings.potential_cutoff_angstrom,
    )
    maximum_fixed_sigma = max((atom.sigma for atom in fixed), default=0.0)
    fixed_index = _PeriodicCellList(
        [atom.coord for atom in fixed], vectors,
        max(5.0 * (maximum_fixed_sigma + 1.0), 1.0),
    )
    timings = {"spatial_indices": perf_counter() - t0}
    pair_count = 0
    neutralization_candidate_count = 0
    placed: list[_PlacedIon] = []
    occupied_candidates: Dict[int, _PlacedIon] = {}
    neutralization: Dict[str, int] = {}
    added_counts: Dict[str, int] = {}

    primary_cation = salts[0].cation if salts else settings.default_cation
    primary_anion = salts[0].anion if salts else settings.default_anion
    t0 = perf_counter()
    if initial_charge != 0:
        identity = primary_cation if initial_charge < 0 else primary_anion
        _name, ion_charge_raw, _sigma = _ion_params(parameters, identity)
        ion_charge = _integer_charge(ion_charge_raw, f"Neutralizing {identity}")
        if (-initial_charge) % ion_charge != 0:
            raise ValueError(
                f"System charge {initial_charge:+d} cannot be exactly neutralized by "
                f"{identity.resname} ({ion_charge:+d})"
            )
        count = (-initial_charge) // ion_charge
        potential, potential_index, initial_pairs, eligible_indices = _initial_potential(
            candidates, fixed, vectors, settings, candidate_index
        )
        # Reuse the equivalent water index made by _initial_potential.  Keeping
        # this construction local to charged systems avoids all electrostatic
        # setup for neutral solutes receiving salt excess only.
        candidate_index = potential_index
        pair_count += initial_pairs
        neutralization_candidate_count = len(eligible_indices)
        if neutralization_candidate_count < count:
            raise RuntimeError(
                "The local neutralization shell contains only "
                f"{neutralization_candidate_count} water positions for {count} ions"
            )
        selector = _NeutralizationSelector(
            identity,
            candidates,
            eligible_indices,
            potential,
            candidate_index,
            fixed_index,
            fixed,
            maximum_fixed_sigma,
            occupied_candidates,
            parameters,
            settings.minimum_ion_separation_angstrom,
        )
        for _ in range(count):
            ion = selector.select()
            placed.append(ion)
            pair_count += selector.update(ion, settings)
        neutralization[f"{identity.mol_type}:{identity.resname}"] = count
        added_counts[f"{identity.mol_type}:{identity.resname}"] = count
    timings["neutralization"] = perf_counter() - t0

    salt_reports: list[Dict[str, Any]] = []
    t0 = perf_counter()
    random_selector = _RandomSaltSelector(
        candidates,
        candidate_index,
        fixed_index,
        fixed,
        maximum_fixed_sigma,
        occupied_candidates,
        parameters,
        settings.salt_random_seed,
        settings.minimum_ion_separation_angstrom,
    )
    for salt in salts:
        _cn, cq_raw, _cs = _ion_params(parameters, salt.cation)
        _an, aq_raw, _as = _ion_params(parameters, salt.anion)
        cq = _integer_charge(cq_raw, f"{salt.cation}")
        aq = _integer_charge(aq_raw, f"{salt.anion}")
        n_cation, n_anion = _formula_stoichiometry(cq, aq)
        units = _formula_units(
            salt.concentration_mol_l, initial_water_count, vectors,
            settings.concentration_mode,
        )
        for _ in range(units):
            # Salt excess samples the complete remaining solvent rather than
            # the extrema of the solute potential.  The shared shuffled order
            # is deterministic, while the persistent local occupancy check
            # prevents accidental ion clusters without an all-ion scan.
            for identity, number in ((salt.cation, n_cation), (salt.anion, n_anion)):
                for _item in range(number):
                    ion = random_selector.select(identity)
                    placed.append(ion)
                    key = f"{identity.mol_type}:{identity.resname}"
                    added_counts[key] = added_counts.get(key, 0) + 1
        salt_reports.append(
            {
                "cation": f"{salt.cation.mol_type}:{salt.cation.resname}",
                "anion": f"{salt.anion.mol_type}:{salt.anion.resname}",
                "requested_concentration_mol_l": salt.concentration_mol_l,
                "formula_units": units,
                "cation_stoichiometry": n_cation,
                "anion_stoichiometry": n_anion,
                "actual_concentration_mol_l": (
                    units * (1000.0 / 18.0) / initial_water_count
                    if settings.concentration_mode == "water_ratio"
                    else units / (abs(float(np.linalg.det(vectors))) * 1.0e-27 * 6.02214076e23)
                ),
            }
        )
    timings["salt_excess"] = perf_counter() - t0

    t0 = perf_counter()
    minimum_ion_distance, close_ion_pairs = _added_ion_distance_diagnostics(
        placed, vectors
    )
    _remove_selected_waters(target, candidates)
    ion_chain_ids = _consolidate_ions(target, placed, parameters)
    timings["finalize"] = perf_counter() - t0
    timings["total"] = perf_counter() - started
    final_charge = initial_charge_raw + sum(ion.charge for ion in placed)
    return target, IonAdditionReport(
        initial_fixed_charge=initial_charge_raw,
        final_system_charge=final_charge,
        initial_waters=initial_water_count,
        final_waters=initial_water_count - len(placed),
        neutralization_ions=neutralization,
        salt_formula_units=salt_reports,
        added_ions=added_counts,
        ion_chain_ids=ion_chain_ids,
        potential_pair_evaluations=pair_count,
        neutralization_candidate_waters=neutralization_candidate_count,
        salt_random_seed=settings.salt_random_seed,
        minimum_added_ion_distance_angstrom=minimum_ion_distance,
        added_ion_pairs_below_3_angstrom=close_ion_pairs,
        timings_seconds=timings,
    )


def ion_report_text(
    cleanup: CrystalIonCleanupReport, addition: IonAdditionReport
) -> str:
    lines = [
        "FORGE ion report", "================",
        f"Input crystal ions: {cleanup.input_ions}",
        f"Removed monovalent crystal ions: {cleanup.removed_monovalent}",
        f"Removed nonstructural multivalent ions: {cleanup.removed_nonstructural_multivalent}",
        f"Retained structural monovalent ions: {cleanup.retained_structural_monovalent}",
        f"Retained structural multivalent ions: {cleanup.retained_structural_multivalent}",
        f"Structural ions replaced by Mg2+: {cleanup.replaced_by_magnesium}",
        f"Initial fixed charge: {addition.initial_fixed_charge:.6f} e",
        f"Final system charge: {addition.final_system_charge:.6f} e",
        f"Waters before/after ion replacement: {addition.initial_waters}/{addition.final_waters}",
        f"Ion chain(s): {', '.join(addition.ion_chain_ids)}",
        f"Potential pair evaluations: {addition.potential_pair_evaluations}",
        f"Neutralization-shell candidate waters: {addition.neutralization_candidate_waters}",
        f"Salt random seed: {addition.salt_random_seed}",
        "Minimum added-ion distance (A): "
        + (
            f"{addition.minimum_added_ion_distance_angstrom:.6f}"
            if addition.minimum_added_ion_distance_angstrom is not None
            else "n/a"
        ),
        f"Added-ion pairs below 3 A: {addition.added_ion_pairs_below_3_angstrom}",
        f"Neutralization ions: {addition.neutralization_ions}",
        f"All added ions: {addition.added_ions}",
        f"Timings (s): {addition.timings_seconds}",
        "Salt formula units:",
    ]
    lines.extend(f"  {entry}" for entry in addition.salt_formula_units)
    return "\n".join(lines) + "\n"
