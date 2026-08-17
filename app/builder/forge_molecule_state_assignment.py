#!/usr/bin/env python3
"""Protonation, tautomer, and covalent-state assignment for FORGE molecules.

The layer handles geometry-based covalent state pairs (currently CYS/CYX
disulfides) and protonation/tautomer assignment from hydrogen-bond geometry.
It operates after parsing and before build planning so that it can use input
coordinates while presenting a normalized residue inventory to the planner.
"""

from __future__ import annotations

import copy
import itertools
import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from forge_molecule_parser import Atom, AtomKey, Molecule, Residue, infer_element
from forge_molecule_builder import (
    AtomID,
    collect_available_atoms,
    compute_dihedral_deg,
    find_matching_torsion_group,
    first_available_member_index,
    member_dihedral_atoms,
    normalize_angle,
    place_atom_from_internal,
    resolve_refs,
    torsion_class,
)


class CovalentStateAssignmentError(ValueError):
    """Raised when geometry does not define an unambiguous covalent state."""


@dataclass(frozen=True)
class CovalentBondAssignment:
    atom1: AtomKey
    atom2: AtomKey
    distance_angstrom: float


@dataclass
class CovalentStateAssignmentReport:
    bonds: List[CovalentBondAssignment] = field(default_factory=list)
    state_changes: List[Tuple[Tuple[str, int, str], str, str]] = field(
        default_factory=list
    )
    missing_bond_atoms: List[Tuple[Tuple[str, int, str], str, str]] = field(
        default_factory=list
    )


@dataclass(frozen=True)
class ResidueIdentityAssignment:
    """Description of one generic residue-identity/inventory normalization."""

    residue_key: Tuple[str, int, str]
    old_mol_type: Optional[str]
    old_resname: str
    new_mol_type: str
    new_resname: str
    restored_observed_atoms: Tuple[str, ...]
    created_missing_atoms: Tuple[str, ...]
    moved_to_observed_extras: Tuple[str, ...]

    @property
    def identity_changed(self) -> bool:
        return (
            self.old_mol_type != self.new_mol_type
            or self.old_resname != self.new_resname
        )


@dataclass(frozen=True)
class HydrogenBondGeometrySettings:
    heavy_atom_cutoff_angstrom: float = 3.3
    max_donor_deviation_deg: float = 40.0


@dataclass(frozen=True)
class FixedSiteEvidence:
    variable_site: AtomKey
    required_role: str
    partner_site: AtomKey
    distance_angstrom: float
    donor_deviation_deg: float


@dataclass(frozen=True)
class VariableSiteContact:
    site1: AtomKey
    site2: AtomKey
    distance_angstrom: float
    allowed_role_pairs: Tuple[Tuple[str, str], ...]
    donor_deviation_site1_to_site2: Optional[float] = None
    donor_deviation_site2_to_site1: Optional[float] = None


@dataclass(frozen=True)
class AmbivalentSiteContact:
    variable_site: AtomKey
    partner_site: AtomKey
    distance_angstrom: float


@dataclass(frozen=True)
class ProtonationConflict:
    kind: str
    message: str
    sites: Tuple[AtomKey, ...] = ()


@dataclass(frozen=True)
class ProtonationResidueAssignment:
    residue_key: Tuple[str, int, str]
    old_mol_type: Optional[str]
    old_resname: str
    new_mol_type: str
    new_resname: str
    default_resname: str
    is_default: bool


@dataclass
class ProtonationAssignmentReport:
    pH: float
    family_defaults: List[Tuple[Tuple[str, ...], str]] = field(default_factory=list)
    fixed_evidence: List[FixedSiteEvidence] = field(default_factory=list)
    variable_contacts: List[VariableSiteContact] = field(default_factory=list)
    ambivalent_contacts: List[AmbivalentSiteContact] = field(default_factory=list)
    site_statuses: Dict[AtomKey, str] = field(default_factory=dict)
    assignments: List[ProtonationResidueAssignment] = field(default_factory=list)
    conflicts: List[ProtonationConflict] = field(default_factory=list)
    unevaluable_sites: List[Tuple[AtomKey, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def has_conflicts(self) -> bool:
        return bool(self.conflicts)


@dataclass(frozen=True)
class StateAssignmentReport:
    covalent: CovalentStateAssignmentReport
    protonation: ProtonationAssignmentReport


@dataclass(frozen=True)
class _StateDefinition:
    family: str
    unbound_mol_type: str
    unbound_resname: str
    bound_mol_type: str
    bound_resname: str
    bond_atom: str
    partners: frozenset[str]


@dataclass(frozen=True)
class _Candidate:
    residue: Residue
    definition: _StateDefinition
    atom_key: AtomKey
    coord: Tuple[float, float, float]


@dataclass(frozen=True)
class _ProtonationStateOption:
    mol_type: str
    resname: str
    level_index: int
    tautomer_index: int
    site_roles: Mapping[str, FrozenSet[str]]
    site_donor_hydrogens: Mapping[str, Tuple[str, ...]]


@dataclass(frozen=True)
class _ProtonationFamily:
    index: int
    options: Tuple[_ProtonationStateOption, ...]
    levels: Tuple[Tuple[_ProtonationStateOption, ...], ...]
    pka_values: Tuple[float, ...]
    variable_sites: Tuple[str, ...]
    default_level_index: int
    default_option: _ProtonationStateOption


@dataclass(frozen=True)
class _TitratableResidue:
    index: int
    residue: Residue
    family: _ProtonationFamily


@dataclass(frozen=True)
class _HBondSite:
    residue: Residue
    atom_name: str
    coord: Tuple[float, float, float]
    acceptor: bool
    donor_hydrogens: Tuple[str, ...]
    variable_owner: Optional[int]

    @property
    def atom_key(self) -> AtomKey:
        return _atom_key(self.residue, self.atom_name)


def _residue_key(residue: Residue) -> Tuple[str, int, str]:
    return residue.chain_id, residue.resseq, residue.icode


def _atom_key(residue: Residue, atom_name: str) -> AtomKey:
    return residue.chain_id, residue.resseq, residue.icode, atom_name


def _canonical_bond(atom1: AtomKey, atom2: AtomKey) -> Tuple[AtomKey, AtomKey]:
    return min((atom1, atom2), (atom2, atom1))


def _load_state_definitions(
    state_data: Mapping[str, Any],
) -> Tuple[List[_StateDefinition], Dict[Tuple[str, str], _StateDefinition]]:
    definitions: List[_StateDefinition] = []
    by_state: Dict[Tuple[str, str], _StateDefinition] = {}
    for family, raw in state_data.get("covalent_state_pairs", {}).items():
        unbound = raw["unbound"]
        bound = raw["bound"]
        definition = _StateDefinition(
            family=family,
            unbound_mol_type=unbound["mol_type"],
            unbound_resname=unbound["resn"],
            bound_mol_type=bound["mol_type"],
            bound_resname=bound["resn"],
            bond_atom=raw["bond_atom"],
            partners=frozenset(raw.get("partners", [])),
        )
        definitions.append(definition)
        for state_key in (
            (definition.unbound_mol_type, definition.unbound_resname),
            (definition.bound_mol_type, definition.bound_resname),
        ):
            if state_key in by_state:
                raise ValueError(
                    f"Covalent state {state_key!r} belongs to more than one family"
                )
            by_state[state_key] = definition
    return definitions, by_state


def _iter_residues(molecule: Molecule):
    for chain in molecule.chains.values():
        yield from chain.residues


def _compatible(candidate1: _Candidate, candidate2: _Candidate) -> bool:
    definition1 = candidate1.definition
    definition2 = candidate2.definition
    return (
        definition2.family in definition1.partners
        and definition1.family in definition2.partners
    )


def _find_unambiguous_pairs(
    candidates: Sequence[_Candidate], cutoff_angstrom: float
) -> List[Tuple[_Candidate, _Candidate, float]]:
    edges: List[Tuple[int, int, float]] = []
    neighbors: Dict[int, List[Tuple[int, float]]] = {
        index: [] for index in range(len(candidates))
    }
    for left in range(len(candidates)):
        for right in range(left + 1, len(candidates)):
            candidate1 = candidates[left]
            candidate2 = candidates[right]
            if not _compatible(candidate1, candidate2):
                continue
            distance = math.dist(candidate1.coord, candidate2.coord)
            if distance <= cutoff_angstrom:
                edges.append((left, right, distance))
                neighbors[left].append((right, distance))
                neighbors[right].append((left, distance))

    ambiguous = [index for index, adjacent in neighbors.items() if len(adjacent) > 1]
    if ambiguous:
        details: List[str] = []
        for index in ambiguous:
            candidate = candidates[index]
            partners = ", ".join(
                f"{candidates[other].atom_key[0]}:"
                f"{candidates[other].atom_key[1]}{candidates[other].atom_key[2]} "
                f"({distance:.3f} A)"
                for other, distance in sorted(neighbors[index], key=lambda item: item[1])
            )
            details.append(
                f"{candidate.atom_key[0]}:{candidate.atom_key[1]}"
                f"{candidate.atom_key[2]} {candidate.atom_key[3]} -> {partners}"
            )
        raise CovalentStateAssignmentError(
            "Ambiguous covalent-bond cluster within "
            f"{cutoff_angstrom:.3f} A cutoff: " + "; ".join(details)
        )

    return [
        (candidates[left], candidates[right], distance)
        for left, right, distance in edges
    ]


def _new_missing_atom(atom_name: str) -> Atom:
    return Atom(
        name=atom_name,
        element=infer_element(atom_name),
        coord=None,
        serial=None,
        built=False,
        build_source=None,
    )


def assign_residue_identity(
    residue: Residue,
    mol_type: str,
    resname: str,
    converting: Mapping[str, Any],
    *,
    assignment_source: str = "state_assignment",
) -> ResidueIdentityAssignment:
    """Assign a converting-dictionary residue identity and reconcile its atoms.

    This operation is intentionally independent of the algorithm that selected
    the state.  Covalent, protonation, and tautomer assignment can therefore use
    the same lossless inventory transition:

    * expected atoms already active on the residue are retained;
    * expected atoms preserved in ``observed_extra_atoms`` are restored;
    * other expected atoms are created with ``coord=None``;
    * no-longer-expected atoms with coordinates are preserved as observed extras.
    """
    if mol_type not in converting or resname not in converting[mol_type]:
        raise KeyError(
            f"Residue identity {mol_type}/{resname} is absent from converting dictionary"
        )

    old_mol_type = residue.group
    old_resname = residue.ff_resname
    expected_names = list(converting[mol_type][resname]["atom"].keys())
    expected_set = set(expected_names)

    previous_active = dict(residue.atoms)
    previous_extras = dict(residue.observed_extra_atoms)
    moved_to_observed_extras = sorted(
        atom_name
        for atom_name, atom in previous_active.items()
        if atom_name not in expected_set and atom.coord is not None
    )
    available: Dict[str, Atom] = dict(previous_extras)
    available.update(previous_active)
    new_atoms: Dict[str, Atom] = {}
    restored_observed_atoms: List[str] = []
    created_missing_atoms: List[str] = []
    for atom_name in expected_names:
        atom = available.pop(atom_name, None)
        if atom is None:
            atom = _new_missing_atom(atom_name)
            created_missing_atoms.append(atom_name)
        elif atom_name in previous_extras and atom_name not in previous_active:
            restored_observed_atoms.append(atom_name)
        new_atoms[atom_name] = atom

    residue.observed_extra_atoms = {
        atom_name: atom
        for atom_name, atom in available.items()
        if atom.coord is not None
    }
    residue.atoms = new_atoms
    residue.group = mol_type
    residue.ff_resname = resname
    residue.state_assignment_source = assignment_source
    residue.connectivity_parts = [
        [atom_name for atom_name in part if atom_name in expected_set]
        for part in residue.connectivity_parts
    ]
    residue.connectivity_parts = [
        part for part in residue.connectivity_parts if part
    ]

    return ResidueIdentityAssignment(
        residue_key=_residue_key(residue),
        old_mol_type=old_mol_type,
        old_resname=old_resname,
        new_mol_type=mol_type,
        new_resname=resname,
        restored_observed_atoms=tuple(restored_observed_atoms),
        created_missing_atoms=tuple(created_missing_atoms),
        moved_to_observed_extras=tuple(moved_to_observed_extras),
    )


def assign_covalent_states(
    molecule: Molecule,
    converting: Mapping[str, Any],
    state_data: Mapping[str, Any],
    *,
    cutoff_angstrom: float = 2.3,
    modify_myself: bool = False,
) -> Tuple[Molecule, CovalentStateAssignmentReport]:
    """Normalize configured covalent residue states from bond-atom geometry.

    Every configured residue is evaluated regardless of its input bound/unbound
    label.  A residue without coordinates for its configured ``bond_atom`` is
    assigned to the unbound state.  Any atom with more than one compatible
    partner inside the cutoff makes the geometry ambiguous and raises
    :class:`CovalentStateAssignmentError`.
    """
    if not math.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0.0:
        raise ValueError("cutoff_angstrom must be a positive finite number")

    target = molecule if modify_myself else copy.deepcopy(molecule)
    report = CovalentStateAssignmentReport()
    _, by_state = _load_state_definitions(state_data)

    configured: List[Tuple[Residue, _StateDefinition]] = []
    candidates: List[_Candidate] = []
    for residue in _iter_residues(target):
        definition = by_state.get((residue.group or "", residue.ff_resname))
        if definition is None:
            continue
        configured.append((residue, definition))
        atom = residue.atoms.get(definition.bond_atom)
        if atom is None:
            atom = residue.observed_extra_atoms.get(definition.bond_atom)
        if atom is None or atom.coord is None:
            report.missing_bond_atoms.append(
                (_residue_key(residue), residue.ff_resname, definition.bond_atom)
            )
            continue
        candidates.append(
            _Candidate(
                residue=residue,
                definition=definition,
                atom_key=_atom_key(residue, definition.bond_atom),
                coord=atom.coord,
            )
        )

    pairs = _find_unambiguous_pairs(candidates, cutoff_angstrom)
    bound_residue_ids = {
        id(candidate.residue)
        for pair in pairs
        for candidate in pair[:2]
    }

    for residue, definition in configured:
        if id(residue) in bound_residue_ids:
            identity_assignment = assign_residue_identity(
                residue,
                definition.bound_mol_type,
                definition.bound_resname,
                converting,
                assignment_source="covalent_state_geometry",
            )
        else:
            identity_assignment = assign_residue_identity(
                residue,
                definition.unbound_mol_type,
                definition.unbound_resname,
                converting,
                assignment_source="covalent_state_geometry",
            )
        if identity_assignment.identity_changed:
            report.state_changes.append(
                (
                    identity_assignment.residue_key,
                    identity_assignment.old_resname,
                    identity_assignment.new_resname,
                )
            )

    assigned_bonds: List[Tuple[AtomKey, AtomKey]] = []
    for candidate1, candidate2, distance in pairs:
        atom1, atom2 = _canonical_bond(candidate1.atom_key, candidate2.atom_key)
        assigned_bonds.append((atom1, atom2))
        report.bonds.append(
            CovalentBondAssignment(
                atom1=atom1,
                atom2=atom2,
                distance_angstrom=distance,
            )
        )
    managed_atom_keys = {
        _atom_key(residue, definition.bond_atom)
        for residue, definition in configured
    }
    preserved_bonds = [
        bond
        for bond in target.covalent_bonds
        if bond[0] not in managed_atom_keys and bond[1] not in managed_atom_keys
    ]
    target.covalent_bonds = sorted(set(preserved_bonds + assigned_bonds))

    return target, report


# -----------------------------------------------------------------------------
# Protonation and tautomer assignment
# -----------------------------------------------------------------------------

def _roles_from_site(raw: Mapping[str, Any]) -> FrozenSet[str]:
    roles: Set[str] = set()
    if bool(raw.get("acceptor", False)):
        roles.add("acceptor")
    if raw.get("donor_hydrogens", []):
        roles.add("donor")
    return frozenset(roles)


def _default_level_for_ph(pH: float, pka_values: Sequence[float]) -> int:
    level = 0
    for pka in pka_values:
        if pH > pka:
            level += 1
        elif math.isclose(pH, pka, rel_tol=0.0, abs_tol=1.0e-12):
            # At equality choose the side approached from physiological pH 7.
            # For an exact pKa of 7, use the less protonated side deterministically.
            if pka <= 7.0:
                level += 1
        else:
            break
    return level


def _load_protonation_families(
    state_data: Mapping[str, Any],
    pH: float,
) -> Tuple[List[_ProtonationFamily], Dict[Tuple[str, str], _ProtonationFamily]]:
    families: List[_ProtonationFamily] = []
    by_state: Dict[Tuple[str, str], _ProtonationFamily] = {}
    for family_index, raw_family in enumerate(
        state_data.get("protonation_families", [])
    ):
        raw_levels = raw_family.get("members", [])
        pka_values = tuple(float(value) for value in raw_family.get("pKa", []))
        if len(pka_values) != max(0, len(raw_levels) - 1):
            raise ValueError(
                f"Protonation family {family_index} must have one fewer pKa values "
                "than protonation levels"
            )
        variable_raw = raw_family.get("variable_sites", {})
        levels: List[Tuple[_ProtonationStateOption, ...]] = []
        options: List[_ProtonationStateOption] = []
        for level_index, raw_level in enumerate(raw_levels):
            level_options: List[_ProtonationStateOption] = []
            for tautomer_index, raw_state in enumerate(raw_level):
                mol_type = str(raw_state["mol_type"])
                resname = str(raw_state["resn"])
                site_roles: Dict[str, FrozenSet[str]] = {}
                donor_hydrogens: Dict[str, Tuple[str, ...]] = {}
                for site_name, states in variable_raw.items():
                    if resname not in states:
                        raise ValueError(
                            f"Variable site {site_name} lacks state {resname} "
                            f"in protonation family {family_index}"
                        )
                    site_data = states[resname]
                    site_roles[site_name] = _roles_from_site(site_data)
                    donor_hydrogens[site_name] = tuple(
                        str(name) for name in site_data.get("donor_hydrogens", [])
                    )
                option = _ProtonationStateOption(
                    mol_type=mol_type,
                    resname=resname,
                    level_index=level_index,
                    tautomer_index=tautomer_index,
                    site_roles=site_roles,
                    site_donor_hydrogens=donor_hydrogens,
                )
                options.append(option)
                level_options.append(option)
            if not level_options:
                raise ValueError(
                    f"Protonation family {family_index} contains an empty level"
                )
            levels.append(tuple(level_options))
        if not levels:
            continue
        default_level_index = _default_level_for_ph(pH, pka_values)
        if default_level_index >= len(levels):
            default_level_index = len(levels) - 1
        family = _ProtonationFamily(
            index=family_index,
            options=tuple(options),
            levels=tuple(levels),
            pka_values=pka_values,
            variable_sites=tuple(str(name) for name in variable_raw),
            default_level_index=default_level_index,
            default_option=levels[default_level_index][0],
        )
        families.append(family)
        for option in options:
            key = (option.mol_type, option.resname)
            if key in by_state:
                raise ValueError(
                    f"Protonation state {key!r} belongs to more than one family"
                )
            by_state[key] = family
    return families, by_state


def _coord_for_atom_key(molecule: Molecule, atom_key: AtomKey):
    chain_id, resseq, icode, atom_name = atom_key
    chain = molecule.chains.get(chain_id)
    if chain is None:
        return None
    for residue in chain.residues:
        if residue.resseq == resseq and residue.icode == icode:
            atom = residue.atoms.get(atom_name)
            if atom is None:
                atom = residue.observed_extra_atoms.get(atom_name)
            return None if atom is None else atom.coord
    return None


def _resolved_coord(molecule: Molecule, resolved):
    residue = molecule.chains[resolved.chain_id].residues[resolved.residue_index]
    atom = residue.atoms[resolved.atom_name]
    return atom.coord


def _internal_values(rule: Mapping[str, Any]) -> Tuple[float, float, float]:
    internal = rule.get("internal", {})
    return (
        float(internal["r"]),
        float(internal["angle"]),
        float(internal["dihedral"]),
    )


def _virtual_hydrogen_from_rule(
    molecule: Molecule,
    residue: Residue,
    residue_template: Mapping[str, Any],
    hydrogen_name: str,
    rule: Mapping[str, Any],
    acceptor_coord: Tuple[float, float, float],
) -> Optional[Tuple[float, float, float]]:
    refs = resolve_refs(molecule, residue, rule.get("refs", []))
    if refs is None or len(refs) != 3:
        return None
    coords = [_resolved_coord(molecule, ref) for ref in refs]
    if any(coord is None for coord in coords):
        return None
    r, angle, template_dihedral = _internal_values(rule)
    cls = torsion_class(rule)
    if cls == "rigid":
        dihedral = template_dihedral
    elif cls == "free_rotor_hydrogen":
        try:
            dihedral = compute_dihedral_deg(
                coords[0], coords[1], coords[2], acceptor_coord
            )
        except RuntimeError:
            return None
    elif cls == "derived_rotatable":
        target_id = AtomID(residue.chain_id, residue.index_in_chain, hydrogen_name)
        previous_atom = residue.atoms.get(hydrogen_name)
        added_target = previous_atom is None
        if added_target:
            residue.atoms[hydrogen_name] = _new_missing_atom(hydrogen_name)
        try:
            available = collect_available_atoms(molecule)
            group_index, target_member_index, group, _reversed = (
                find_matching_torsion_group(
                    molecule,
                    residue,
                    residue_template,
                    rule,
                    target_id,
                )
            )
            del group_index
            observed_index = first_available_member_index(
                molecule, residue, group, available
            )
            if observed_index is None:
                return None
            observed_atoms = member_dihedral_atoms(
                molecule, residue, group, observed_index
            )
            if observed_atoms is None:
                return None
            observed_coords = [_resolved_coord(molecule, ref) for ref in observed_atoms]
            if any(coord is None for coord in observed_coords):
                return None
            observed_phi = compute_dihedral_deg(*observed_coords)
            members = group.get("members", [])
            observed_phase = float(members[observed_index].get("phase", 0.0))
            target_phase = float(members[target_member_index].get("phase", 0.0))
            dihedral = normalize_angle(
                observed_phi + target_phase - observed_phase
            )
        except (RuntimeError, KeyError, IndexError, ValueError):
            return None
        finally:
            if added_target:
                residue.atoms.pop(hydrogen_name, None)
    else:
        return None
    try:
        return place_atom_from_internal(
            coords[0], coords[1], coords[2], r, angle, dihedral
        )
    except RuntimeError:
        return None


def _virtual_hydrogen_coord(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    residue: Residue,
    mol_type: str,
    resname: str,
    hydrogen_name: str,
    acceptor_coord: Tuple[float, float, float],
) -> Optional[Tuple[float, float, float]]:
    existing = residue.atoms.get(hydrogen_name)
    if existing is None:
        existing = residue.observed_extra_atoms.get(hydrogen_name)
    if existing is not None and existing.coord is not None:
        return existing.coord
    try:
        rt = building_template[mol_type][resname]
    except KeyError:
        return None
    for rule in rt.get("build_rules", {}).get(hydrogen_name, []):
        coord = _virtual_hydrogen_from_rule(
            molecule,
            residue,
            rt,
            hydrogen_name,
            rule,
            acceptor_coord,
        )
        if coord is not None:
            return coord
    return None


def _donor_deviation_deg(
    hydrogen_coord: Tuple[float, float, float],
    donor_coord: Tuple[float, float, float],
    acceptor_coord: Tuple[float, float, float],
) -> Optional[float]:
    hx = tuple(h - x for h, x in zip(hydrogen_coord, donor_coord))
    yx = tuple(y - x for y, x in zip(acceptor_coord, donor_coord))
    nh = math.sqrt(sum(value * value for value in hx))
    ny = math.sqrt(sum(value * value for value in yx))
    if nh < 1.0e-12 or ny < 1.0e-12:
        return None
    cosine = sum(a * b for a, b in zip(hx, yx)) / (nh * ny)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def _best_donor_angle_for_options(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    titratable: _TitratableResidue,
    site_name: str,
    acceptor_coord: Tuple[float, float, float],
) -> Optional[float]:
    donor_atom = titratable.residue.atoms.get(site_name)
    if donor_atom is None or donor_atom.coord is None:
        return None
    angles: List[float] = []
    seen: Set[Tuple[str, str, str]] = set()
    for option in titratable.family.options:
        if "donor" not in option.site_roles.get(site_name, frozenset()):
            continue
        for hydrogen_name in option.site_donor_hydrogens.get(site_name, ()):
            key = (option.mol_type, option.resname, hydrogen_name)
            if key in seen:
                continue
            seen.add(key)
            coord = _virtual_hydrogen_coord(
                molecule,
                building_template,
                titratable.residue,
                option.mol_type,
                option.resname,
                hydrogen_name,
                acceptor_coord,
            )
            if coord is None:
                continue
            angle = _donor_deviation_deg(
                coord, donor_atom.coord, acceptor_coord
            )
            if angle is not None:
                angles.append(angle)
    return min(angles) if angles else None


def _best_fixed_donor_angle(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    site: _HBondSite,
    acceptor_coord: Tuple[float, float, float],
) -> Optional[float]:
    angles: List[float] = []
    for hydrogen_name in site.donor_hydrogens:
        coord = _virtual_hydrogen_coord(
            molecule,
            building_template,
            site.residue,
            site.residue.group or "",
            site.residue.ff_resname,
            hydrogen_name,
            acceptor_coord,
        )
        if coord is None:
            continue
        angle = _donor_deviation_deg(coord, site.coord, acceptor_coord)
        if angle is not None:
            angles.append(angle)
    return min(angles) if angles else None


def _collect_hbond_sites(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    titratables: Sequence[_TitratableResidue],
) -> Tuple[List[_HBondSite], Dict[AtomKey, _HBondSite]]:
    owner_by_residue = {id(item.residue): item for item in titratables}
    sites: List[_HBondSite] = []
    by_key: Dict[AtomKey, _HBondSite] = {}
    for residue in _iter_residues(molecule):
        if residue.group not in {"R", "D", "P"}:
            continue
        try:
            rt = building_template[residue.group][residue.ff_resname]
        except KeyError:
            continue
        owner = owner_by_residue.get(id(residue))
        variable_names = set(owner.family.variable_sites) if owner else set()
        for atom_name, site_data in rt.get("hbond_sites", {}).items():
            atom = residue.atoms.get(atom_name)
            if atom is None or atom.coord is None:
                continue
            site = _HBondSite(
                residue=residue,
                atom_name=atom_name,
                coord=atom.coord,
                acceptor=bool(site_data.get("acceptor", False)),
                donor_hydrogens=tuple(site_data.get("donor_hydrogens", [])),
                variable_owner=(owner.index if atom_name in variable_names else None),
            )
            sites.append(site)
            by_key[site.atom_key] = site
    return sites, by_key


def _cell_key(coord: Tuple[float, float, float], size: float) -> Tuple[int, int, int]:
    return tuple(math.floor(value / size) for value in coord)  # type: ignore[return-value]


def _nearby_sites(
    variable_sites: Sequence[_HBondSite],
    all_sites: Sequence[_HBondSite],
    cutoff: float,
) -> Iterable[Tuple[_HBondSite, _HBondSite, float]]:
    cells: Dict[Tuple[int, int, int], List[_HBondSite]] = defaultdict(list)
    for site in all_sites:
        cells[_cell_key(site.coord, cutoff)].append(site)
    seen_variable_pairs: Set[Tuple[AtomKey, AtomKey]] = set()
    for variable in variable_sites:
        center = _cell_key(variable.coord, cutoff)
        for shift in itertools.product((-1, 0, 1), repeat=3):
            neighbor_cell = tuple(center[i] + shift[i] for i in range(3))
            for partner in cells.get(neighbor_cell, []):
                if partner.atom_key == variable.atom_key:
                    continue
                if partner.residue is variable.residue:
                    continue
                distance = math.dist(variable.coord, partner.coord)
                if distance > cutoff:
                    continue
                if partner.variable_owner is not None:
                    pair_key = _canonical_bond(variable.atom_key, partner.atom_key)
                    if pair_key in seen_variable_pairs:
                        continue
                    seen_variable_pairs.add(pair_key)
                yield variable, partner, distance


def _option_roles(option: _ProtonationStateOption, site_name: str) -> FrozenSet[str]:
    return option.site_roles.get(site_name, frozenset())


def _analyze_hbond_constraints(
    molecule: Molecule,
    building_template: Mapping[str, Any],
    titratables: Sequence[_TitratableResidue],
    settings: HydrogenBondGeometrySettings,
    report: ProtonationAssignmentReport,
) -> Tuple[Dict[AtomKey, FixedSiteEvidence], List[VariableSiteContact]]:
    all_sites, _site_by_key = _collect_hbond_sites(
        molecule, building_template, titratables
    )
    titratable_by_index = {item.index: item for item in titratables}
    variable_sites = [site for site in all_sites if site.variable_owner is not None]
    found_variable_keys = {site.atom_key for site in variable_sites}
    for item in titratables:
        for site_name in item.family.variable_sites:
            key = _atom_key(item.residue, site_name)
            if key not in found_variable_keys:
                report.unevaluable_sites.append((key, "missing variable-site coordinates"))

    evidence_by_site_role: Dict[Tuple[AtomKey, str], List[FixedSiteEvidence]] = defaultdict(list)
    variable_contacts: List[VariableSiteContact] = []
    ambivalent_keys: Set[Tuple[AtomKey, AtomKey]] = set()

    for variable, partner, distance in _nearby_sites(
        variable_sites,
        all_sites,
        settings.heavy_atom_cutoff_angstrom,
    ):
        owner = titratable_by_index[variable.variable_owner]  # type: ignore[index]
        if partner.variable_owner is not None:
            other = titratable_by_index[partner.variable_owner]
            allowed: List[Tuple[str, str]] = []
            angle_forward = _best_donor_angle_for_options(
                molecule, building_template, owner, variable.atom_name, partner.coord
            )
            if (
                angle_forward is not None
                and angle_forward <= settings.max_donor_deviation_deg
                and any(
                    "acceptor" in _option_roles(option, partner.atom_name)
                    for option in other.family.options
                )
            ):
                allowed.append(("donor", "acceptor"))
            angle_reverse = _best_donor_angle_for_options(
                molecule, building_template, other, partner.atom_name, variable.coord
            )
            if (
                angle_reverse is not None
                and angle_reverse <= settings.max_donor_deviation_deg
                and any(
                    "acceptor" in _option_roles(option, variable.atom_name)
                    for option in owner.family.options
                )
            ):
                allowed.append(("acceptor", "donor"))
            if allowed:
                contact = VariableSiteContact(
                    site1=variable.atom_key,
                    site2=partner.atom_key,
                    distance_angstrom=distance,
                    allowed_role_pairs=tuple(allowed),
                    donor_deviation_site1_to_site2=angle_forward,
                    donor_deviation_site2_to_site1=angle_reverse,
                )
                variable_contacts.append(contact)
                report.variable_contacts.append(contact)
            continue

        partner_is_donor = bool(partner.donor_hydrogens)
        partner_is_acceptor = partner.acceptor
        if partner_is_acceptor and not partner_is_donor:
            angle = _best_donor_angle_for_options(
                molecule, building_template, owner, variable.atom_name, partner.coord
            )
            if angle is not None and angle <= settings.max_donor_deviation_deg:
                evidence = FixedSiteEvidence(
                    variable_site=variable.atom_key,
                    required_role="donor",
                    partner_site=partner.atom_key,
                    distance_angstrom=distance,
                    donor_deviation_deg=angle,
                )
                evidence_by_site_role[(variable.atom_key, "donor")].append(evidence)
                report.fixed_evidence.append(evidence)
        elif partner_is_donor and not partner_is_acceptor:
            angle = _best_fixed_donor_angle(
                molecule, building_template, partner, variable.coord
            )
            if angle is not None and angle <= settings.max_donor_deviation_deg:
                evidence = FixedSiteEvidence(
                    variable_site=variable.atom_key,
                    required_role="acceptor",
                    partner_site=partner.atom_key,
                    distance_angstrom=distance,
                    donor_deviation_deg=angle,
                )
                evidence_by_site_role[(variable.atom_key, "acceptor")].append(evidence)
                report.fixed_evidence.append(evidence)
        elif partner_is_donor and partner_is_acceptor:
            target_angle = _best_donor_angle_for_options(
                molecule, building_template, owner, variable.atom_name, partner.coord
            )
            partner_angle = _best_fixed_donor_angle(
                molecule, building_template, partner, variable.coord
            )
            if (
                target_angle is not None
                and target_angle <= settings.max_donor_deviation_deg
            ) or (
                partner_angle is not None
                and partner_angle <= settings.max_donor_deviation_deg
            ):
                key = _canonical_bond(variable.atom_key, partner.atom_key)
                if key not in ambivalent_keys:
                    ambivalent_keys.add(key)
                    report.ambivalent_contacts.append(
                        AmbivalentSiteContact(
                            variable_site=variable.atom_key,
                            partner_site=partner.atom_key,
                            distance_angstrom=distance,
                        )
                    )

    selected_evidence: Dict[AtomKey, FixedSiteEvidence] = {}
    all_variable_keys = {
        _atom_key(item.residue, site_name)
        for item in titratables
        for site_name in item.family.variable_sites
    }
    variable_contact_keys = {
        key for contact in variable_contacts for key in (contact.site1, contact.site2)
    }
    ambivalent_variable_keys = {
        contact.variable_site for contact in report.ambivalent_contacts
    }
    for site_key in all_variable_keys:
        donor = min(
            evidence_by_site_role.get((site_key, "donor"), []),
            key=lambda item: (item.distance_angstrom, item.donor_deviation_deg),
            default=None,
        )
        acceptor = min(
            evidence_by_site_role.get((site_key, "acceptor"), []),
            key=lambda item: (item.distance_angstrom, item.donor_deviation_deg),
            default=None,
        )
        if donor is not None and acceptor is not None:
            selected = min(
                (donor, acceptor),
                key=lambda item: (
                    item.distance_angstrom,
                    item.donor_deviation_deg,
                    0 if item.required_role == "donor" else 1,
                ),
            )
            selected_evidence[site_key] = selected
            report.conflicts.append(
                ProtonationConflict(
                    kind="opposing_fixed_evidence",
                    message=(
                        f"{site_key} has donor and acceptor evidence; selected "
                        f"{selected.required_role} from the shorter contact"
                    ),
                    sites=(site_key,),
                )
            )
            report.site_statuses[site_key] = selected.required_role
        elif donor is not None or acceptor is not None:
            selected = donor if donor is not None else acceptor
            assert selected is not None
            selected_evidence[site_key] = selected
            report.site_statuses[site_key] = selected.required_role
        elif site_key in variable_contact_keys:
            report.site_statuses[site_key] = "variable"
        elif site_key in ambivalent_variable_keys:
            report.site_statuses[site_key] = "ambivalent"
        else:
            report.site_statuses[site_key] = "unknown"
    return selected_evidence, variable_contacts


def _state_supports_role(
    option: _ProtonationStateOption,
    site_name: str,
    role: str,
) -> bool:
    return role in option.site_roles.get(site_name, frozenset())


def _component_indices(
    count: int,
    contacts: Sequence[VariableSiteContact],
    owner_by_site: Mapping[AtomKey, int],
) -> List[List[int]]:
    graph: Dict[int, Set[int]] = {index: set() for index in range(count)}
    for contact in contacts:
        left = owner_by_site[contact.site1]
        right = owner_by_site[contact.site2]
        graph[left].add(right)
        graph[right].add(left)
    components: List[List[int]] = []
    remaining = set(range(count))
    while remaining:
        start = min(remaining)
        stack = [start]
        component: List[int] = []
        remaining.remove(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in sorted(graph[node]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def _solve_protonation_states(
    titratables: Sequence[_TitratableResidue],
    fixed_requirements: Mapping[AtomKey, FixedSiteEvidence],
    contacts: Sequence[VariableSiteContact],
    report: ProtonationAssignmentReport,
) -> Dict[int, _ProtonationStateOption]:
    owner_by_site: Dict[AtomKey, int] = {}
    site_name_by_key: Dict[AtomKey, str] = {}
    for item in titratables:
        for site_name in item.family.variable_sites:
            key = _atom_key(item.residue, site_name)
            owner_by_site[key] = item.index
            site_name_by_key[key] = site_name

    chosen: Dict[int, _ProtonationStateOption] = {}
    for component in _component_indices(len(titratables), contacts, owner_by_site):
        component_set = set(component)
        component_contacts = [
            contact
            for contact in contacts
            if owner_by_site[contact.site1] in component_set
            and owner_by_site[contact.site2] in component_set
        ]
        option_lists = [titratables[index].family.options for index in component]
        best_score = None
        best_deterministic = None
        best_combo = None
        for combo in itertools.product(*option_lists):
            option_by_index = dict(zip(component, combo))
            fixed_satisfied = 0
            fixed_distance = 0.0
            for site_key, evidence in fixed_requirements.items():
                owner_index = owner_by_site[site_key]
                if owner_index not in component_set:
                    continue
                option = option_by_index[owner_index]
                if _state_supports_role(
                    option, site_name_by_key[site_key], evidence.required_role
                ):
                    fixed_satisfied += 1
                    fixed_distance += evidence.distance_angstrom
            variable_satisfied = 0
            variable_distance = 0.0
            for contact in component_contacts:
                left_index = owner_by_site[contact.site1]
                right_index = owner_by_site[contact.site2]
                left_roles = option_by_index[left_index].site_roles.get(
                    site_name_by_key[contact.site1], frozenset()
                )
                right_roles = option_by_index[right_index].site_roles.get(
                    site_name_by_key[contact.site2], frozenset()
                )
                if any(
                    left_role in left_roles and right_role in right_roles
                    for left_role, right_role in contact.allowed_role_pairs
                ):
                    variable_satisfied += 1
                    variable_distance += contact.distance_angstrom
            defaults = sum(
                option.resname == titratables[index].family.default_option.resname
                and option.mol_type == titratables[index].family.default_option.mol_type
                for index, option in option_by_index.items()
            )
            level_penalty = sum(
                abs(option.level_index - titratables[index].family.default_level_index)
                for index, option in option_by_index.items()
            )
            tautomer_penalty = sum(option.tautomer_index for option in combo)
            score = (
                fixed_satisfied,
                variable_satisfied,
                -(fixed_distance + variable_distance),
                defaults,
                -level_penalty,
                -tautomer_penalty,
            )
            deterministic = tuple(
                (option.level_index, option.tautomer_index, option.resname)
                for option in combo
            )
            if (
                best_score is None
                or score > best_score
                or (score == best_score and deterministic < best_deterministic)
            ):
                best_score = score
                best_deterministic = deterministic
                best_combo = combo
        assert best_combo is not None
        chosen.update(dict(zip(component, best_combo)))

    for site_key, evidence in fixed_requirements.items():
        option = chosen[owner_by_site[site_key]]
        if not _state_supports_role(
            option, site_name_by_key[site_key], evidence.required_role
        ):
            report.conflicts.append(
                ProtonationConflict(
                    kind="unsatisfied_fixed_evidence",
                    message=(
                        f"Selected {option.resname} cannot make {site_key} "
                        f"a required {evidence.required_role}"
                    ),
                    sites=(site_key,),
                )
            )
    for contact in contacts:
        left_index = owner_by_site[contact.site1]
        right_index = owner_by_site[contact.site2]
        left = chosen[left_index]
        right = chosen[right_index]
        left_roles = left.site_roles.get(site_name_by_key[contact.site1], frozenset())
        right_roles = right.site_roles.get(site_name_by_key[contact.site2], frozenset())
        if not any(
            role1 in left_roles and role2 in right_roles
            for role1, role2 in contact.allowed_role_pairs
        ):
            report.conflicts.append(
                ProtonationConflict(
                    kind="unsatisfied_variable_contact",
                    message=(
                        f"Selected {left.resname}/{right.resname} do not satisfy "
                        f"the variable contact {contact.site1} -- {contact.site2}"
                    ),
                    sites=(contact.site1, contact.site2),
                )
            )
    return chosen


def assign_protonation_states(
    molecule: Molecule,
    converting: Mapping[str, Any],
    building_template: Mapping[str, Any],
    state_data: Mapping[str, Any],
    *,
    pH: float = 7.0,
    geometry_settings: Optional[HydrogenBondGeometrySettings] = None,
    modify_myself: bool = False,
) -> Tuple[Molecule, ProtonationAssignmentReport]:
    """Assign protonation levels and tautomers from pH and H-bond geometry."""
    if not math.isfinite(pH):
        raise ValueError("pH must be finite")
    settings = geometry_settings or HydrogenBondGeometrySettings()
    if settings.heavy_atom_cutoff_angstrom <= 0.0:
        raise ValueError("heavy_atom_cutoff_angstrom must be positive")
    if not 0.0 <= settings.max_donor_deviation_deg <= 180.0:
        raise ValueError("max_donor_deviation_deg must be between 0 and 180")

    target = molecule if modify_myself else copy.deepcopy(molecule)
    report = ProtonationAssignmentReport(pH=float(pH))
    families, by_state = _load_protonation_families(state_data, pH)
    for family in families:
        family_names = tuple(option.resname for option in family.options)
        report.family_defaults.append(
            (family_names, family.default_option.resname)
        )

    original_identity: Dict[int, Tuple[Optional[str], str]] = {}
    titratables: List[_TitratableResidue] = []
    for residue in _iter_residues(target):
        family = by_state.get((residue.group or "", residue.ff_resname))
        if family is None:
            continue
        index = len(titratables)
        original_identity[index] = (residue.group, residue.ff_resname)
        assign_residue_identity(
            residue,
            family.default_option.mol_type,
            family.default_option.resname,
            converting,
            assignment_source="protonation_ph_default",
        )
        titratables.append(
            _TitratableResidue(index=index, residue=residue, family=family)
        )

    if not titratables:
        return target, report

    fixed_requirements, variable_contacts = _analyze_hbond_constraints(
        target,
        building_template,
        titratables,
        settings,
        report,
    )
    chosen = _solve_protonation_states(
        titratables,
        fixed_requirements,
        variable_contacts,
        report,
    )

    for item in titratables:
        option = chosen[item.index]
        old_mol_type, old_resname = original_identity[item.index]
        assign_residue_identity(
            item.residue,
            option.mol_type,
            option.resname,
            converting,
            assignment_source="protonation_hbond_assignment",
        )
        default = item.family.default_option
        report.assignments.append(
            ProtonationResidueAssignment(
                residue_key=_residue_key(item.residue),
                old_mol_type=old_mol_type,
                old_resname=old_resname,
                new_mol_type=option.mol_type,
                new_resname=option.resname,
                default_resname=default.resname,
                is_default=(
                    option.mol_type == default.mol_type
                    and option.resname == default.resname
                ),
            )
        )

    if report.conflicts:
        warning = (
            f"Protonation assignment completed with {len(report.conflicts)} "
            "conflicting or unsatisfied hydrogen-bond constraints"
        )
        report.warnings.append(warning)
        target.warnings.append(warning)
    return target, report


def assign_molecule_states(
    molecule: Molecule,
    converting: Mapping[str, Any],
    building_template: Mapping[str, Any],
    state_data: Mapping[str, Any],
    *,
    pH: float = 7.0,
    covalent_cutoff_angstrom: float = 2.3,
    hydrogen_bond_settings: Optional[HydrogenBondGeometrySettings] = None,
    modify_myself: bool = False,
) -> Tuple[Molecule, StateAssignmentReport]:
    """Run covalent assignment followed by protonation/tautomer assignment."""
    target = molecule if modify_myself else copy.deepcopy(molecule)
    target, covalent_report = assign_covalent_states(
        target,
        converting,
        state_data,
        cutoff_angstrom=covalent_cutoff_angstrom,
        modify_myself=True,
    )
    target, protonation_report = assign_protonation_states(
        target,
        converting,
        building_template,
        state_data,
        pH=pH,
        geometry_settings=hydrogen_bond_settings,
        modify_myself=True,
    )
    return target, StateAssignmentReport(
        covalent=covalent_report,
        protonation=protonation_report,
    )


def summarize_covalent_assignment(
    report: CovalentStateAssignmentReport,
) -> Dict[str, int]:
    """Return a compact covalent-assignment summary."""
    return {
        "covalent_bonds": len(report.bonds),
        "state_changes": len(report.state_changes),
        "missing_bond_atoms": len(report.missing_bond_atoms),
    }


def summarize_state_assignment(
    report: CovalentStateAssignmentReport,
) -> Dict[str, int]:
    """Backward-compatible alias for the original covalent-only summary."""
    return summarize_covalent_assignment(report)


def summarize_protonation_assignment(
    report: ProtonationAssignmentReport,
) -> Dict[str, int]:
    return {
        "assigned_residues": len(report.assignments),
        "nondefault_residues": sum(not item.is_default for item in report.assignments),
        "fixed_evidence": len(report.fixed_evidence),
        "variable_contacts": len(report.variable_contacts),
        "ambivalent_contacts": len(report.ambivalent_contacts),
        "unevaluable_sites": len(report.unevaluable_sites),
        "conflicts": len(report.conflicts),
        "warnings": len(report.warnings),
    }


def summarize_molecule_state_assignment(
    report: StateAssignmentReport,
) -> Dict[str, Dict[str, int]]:
    return {
        "covalent": summarize_covalent_assignment(report.covalent),
        "protonation": summarize_protonation_assignment(report.protonation),
    }
