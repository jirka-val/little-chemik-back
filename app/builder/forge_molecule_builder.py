#!/usr/bin/env python3
"""Planning and coordinate execution for deterministic missing-atom building.

The module contains:
  * a symbolic planner that converts missing atoms into build steps, including
    missing-DOF boundaries;
  * a coordinate executor that applies build steps until the next missing DOF;
  * molecule/build-plan adapters for the reusable local MM implementation in
    ``forge_molecule_mm`` (free-rotor H and missing-sidechain optimization).
"""

from __future__ import annotations

import heapq
import time
from collections import defaultdict, deque
from dataclasses import dataclass, field, replace
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

from forge_molecule_parser import Molecule, Residue, Chain
from forge_molecule_mm import (
    COULOMB_KJ_MOL_NM,
    FreeRotorSearchSettings,
    MMAtomParams,
    MMParameterProvider,
    MMTorsionTerm,
    MMSpatialIndex,
    PeriodicOptimizationResult,
    SidechainOptimizationSettings,
    SpatialAtom,
    optimize_periodic_dofs,
    optimize_symmetric_periodic_1d,
    pair_nonbonded_energy,
    proper_torsion_energy,
    refine_periodic_dofs,
    switch_weights,
    topological_shells_upto_three,
    vectorized_nonbonded_energy,
)


TORSION_CLASS_PRIORITY = {
    "rigid": 0,
    "derived_rotatable": 1,
    "free_rotor_hydrogen": 2,
}


@dataclass(frozen=True, order=True)
class AtomID:
    chain_id: str
    residue_index: int
    atom_name: str


@dataclass(frozen=True)
class ResolvedRef:
    chain_id: str
    residue_index: int
    atom_name: str

    def atom_id(self) -> AtomID:
        return AtomID(self.chain_id, self.residue_index, self.atom_name)


@dataclass(frozen=True)
class DOFKey:
    chain_id: str
    residue_index: int
    atom_name: str
    rule_index: int


@dataclass(frozen=True, order=True)
class ResidueID:
    chain_id: str
    residue_index: int


@dataclass
class PlannedMissingDOFStep:
    dof_key: DOFKey
    central_bond: Tuple[ResolvedRef, ResolvedRef]
    requested_dihedral_atoms: Tuple[ResolvedRef, ResolvedRef, ResolvedRef, ResolvedRef]
    torsion_group_index: int
    requested_member_index: int
    reason_atom: AtomID
    reason_rule_index: int
    prerequisite_dofs: Tuple[DOFKey, ...] = ()
    local_completion_group: Optional[ResidueID] = None
    completion_classification: str = "unclassified"


@dataclass
class PlannedAtomBuildStep:
    atom_key: AtomID
    rule_index: int
    torsion_class: str
    torsion_source: str
    # torsion_source values currently planned:
    #   internal                 rigid or ordinary fallback
    #   observed_member          derived_rotatable from a template torsion-group member
    #   supplied_dof             derived_rotatable from a PlannedMissingDOFStep
    #   free_rotor_search        first H in a free-rotor group
    #   free_rotor_group_phase   subsequent H in a free-rotor group
    observed_member_index: Optional[int] = None
    dof_key: Optional[DOFKey] = None
    free_rotor_anchor: Optional[AtomID] = None
    phase_offset: Optional[float] = None
    free_rotor_group_atoms: Tuple[str, ...] = ()
    required_dofs: Tuple[DOFKey, ...] = ()
    local_completion_group: Optional[ResidueID] = None
    torsion_group_index: Optional[int] = None
    target_member_index: Optional[int] = None


PlanStep = Union[PlannedAtomBuildStep, PlannedMissingDOFStep]


@dataclass(frozen=True)
class ResidueLocalCompletionGroup:
    residue_id: ResidueID
    classification: str
    dof_keys: Tuple[DOFKey, ...]
    atom_keys: Tuple[AtomID, ...]
    anchor_atoms: Tuple[AtomID, ...]
    non_bridge_dofs: Tuple[DOFKey, ...] = ()
    external_pending_atoms: Tuple[AtomID, ...] = ()


@dataclass
class BuildableCandidate:
    atom_key: AtomID
    rule_index: int
    torsion_class: str
    torsion_source: str
    observed_member_index: Optional[int] = None
    dof_key: Optional[DOFKey] = None
    torsion_group_index: Optional[int] = None
    target_member_index: Optional[int] = None


@dataclass
class MissingDOFCandidate:
    atom_key: AtomID
    rule_index: int
    dof_key: DOFKey
    central_bond: Tuple[ResolvedRef, ResolvedRef]
    requested_dihedral_atoms: Tuple[ResolvedRef, ResolvedRef, ResolvedRef, ResolvedRef]
    torsion_group_index: int
    requested_member_index: int


@dataclass
class AtomAvailability:
    atom_key: AtomID
    buildable: Optional[BuildableCandidate] = None
    missing_dof: Optional[MissingDOFCandidate] = None
    missing_refs: List[ResolvedRef] = field(default_factory=list)
    no_rules: bool = False
    error: Optional[str] = None


@dataclass
class BuildPlan:
    steps: List[PlanStep]
    unresolved_atoms: Set[AtomID]
    requirements: List[PlannedMissingDOFStep]
    local_completions: Tuple[ResidueLocalCompletionGroup, ...] = ()


@dataclass(frozen=True)
class SidechainExecutionIndex:
    """Immutable execution view distilled from one remaining build plan."""

    residue_groups: Tuple[ResidueID, ...]
    dof_steps: Tuple[PlannedMissingDOFStep, ...]
    atom_steps: Tuple[PlannedAtomBuildStep, ...]
    sidechain_steps: Tuple[PlanStep, ...]
    residual_steps: Tuple[PlanStep, ...]
    affected_atom_steps_by_dof: Mapping[
        DOFKey,
        Tuple[PlannedAtomBuildStep, ...],
    ]


@dataclass
class SidechainMMOptimizationStats:
    """Optional diagnostics populated by ``optimize_missing_sidechain_dofs``."""

    branch_count: int = 0
    cluster_count: int = 0
    coupled_branch_pairs: int = 0
    coarse_conformations: int = 0
    energy_evaluations: int = 0
    refinement_sweeps: int = 0
    neighbor_rebuilds: int = 0
    preparation_seconds: float = 0.0
    coarse_sampling_seconds: float = 0.0
    clustering_seconds: float = 0.0
    optimization_seconds: float = 0.0
    total_seconds: float = 0.0


@dataclass(frozen=True)
class SidechainLocalOptimizationContext:
    """Opaque reusable MM context for the GUI ``Opt`` action of one residue.

    The selected residue is dynamic; every other atom is frozen at the
    coordinates present when the context is prepared.  Reuse is therefore safe
    across slider changes of this residue, but the context must be rebuilt
    after coordinates outside this residue change.
    """

    residue_id: ResidueID
    dof_order: Tuple[DOFKey, ...]
    affected_atoms_by_dof: Mapping[DOFKey, Tuple[AtomID, ...]]
    coordinates: Any = field(repr=False, compare=False)
    energy_model: Any = field(repr=False, compare=False)
    settings: SidechainOptimizationSettings = field(repr=False, compare=False)


@dataclass(frozen=True, order=True)
class TargetRuleID:
    target_atom: AtomID
    rule_index: int


@dataclass
class PlannerDependencyIndex:
    atom_to_rules: Dict[AtomID, Set[TargetRuleID]]
    rigid_atom_to_rules: Dict[AtomID, Set[TargetRuleID]]
    target_to_rigid_rules: Dict[AtomID, Tuple[TargetRuleID, ...]]
    dof_to_rule: Dict[DOFKey, TargetRuleID]

    def affected_targets(self, atom_id: AtomID) -> Set[AtomID]:
        return {
            rule_id.target_atom
            for rule_id in self.atom_to_rules.get(atom_id, set())
        }

    def rigid_affected_targets(self, atom_id: AtomID) -> Set[AtomID]:
        return {
            rule_id.target_atom
            for rule_id in self.rigid_atom_to_rules.get(atom_id, set())
        }


@dataclass
class PlannerStats:
    initial_pending_atoms: int = 0
    dependency_edges: int = 0
    full_atom_evaluations: int = 0
    rigid_rule_evaluations: int = 0
    dependency_target_updates: int = 0
    raw_planning_seconds: float = 0.0
    local_completion_analysis_seconds: float = 0.0
    total_planning_seconds: float = 0.0


class PlanningError(RuntimeError):
    pass


# -----------------------------------------------------------------------------
# Molecule helpers
# -----------------------------------------------------------------------------

def chain_order(mol: Molecule) -> Dict[str, int]:
    return {cid: i for i, cid in enumerate(mol.chains.keys())}


def get_residue(mol: Molecule, atom_id: AtomID) -> Residue:
    return mol.chains[atom_id.chain_id].residues[atom_id.residue_index]


def atom_order_in_residue(residue: Residue) -> Dict[str, int]:
    return {name: i for i, name in enumerate(residue.atoms.keys())}


def atom_sort_key(mol: Molecule, atom_id: AtomID) -> Tuple[int, int, int, str]:
    c_order = chain_order(mol).get(atom_id.chain_id, 10**9)
    residue = get_residue(mol, atom_id)
    a_order = atom_order_in_residue(residue).get(atom_id.atom_name, 10**9)
    return (c_order, atom_id.residue_index, a_order, atom_id.atom_name)


def collect_available_atoms(mol: Molecule) -> Set[AtomID]:
    available: Set[AtomID] = set()
    for cid, chain in mol.chains.items():
        for res in chain.residues:
            for atom_name, atom in res.atoms.items():
                if atom.coord is not None:
                    available.add(AtomID(cid, res.index_in_chain, atom_name))
    return available


def collect_pending_atoms(mol: Molecule) -> Set[AtomID]:
    pending: Set[AtomID] = set()
    for cid, chain in mol.chains.items():
        for res in chain.residues:
            for atom_name, atom in res.atoms.items():
                if atom.coord is None:
                    pending.add(AtomID(cid, res.index_in_chain, atom_name))
    return pending


def torsion_class(rule: Mapping[str, Any]) -> str:
    torsion = rule.get("torsion", {})
    if isinstance(torsion, str):
        return torsion
    return str(torsion.get("class"))


def torsion_data(rule: Mapping[str, Any]) -> Mapping[str, Any]:
    torsion = rule.get("torsion", {})
    return torsion if isinstance(torsion, Mapping) else {"class": torsion}


def normalize_angle(angle: float) -> float:
    out = ((angle + 180.0) % 360.0) - 180.0
    # Prefer +180 over -180 only if it came from positive 180-ish input; irrelevant for planning.
    if out == -180.0:
        return 180.0
    return out


# -----------------------------------------------------------------------------
# Template/reference helpers
# -----------------------------------------------------------------------------

def residue_template(template: Mapping[str, Any], residue: Residue) -> Mapping[str, Any]:
    if not residue.group or residue.group not in template:
        raise PlanningError(f"No template mol_type {residue.group!r} for residue {residue.chain_id}:{residue.resseq}{residue.icode} {residue.ff_resname}")
    if residue.ff_resname not in template[residue.group]:
        raise PlanningError(f"No template for residue {residue.chain_id}:{residue.resseq}{residue.icode} {residue.ff_resname} in mol_type {residue.group}")
    return template[residue.group][residue.ff_resname]


def resolve_ref(mol: Molecule, residue: Residue, ref: Mapping[str, Any]) -> Optional[ResolvedRef]:
    offset = int(ref.get("residue_offset", 0))
    atom_name = str(ref.get("atom"))
    chain = mol.chains[residue.chain_id]
    idx = residue.index_in_chain + offset
    if idx < 0 or idx >= len(chain.residues):
        return None
    target_res = chain.residues[idx]
    if atom_name not in target_res.atoms:
        return None
    return ResolvedRef(residue.chain_id, idx, atom_name)


def require_resolved_ref(mol: Molecule, residue: Residue, ref: Mapping[str, Any], context: str) -> ResolvedRef:
    resolved = resolve_ref(mol, residue, ref)
    if resolved is None:
        raise PlanningError(f"Template reference cannot be resolved in {context}: offset={ref.get('residue_offset',0)} atom={ref.get('atom')}")
    return resolved


def resolve_refs(
    mol: Molecule,
    residue: Residue,
    refs: Sequence[Mapping[str, Any]],
) -> Optional[List[ResolvedRef]]:
    """Resolve a complete optional reference set in one concrete sequence.

    Protein templates intentionally contain rules for several possible neighbor
    residue contexts (for example, generic and proline-specific peptide rules).
    A rule is contextually usable only when all of its references resolve.
    """
    resolved = [resolve_ref(mol, residue, ref) for ref in refs]
    if any(ref is None for ref in resolved):
        return None
    return [ref for ref in resolved if ref is not None]


def _same_atom(a: ResolvedRef, b: ResolvedRef) -> bool:
    return a == b


def find_matching_torsion_group(
    mol: Molecule,
    residue: Residue,
    res_template: Mapping[str, Any],
    rule: Mapping[str, Any],
    target: AtomID,
) -> Tuple[int, int, Mapping[str, Any], bool]:
    """Return (group_index, target_member_index, group, reversed_central)."""
    refs = resolve_refs(mol, residue, rule.get("refs", []))
    if refs is None:
        raise PlanningError(
            f"Contextually unavailable derived rule was selected for {target}"
        )
    if len(refs) != 3:
        raise PlanningError(f"Rule for {target} does not have exactly three refs")
    ref0, ref1, ref2 = refs
    target_ref = ResolvedRef(target.chain_id, target.residue_index, target.atom_name)

    for gi, group in enumerate(res_template.get("torsion_groups", [])):
        cb = resolve_refs(mol, residue, group.get("central_bond", []))
        if cb is None or len(cb) != 2:
            continue
        same = _same_atom(cb[0], ref1) and _same_atom(cb[1], ref2)
        rev = _same_atom(cb[0], ref2) and _same_atom(cb[1], ref1)
        if not (same or rev):
            continue
        for mi, member in enumerate(group.get("members", [])):
            terms = resolve_refs(mol, residue, member.get("terminal_atoms", []))
            if terms is None or len(terms) != 2:
                continue
            if same and _same_atom(terms[0], ref0) and _same_atom(terms[1], target_ref):
                return gi, mi, group, False
            if rev and _same_atom(terms[0], target_ref) and _same_atom(terms[1], ref0):
                return gi, mi, group, True
    raise PlanningError(f"No matching torsion group/member for derived rule {target} in {residue.ff_resname}")


def member_dihedral_atoms(
    mol: Molecule,
    residue: Residue,
    group: Mapping[str, Any],
    member_index: int,
) -> Optional[Tuple[ResolvedRef, ResolvedRef, ResolvedRef, ResolvedRef]]:
    cb = resolve_refs(mol, residue, group.get("central_bond", []))
    if cb is None or len(cb) != 2:
        return None
    member = group.get("members", [])[member_index]
    terms = resolve_refs(mol, residue, member.get("terminal_atoms", []))
    if terms is None or len(terms) != 2:
        return None
    return (terms[0], cb[0], cb[1], terms[1])


def first_available_member_index(
    mol: Molecule,
    residue: Residue,
    group: Mapping[str, Any],
    available_atoms: Set[AtomID],
) -> Optional[int]:
    for mi, _member in enumerate(group.get("members", [])):
        atoms = member_dihedral_atoms(mol, residue, group, mi)
        if atoms is not None and all(a.atom_id() in available_atoms for a in atoms):
            return mi
    return None


# -----------------------------------------------------------------------------
# Rule evaluation
# -----------------------------------------------------------------------------

def evaluate_rule(
    mol: Molecule,
    template: Mapping[str, Any],
    atom_key: AtomID,
    rule: Mapping[str, Any],
    rule_index: int,
    available_atoms: Set[AtomID],
    available_dofs: Set[DOFKey],
) -> Tuple[str, Optional[BuildableCandidate], Optional[MissingDOFCandidate], List[ResolvedRef]]:
    """Evaluate one rule.

    Returns status, buildable, missing_dof, missing_refs.
    status values: buildable, missing_dof, missing_refs.
    """
    residue = get_residue(mol, atom_key)
    res_template = residue_template(template, residue)
    refs = resolve_refs(mol, residue, rule.get("refs", []))
    if refs is None:
        return "missing_refs", None, None, []
    if len(refs) != 3:
        raise PlanningError(f"Rule {atom_key}[{rule_index}] does not have exactly three refs")

    missing = [r for r in refs if r.atom_id() not in available_atoms]
    if missing:
        return "missing_refs", None, None, missing

    cls = torsion_class(rule)
    if cls == "rigid":
        return "buildable", BuildableCandidate(atom_key, rule_index, cls, "internal"), None, []

    if cls == "free_rotor_hydrogen":
        return "buildable", BuildableCandidate(atom_key, rule_index, cls, "free_rotor_search"), None, []

    if cls != "derived_rotatable":
        raise PlanningError(f"Unsupported torsion class {cls!r} for {atom_key}[{rule_index}]")

    gi, target_mi, group, _rev = find_matching_torsion_group(mol, residue, res_template, rule, atom_key)
    observed_mi = first_available_member_index(mol, residue, group, available_atoms)
    if observed_mi is not None:
        return "buildable", BuildableCandidate(
            atom_key=atom_key,
            rule_index=rule_index,
            torsion_class=cls,
            torsion_source="observed_member",
            observed_member_index=observed_mi,
            torsion_group_index=gi,
            target_member_index=target_mi,
        ), None, []

    dof_key = DOFKey(atom_key.chain_id, atom_key.residue_index, atom_key.atom_name, rule_index)
    if dof_key in available_dofs:
        return "buildable", BuildableCandidate(
            atom_key=atom_key,
            rule_index=rule_index,
            torsion_class=cls,
            torsion_source="supplied_dof",
            dof_key=dof_key,
            torsion_group_index=gi,
            target_member_index=target_mi,
        ), None, []

    requested_atoms = member_dihedral_atoms(mol, residue, group, target_mi)
    if requested_atoms is None:
        raise PlanningError(
            f"Target torsion member cannot be resolved for {atom_key}[{rule_index}]"
        )
    cb = (requested_atoms[1], requested_atoms[2])
    missing_dof = MissingDOFCandidate(
        atom_key=atom_key,
        rule_index=rule_index,
        dof_key=dof_key,
        central_bond=cb,
        requested_dihedral_atoms=requested_atoms,
        torsion_group_index=gi,
        requested_member_index=target_mi,
    )
    return "missing_dof", None, missing_dof, []


def evaluate_atom(
    mol: Molecule,
    template: Mapping[str, Any],
    atom_key: AtomID,
    available_atoms: Set[AtomID],
    available_dofs: Set[DOFKey],
) -> AtomAvailability:
    residue = get_residue(mol, atom_key)
    rt = residue_template(template, residue)
    rules = rt.get("build_rules", {}).get(atom_key.atom_name, [])
    if not rules:
        return AtomAvailability(atom_key=atom_key, no_rules=True, error="no_rules")

    buildable_candidates: List[BuildableCandidate] = []
    dof_candidates: List[MissingDOFCandidate] = []
    all_missing_refs: List[ResolvedRef] = []

    for ri, rule in enumerate(rules):
        status, buildable, missing_dof, missing_refs = evaluate_rule(
            mol, template, atom_key, rule, ri, available_atoms, available_dofs
        )
        if status == "buildable" and buildable is not None:
            buildable_candidates.append(buildable)
        elif status == "missing_dof" and missing_dof is not None:
            dof_candidates.append(missing_dof)
        elif status == "missing_refs":
            all_missing_refs.extend(missing_refs)

    if buildable_candidates:
        best = min(buildable_candidates, key=lambda c: (
            TORSION_CLASS_PRIORITY.get(c.torsion_class, 999),
            c.rule_index,
        ))
        return AtomAvailability(atom_key=atom_key, buildable=best, missing_refs=all_missing_refs)

    if dof_candidates:
        # Missing DOF priority follows torsion-group member priority, then rule priority.
        best_dof = min(dof_candidates, key=lambda d: (d.requested_member_index, d.rule_index))
        return AtomAvailability(atom_key=atom_key, missing_dof=best_dof, missing_refs=all_missing_refs)

    return AtomAvailability(atom_key=atom_key, missing_refs=all_missing_refs)


def evaluate_pending_atoms(
    mol: Molecule,
    template: Mapping[str, Any],
    pending_atoms: Set[AtomID],
    available_atoms: Set[AtomID],
    available_dofs: Set[DOFKey],
) -> Dict[AtomID, AtomAvailability]:
    return {
        atom_key: evaluate_atom(mol, template, atom_key, available_atoms, available_dofs)
        for atom_key in pending_atoms
    }


def select_best_buildable(
    mol: Molecule,
    availability: Mapping[AtomID, AtomAvailability],
) -> Optional[BuildableCandidate]:
    candidates = [a.buildable for a in availability.values() if a.buildable is not None]
    if not candidates:
        return None
    return min(candidates, key=lambda c: (
        TORSION_CLASS_PRIORITY.get(c.torsion_class, 999),
        c.rule_index,
        atom_sort_key(mol, c.atom_key),
    ))


def select_best_missing_dof(
    mol: Molecule,
    availability: Mapping[AtomID, AtomAvailability],
) -> Optional[MissingDOFCandidate]:
    candidates = [a.missing_dof for a in availability.values() if a.missing_dof is not None]
    if not candidates:
        return None
    return min(candidates, key=lambda d: (
        d.requested_member_index,
        d.rule_index,
        atom_sort_key(mol, d.atom_key),
    ))


def build_planner_dependency_index(
    mol: Molecule,
    template: Mapping[str, Any],
    pending_atoms: Set[AtomID],
) -> PlannerDependencyIndex:
    """Resolve reverse rule dependencies for one concrete molecule.

    Direct rule references are indexed for all torsion classes. Derived rules
    additionally depend on every atom of every member in their matching torsion
    group, because any one of those atoms can make an observed member available.
    """
    atom_to_rules: Dict[AtomID, Set[TargetRuleID]] = defaultdict(set)
    rigid_atom_to_rules: Dict[AtomID, Set[TargetRuleID]] = defaultdict(set)
    target_to_rigid_rules: Dict[AtomID, List[TargetRuleID]] = defaultdict(list)
    dof_to_rule: Dict[DOFKey, TargetRuleID] = {}

    for target in sorted(pending_atoms, key=lambda atom: atom_sort_key(mol, atom)):
        residue = get_residue(mol, target)
        rt = residue_template(template, residue)
        rules = rt.get("build_rules", {}).get(target.atom_name, [])
        for rule_index, rule in enumerate(rules):
            rule_id = TargetRuleID(target, rule_index)
            cls = torsion_class(rule)
            refs = resolve_refs(mol, residue, rule.get("refs", []))
            if refs is None:
                # Optional rule for another valid neighbor-residue context.
                continue
            if len(refs) != 3:
                raise PlanningError(
                    f"Rule {target}[{rule_index}] does not have exactly three refs"
                )
            for ref in refs:
                atom_to_rules[ref.atom_id()].add(rule_id)
                if cls == "rigid":
                    rigid_atom_to_rules[ref.atom_id()].add(rule_id)

            if cls == "rigid":
                target_to_rigid_rules[target].append(rule_id)
            elif cls == "derived_rotatable":
                group_index, _target_member, group, _reversed = (
                    find_matching_torsion_group(
                        mol,
                        residue,
                        rt,
                        rule,
                        target,
                    )
                )
                for member_index, _member in enumerate(group.get("members", [])):
                    member_atoms = member_dihedral_atoms(
                        mol,
                        residue,
                        group,
                        member_index,
                    )
                    if member_atoms is None:
                        continue
                    for member_atom in member_atoms:
                        atom_to_rules[member_atom.atom_id()].add(rule_id)
                dof_to_rule[
                    DOFKey(
                        target.chain_id,
                        target.residue_index,
                        target.atom_name,
                        rule_index,
                    )
                ] = rule_id
            elif cls != "free_rotor_hydrogen":
                raise PlanningError(
                    f"Unsupported torsion class {cls!r} for "
                    f"{target}[{rule_index}]"
                )

    return PlannerDependencyIndex(
        atom_to_rules=dict(atom_to_rules),
        rigid_atom_to_rules=dict(rigid_atom_to_rules),
        target_to_rigid_rules={
            target: tuple(sorted(rule_ids, key=lambda item: item.rule_index))
            for target, rule_ids in target_to_rigid_rules.items()
        },
        dof_to_rule=dof_to_rule,
    )


# -----------------------------------------------------------------------------
# Free-rotor hydrogen group planning
# -----------------------------------------------------------------------------

def _free_rotor_group_atoms(rule: Mapping[str, Any]) -> Tuple[str, ...]:
    data = torsion_data(rule)
    atoms = data.get("group_atoms") or []
    return tuple(str(a) for a in atoms)


def _find_buildable_free_rotor_rule(
    mol: Molecule,
    template: Mapping[str, Any],
    atom_key: AtomID,
    available_atoms: Set[AtomID],
    available_dofs: Set[DOFKey],
) -> Optional[int]:
    residue = get_residue(mol, atom_key)
    rt = residue_template(template, residue)
    for ri, rule in enumerate(rt.get("build_rules", {}).get(atom_key.atom_name, [])):
        if torsion_class(rule) != "free_rotor_hydrogen":
            continue
        status, buildable, _dof, _missing = evaluate_rule(mol, template, atom_key, rule, ri, available_atoms, available_dofs)
        if status == "buildable" and buildable is not None:
            return ri
    return None


def plan_free_rotor_group_steps(
    mol: Molecule,
    template: Mapping[str, Any],
    selected: BuildableCandidate,
    pending_atoms: Set[AtomID],
    available_atoms: Set[AtomID],
    available_dofs: Set[DOFKey],
) -> List[PlannedAtomBuildStep]:
    residue = get_residue(mol, selected.atom_key)
    rt = residue_template(template, residue)
    rule = rt.get("build_rules", {})[selected.atom_key.atom_name][selected.rule_index]
    group_names = _free_rotor_group_atoms(rule) or (selected.atom_key.atom_name,)
    n = len(group_names)
    step_angle = 360.0 / n if n else 0.0

    group_keys: List[AtomID] = [
        AtomID(selected.atom_key.chain_id, selected.atom_key.residue_index, name)
        for name in group_names
        if name in residue.atoms
    ]
    group_index = {key: i for i, key in enumerate(group_keys)}
    missing_group = [key for key in group_keys if key in pending_atoms]
    available_group = [key for key in group_keys if key in available_atoms]

    if not missing_group:
        return []

    steps: List[PlannedAtomBuildStep] = []

    if available_group:
        anchor = min(available_group, key=lambda k: group_index[k])
        anchor_i = group_index[anchor]
        for key in sorted(missing_group, key=lambda k: group_index[k]):
            ri = _find_buildable_free_rotor_rule(mol, template, key, available_atoms, available_dofs)
            if ri is None:
                raise PlanningError(f"Free-rotor group member {key} has no buildable free-rotor rule")
            phase = normalize_angle((group_index[key] - anchor_i) * step_angle)
            steps.append(PlannedAtomBuildStep(
                atom_key=key,
                rule_index=ri,
                torsion_class="free_rotor_hydrogen",
                torsion_source="free_rotor_group_phase",
                free_rotor_anchor=anchor,
                phase_offset=phase,
                free_rotor_group_atoms=group_names,
            ))
        return steps

    # No H from the group exists yet: choose the first missing H in group order
    # that has a buildable free-rotor rule. This H is placed by search/fallback;
    # the rest are placed by phase offsets relative to it.
    anchor: Optional[AtomID] = None
    anchor_rule_index: Optional[int] = None
    for key in sorted(missing_group, key=lambda k: group_index[k]):
        ri = _find_buildable_free_rotor_rule(mol, template, key, available_atoms, available_dofs)
        if ri is not None:
            anchor = key
            anchor_rule_index = ri
            break
    if anchor is None or anchor_rule_index is None:
        raise PlanningError(f"No buildable anchor found for free-rotor group {group_names} in {residue.ff_resname}")

    anchor_i = group_index[anchor]
    steps.append(PlannedAtomBuildStep(
        atom_key=anchor,
        rule_index=anchor_rule_index,
        torsion_class="free_rotor_hydrogen",
        torsion_source="free_rotor_search",
        free_rotor_anchor=None,
        phase_offset=0.0,
        free_rotor_group_atoms=group_names,
    ))

    # Treat anchor as available for the purpose of defining the local group phase.
    for key in sorted(missing_group, key=lambda k: group_index[k]):
        if key == anchor:
            continue
        ri = _find_buildable_free_rotor_rule(mol, template, key, available_atoms, available_dofs)
        if ri is None:
            raise PlanningError(f"Free-rotor group member {key} has no buildable free-rotor rule")
        phase = normalize_angle((group_index[key] - anchor_i) * step_angle)
        steps.append(PlannedAtomBuildStep(
            atom_key=key,
            rule_index=ri,
            torsion_class="free_rotor_hydrogen",
            torsion_source="free_rotor_group_phase",
            free_rotor_anchor=anchor,
            phase_offset=phase,
            free_rotor_group_atoms=group_names,
        ))
    return steps


# -----------------------------------------------------------------------------
# Missing-DOF provenance and residue-local completion analysis
# -----------------------------------------------------------------------------

def _dof_sort_key(dof_key: DOFKey) -> Tuple[str, int, str, int]:
    return (
        dof_key.chain_id,
        dof_key.residue_index,
        dof_key.atom_name,
        dof_key.rule_index,
    )


def _sorted_dofs(dofs: Iterable[DOFKey]) -> Tuple[DOFKey, ...]:
    return tuple(sorted(set(dofs), key=_dof_sort_key))


def _selected_rule(
    mol: Molecule,
    template: Mapping[str, Any],
    atom_key: AtomID,
    rule_index: int,
) -> Tuple[Residue, Mapping[str, Any], Mapping[str, Any]]:
    residue = get_residue(mol, atom_key)
    rt = residue_template(template, residue)
    rules = rt.get("build_rules", {}).get(atom_key.atom_name, [])
    try:
        rule = rules[rule_index]
    except IndexError as exc:
        raise PlanningError(
            f"Plan references missing rule {atom_key}[{rule_index}]"
        ) from exc
    return residue, rt, rule


def _annotate_dof_dependencies(
    mol: Molecule,
    template: Mapping[str, Any],
    steps: Sequence[PlanStep],
) -> Tuple[
    Dict[AtomID, frozenset[DOFKey]],
    Dict[DOFKey, Set[AtomID]],
    Dict[AtomID, PlannedAtomBuildStep],
    Dict[DOFKey, PlannedMissingDOFStep],
]:
    """Attach provenance only to the suffix following the first missing DOF."""
    atom_dependencies: Dict[AtomID, frozenset[DOFKey]] = {}
    dof_prerequisites: Dict[DOFKey, frozenset[DOFKey]] = {}

    first_dof_index = next(
        (
            index
            for index, step in enumerate(steps)
            if isinstance(step, PlannedMissingDOFStep)
        ),
        len(steps),
    )
    suffix = steps[first_dof_index:]
    suffix_atom_steps = {
        step.atom_key: step
        for step in suffix
        if isinstance(step, PlannedAtomBuildStep)
    }
    suffix_dof_steps = {
        step.dof_key: step
        for step in suffix
        if isinstance(step, PlannedMissingDOFStep)
    }
    affected_by_dof: Dict[DOFKey, Set[AtomID]] = {
        dof_key: set()
        for dof_key in suffix_dof_steps
    }

    def dependencies_for_atoms(
        atom_ids: Iterable[AtomID],
        context: str,
    ) -> Set[DOFKey]:
        dependencies: Set[DOFKey] = set()
        for atom_id in atom_ids:
            if (
                atom_id in suffix_atom_steps
                and atom_id not in atom_dependencies
            ):
                raise PlanningError(
                    f"Internal plan-order error: {atom_id} is unavailable in {context}"
                )
            # Atoms absent from this deliberately sparse map either had input
            # coordinates or were built before the first missing DOF. Both
            # cases have an empty missing-DOF dependency set by construction.
            dependencies.update(atom_dependencies.get(atom_id, ()))
        return dependencies

    for step in suffix:
        if isinstance(step, PlannedMissingDOFStep):
            residue, _rt, rule = _selected_rule(
                mol,
                template,
                step.reason_atom,
                step.reason_rule_index,
            )
            refs = resolve_refs(mol, residue, rule.get("refs", []))
            if refs is None or len(refs) != 3:
                raise PlanningError(
                    f"Cannot resolve reason-rule refs for missing DOF {step.dof_key}"
                )
            prerequisites = dependencies_for_atoms(
                (ref.atom_id() for ref in refs),
                f"missing DOF {step.dof_key}",
            )
            step.prerequisite_dofs = _sorted_dofs(prerequisites)
            dof_prerequisites[step.dof_key] = frozenset(prerequisites)
            continue

        residue, rt, rule = _selected_rule(
            mol,
            template,
            step.atom_key,
            step.rule_index,
        )
        refs = resolve_refs(mol, residue, rule.get("refs", []))
        if refs is None or len(refs) != 3:
            raise PlanningError(
                f"Cannot resolve build refs for planned atom {step.atom_key}"
            )
        dependencies = dependencies_for_atoms(
            (ref.atom_id() for ref in refs),
            f"build step {step.atom_key}",
        )

        if step.torsion_source == "observed_member":
            if step.observed_member_index is None:
                raise PlanningError(
                    f"Observed-member step lacks a member index: {step.atom_key}"
                )
            if step.torsion_group_index is None:
                _gi, _target_mi, group, _reversed = find_matching_torsion_group(
                    mol,
                    residue,
                    rt,
                    rule,
                    step.atom_key,
                )
            else:
                group = rt.get("torsion_groups", [])[step.torsion_group_index]
            member_atoms = member_dihedral_atoms(
                mol,
                residue,
                group,
                step.observed_member_index,
            )
            if member_atoms is None:
                raise PlanningError(
                    f"Cannot resolve observed member for {step.atom_key}"
                )
            dependencies.update(
                dependencies_for_atoms(
                    (ref.atom_id() for ref in member_atoms),
                    f"observed torsion member for {step.atom_key}",
                )
            )
        elif step.torsion_source == "supplied_dof":
            if step.dof_key is None or step.dof_key not in dof_prerequisites:
                raise PlanningError(
                    f"Supplied-DOF step has no preceding DOF: {step.atom_key}"
                )
            dependencies.update(dof_prerequisites[step.dof_key])
            dependencies.add(step.dof_key)
        elif step.torsion_source == "free_rotor_group_phase":
            if step.free_rotor_anchor is None:
                raise PlanningError(
                    f"Free-rotor group phase lacks an anchor: {step.atom_key}"
                )
            dependencies.update(
                dependencies_for_atoms(
                    (step.free_rotor_anchor,),
                    f"free-rotor group phase for {step.atom_key}",
                )
            )

        step.required_dofs = _sorted_dofs(dependencies)
        atom_dependencies[step.atom_key] = frozenset(dependencies)
        for dof_key in dependencies:
            if dof_key in affected_by_dof:
                affected_by_dof[dof_key].add(step.atom_key)

    return (
        atom_dependencies,
        affected_by_dof,
        suffix_atom_steps,
        suffix_dof_steps,
    )


@dataclass(frozen=True)
class _ResidueCovalentTopology:
    intra_neighbors: Mapping[str, Tuple[str, ...]]
    outgoing_inter: Mapping[str, Tuple[Tuple[int, str], ...]]
    incoming_inter: Mapping[Tuple[int, str], Tuple[str, ...]]


class TemplateCovalentGraphView:
    """Lazy molecule-specific view over template-derived covalent topology."""

    def __init__(self, mol: Molecule, template: Mapping[str, Any]):
        self.mol = mol
        self.template = template
        self._topologies: Dict[Tuple[str, str], _ResidueCovalentTopology] = {}
        self._neighbor_cache: Dict[AtomID, frozenset[AtomID]] = {}

    def _topology(self, residue: Residue) -> Optional[_ResidueCovalentTopology]:
        if not residue.group:
            return None
        cache_key = (str(residue.group), residue.ff_resname)
        if cache_key in self._topologies:
            return self._topologies[cache_key]
        if residue.group not in self.template:
            return None
        residues = self.template[residue.group]
        if residue.ff_resname not in residues:
            return None
        rt = residues[residue.ff_resname]
        intra: Dict[str, Set[str]] = defaultdict(set)
        outgoing: Dict[str, Set[Tuple[int, str]]] = defaultdict(set)
        incoming: Dict[Tuple[int, str], Set[str]] = defaultdict(set)
        for target_name, rules in rt.get("build_rules", {}).items():
            for rule in rules:
                refs = rule.get("refs", [])
                if len(refs) != 3:
                    continue
                bonded_ref = refs[2]
                offset = int(bonded_ref.get("residue_offset", 0))
                ref_name = str(bonded_ref.get("atom"))
                if offset == 0:
                    intra[str(target_name)].add(ref_name)
                    intra[ref_name].add(str(target_name))
                else:
                    if abs(offset) != 1:
                        raise PlanningError(
                            "Template covalent topology supports only previous/next "
                            f"residue offsets, got {offset} in {cache_key}"
                        )
                    outgoing[str(target_name)].add((offset, ref_name))
                    incoming[(offset, ref_name)].add(str(target_name))
        result = _ResidueCovalentTopology(
            intra_neighbors={
                atom: tuple(sorted(neighbors))
                for atom, neighbors in intra.items()
            },
            outgoing_inter={
                atom: tuple(sorted(neighbors))
                for atom, neighbors in outgoing.items()
            },
            incoming_inter={
                key: tuple(sorted(targets))
                for key, targets in incoming.items()
            },
        )
        self._topologies[cache_key] = result
        return result

    def neighbors(self, atom_id: AtomID) -> frozenset[AtomID]:
        cached = self._neighbor_cache.get(atom_id)
        if cached is not None:
            return cached

        residue = get_residue(self.mol, atom_id)
        topology = self._topology(residue)
        if topology is None or atom_id.atom_name not in residue.atoms:
            result = frozenset()
            self._neighbor_cache[atom_id] = result
            return result

        chain = self.mol.chains[atom_id.chain_id]
        neighbors: Set[AtomID] = set()
        for neighbor_name in topology.intra_neighbors.get(atom_id.atom_name, ()):
            if neighbor_name in residue.atoms:
                neighbors.add(
                    AtomID(
                        atom_id.chain_id,
                        atom_id.residue_index,
                        neighbor_name,
                    )
                )

        for offset, neighbor_name in topology.outgoing_inter.get(
            atom_id.atom_name,
            (),
        ):
            neighbor_index = atom_id.residue_index + offset
            if 0 <= neighbor_index < len(chain.residues):
                neighbor_residue = chain.residues[neighbor_index]
                if neighbor_name in neighbor_residue.atoms:
                    neighbors.add(
                        AtomID(atom_id.chain_id, neighbor_index, neighbor_name)
                    )

        # Add reverse inter-residue edges without materializing other atoms.
        for origin_index in (
            atom_id.residue_index - 1,
            atom_id.residue_index + 1,
        ):
            if not 0 <= origin_index < len(chain.residues):
                continue
            origin_residue = chain.residues[origin_index]
            origin_topology = self._topology(origin_residue)
            if origin_topology is None:
                continue
            offset = atom_id.residue_index - origin_index
            for target_name in origin_topology.incoming_inter.get(
                (offset, atom_id.atom_name),
                (),
            ):
                if target_name in origin_residue.atoms:
                    neighbors.add(
                        AtomID(atom_id.chain_id, origin_index, target_name)
                    )

        result = frozenset(neighbors)
        self._neighbor_cache[atom_id] = result
        return result


def build_template_covalent_graph(
    mol: Molecule,
    template: Mapping[str, Any],
) -> Dict[AtomID, Set[AtomID]]:
    """Materialize the lazy template graph for diagnostics/backward use."""
    view = TemplateCovalentGraphView(mol, template)
    graph: Dict[AtomID, Set[AtomID]] = {}
    for chain_id, chain in mol.chains.items():
        for residue in chain.residues:
            for atom_name in residue.atoms:
                atom_id = AtomID(chain_id, residue.index_in_chain, atom_name)
                graph[atom_id] = set(view.neighbors(atom_id))
    return graph


def _is_covalent_bridge(
    graph: TemplateCovalentGraphView,
    first: AtomID,
    second: AtomID,
) -> bool:
    """Test one edge lazily, expanding fixed topology only when necessary."""
    if second not in graph.neighbors(first):
        return False

    blocked_edge = frozenset((first, second))
    visited: Set[AtomID] = {second}
    queue = deque((second,))
    while queue:
        atom_id = queue.popleft()
        for neighbor in graph.neighbors(atom_id):
            if frozenset((atom_id, neighbor)) == blocked_edge:
                continue
            if neighbor == first:
                return False
            if neighbor in visited:
                continue
            visited.add(neighbor)
            queue.append(neighbor)
    return True


def _classify_local_completion_groups(
    mol: Molecule,
    atom_dependencies: Mapping[AtomID, frozenset[DOFKey]],
    affected_by_dof: Mapping[DOFKey, Set[AtomID]],
    atom_steps: Mapping[AtomID, PlannedAtomBuildStep],
    dof_steps: Mapping[DOFKey, PlannedMissingDOFStep],
    covalent_graph: TemplateCovalentGraphView,
) -> Tuple[ResidueLocalCompletionGroup, ...]:
    local_dofs: Dict[ResidueID, Set[DOFKey]] = defaultdict(set)
    for dof_key, affected_atoms in affected_by_dof.items():
        residues = {
            ResidueID(atom.chain_id, atom.residue_index)
            for atom in affected_atoms
        }
        if len(residues) == 1:
            local_dofs[next(iter(residues))].add(dof_key)
        else:
            dof_steps[dof_key].completion_classification = "nonlocal"

    groups: List[ResidueLocalCompletionGroup] = []
    atom_to_group: Dict[AtomID, ResidueID] = {}

    for residue_id in sorted(local_dofs):
        group_dofs = set(local_dofs[residue_id])
        group_atoms: Set[AtomID] = set()
        for dof_key in group_dofs:
            group_atoms.update(affected_by_dof[dof_key])

        anchors: Set[AtomID] = set()
        external_pending: Set[AtomID] = set()
        for atom_key in group_atoms:
            for neighbor in covalent_graph.neighbors(atom_key):
                if neighbor in group_atoms:
                    continue
                dependencies = atom_dependencies.get(neighbor)
                if not dependencies:
                    anchors.add(neighbor)
                else:
                    external_pending.add(neighbor)

        external_dof_dependencies: Set[DOFKey] = set()
        for dof_key in group_dofs:
            external_dof_dependencies.update(
                set(dof_steps[dof_key].prerequisite_dofs) - group_dofs
            )
        for atom_key in group_atoms:
            external_dof_dependencies.update(
                set(atom_dependencies.get(atom_key, frozenset())) - group_dofs
            )

        non_bridge_dofs: Set[DOFKey] = set()
        if external_pending or external_dof_dependencies:
            classification = "nonlocal"
        elif not anchors:
            classification = "disconnected"
        elif len(anchors) > 1:
            classification = "residue_local_multi_anchor"
        else:
            # Only a one-anchor, otherwise self-contained candidate needs the
            # more expensive bridge proof. The search expands lazily from the
            # branch-facing second central atom and normally visits only the
            # small distal component.
            for dof_key in group_dofs:
                central = dof_steps[dof_key].central_bond
                if not _is_covalent_bridge(
                    covalent_graph,
                    central[0].atom_id(),
                    central[1].atom_id(),
                ):
                    non_bridge_dofs.add(dof_key)
            classification = (
                "residue_local_cyclic_dof"
                if non_bridge_dofs
                else "residue_local_open_branch"
            )

        group = ResidueLocalCompletionGroup(
            residue_id=residue_id,
            classification=classification,
            dof_keys=_sorted_dofs(group_dofs),
            atom_keys=tuple(
                sorted(group_atoms, key=lambda atom: atom_sort_key(mol, atom))
            ),
            anchor_atoms=tuple(
                sorted(anchors, key=lambda atom: atom_sort_key(mol, atom))
            ),
            non_bridge_dofs=_sorted_dofs(non_bridge_dofs),
            external_pending_atoms=tuple(
                sorted(
                    external_pending,
                    key=lambda atom: atom_sort_key(mol, atom),
                )
            ),
        )
        groups.append(group)
        for dof_key in group_dofs:
            dof_step = dof_steps[dof_key]
            dof_step.local_completion_group = residue_id
            dof_step.completion_classification = classification
        for atom_key in group_atoms:
            atom_to_group[atom_key] = residue_id

    for atom_key, residue_id in atom_to_group.items():
        atom_steps[atom_key].local_completion_group = residue_id

    return tuple(groups)


def _reorder_local_completion_blocks(
    steps: Sequence[PlanStep],
    groups: Sequence[ResidueLocalCompletionGroup],
) -> List[PlanStep]:
    """Move safe residue-local branches directly behind the initial prefix."""
    first_dof_index = next(
        (
            index
            for index, step in enumerate(steps)
            if isinstance(step, PlannedMissingDOFStep)
        ),
        len(steps),
    )
    if first_dof_index == len(steps):
        return list(steps)

    prefix_indices = set(range(first_dof_index))
    indexed_steps = list(enumerate(steps))
    blocks: List[Tuple[int, ResidueID, List[Tuple[int, PlanStep]]]] = []
    selected_indices = set(prefix_indices)

    for group in groups:
        if group.classification != "residue_local_open_branch":
            continue
        dof_keys = set(group.dof_keys)
        atom_keys = set(group.atom_keys)
        members = [
            (index, step)
            for index, step in indexed_steps[first_dof_index:]
            if (
                isinstance(step, PlannedMissingDOFStep)
                and step.dof_key in dof_keys
            ) or (
                isinstance(step, PlannedAtomBuildStep)
                and step.atom_key in atom_keys
            )
        ]
        if not members:
            continue
        blocks.append((members[0][0], group.residue_id, members))
        selected_indices.update(index for index, _step in members)

    blocks.sort(key=lambda item: (item[0], item[1]))
    reordered = list(steps[:first_dof_index])
    for _first_index, _residue_id, members in blocks:
        reordered.extend(step for _index, step in members)
    reordered.extend(
        step
        for index, step in indexed_steps[first_dof_index:]
        if index not in selected_indices
    )
    return reordered


def _validate_planned_dependencies(steps: Sequence[PlanStep]) -> None:
    available_dofs: Set[DOFKey] = set()
    for step in steps:
        if isinstance(step, PlannedMissingDOFStep):
            missing = set(step.prerequisite_dofs) - available_dofs
            if missing:
                raise PlanningError(
                    f"Reordered missing DOF {step.dof_key} precedes prerequisites: "
                    f"{sorted(missing, key=_dof_sort_key)}"
                )
            available_dofs.add(step.dof_key)
        else:
            missing = set(step.required_dofs) - available_dofs
            if missing:
                raise PlanningError(
                    f"Reordered build step {step.atom_key} precedes required DOFs: "
                    f"{sorted(missing, key=_dof_sort_key)}"
                )


def analyze_and_reorder_missing_dof_branches(
    mol: Molecule,
    template: Mapping[str, Any],
    steps: Sequence[PlanStep],
) -> Tuple[List[PlanStep], Tuple[ResidueLocalCompletionGroup, ...]]:
    """Annotate missing-DOF provenance, classify branches, and stably reorder."""
    if not any(isinstance(step, PlannedMissingDOFStep) for step in steps):
        return list(steps), ()
    (
        atom_dependencies,
        affected_by_dof,
        atom_steps,
        dof_steps,
    ) = _annotate_dof_dependencies(mol, template, steps)
    covalent_graph = TemplateCovalentGraphView(mol, template)
    groups = _classify_local_completion_groups(
        mol,
        atom_dependencies,
        affected_by_dof,
        atom_steps,
        dof_steps,
        covalent_graph,
    )
    reordered = _reorder_local_completion_blocks(steps, groups)
    _validate_planned_dependencies(reordered)
    return reordered, groups


# -----------------------------------------------------------------------------
# Planner
# -----------------------------------------------------------------------------

def plan_build_steps(
    mol: Molecule,
    template: Mapping[str, Any],
    *,
    stats: Optional[PlannerStats] = None,
) -> BuildPlan:
    planning_started = time.perf_counter()
    pending_atoms = collect_pending_atoms(mol)
    available_atoms = collect_available_atoms(mol)
    available_dofs: Set[DOFKey] = set()
    steps: List[PlanStep] = []
    dependency_index = build_planner_dependency_index(
        mol,
        template,
        pending_atoms,
    )
    planner_stats = stats if stats is not None else PlannerStats()
    planner_stats.initial_pending_atoms = len(pending_atoms)
    planner_stats.dependency_edges = sum(
        len(rule_ids)
        for rule_ids in dependency_index.atom_to_rules.values()
    )

    availability: Dict[AtomID, AtomAvailability] = {}
    full_versions: Dict[AtomID, int] = defaultdict(int)
    rigid_versions: Dict[AtomID, int] = defaultdict(int)
    current_rigid: Dict[AtomID, Optional[BuildableCandidate]] = {}
    rigid_heap: List[Tuple[Any, ...]] = []
    derived_heap: List[Tuple[Any, ...]] = []
    free_heap: List[Tuple[Any, ...]] = []
    dof_heap: List[Tuple[Any, ...]] = []
    heap_serial = 0

    def next_serial() -> int:
        nonlocal heap_serial
        heap_serial += 1
        return heap_serial

    def set_rigid_candidate(
        target: AtomID,
        candidate: Optional[BuildableCandidate],
    ) -> None:
        rigid_versions[target] += 1
        current_rigid[target] = candidate
        if candidate is not None:
            heapq.heappush(
                rigid_heap,
                (
                    candidate.rule_index,
                    atom_sort_key(mol, target),
                    next_serial(),
                    target,
                    rigid_versions[target],
                    candidate,
                ),
            )

    def push_full_availability(
        target: AtomID,
        evaluation: AtomAvailability,
    ) -> None:
        version = full_versions[target]
        candidate = evaluation.buildable
        if candidate is not None:
            if candidate.torsion_class == "rigid":
                set_rigid_candidate(target, candidate)
            else:
                set_rigid_candidate(target, None)
                heap = (
                    derived_heap
                    if candidate.torsion_class == "derived_rotatable"
                    else free_heap
                )
                heapq.heappush(
                    heap,
                    (
                        candidate.rule_index,
                        atom_sort_key(mol, target),
                        next_serial(),
                        target,
                        version,
                        candidate,
                    ),
                )
        else:
            set_rigid_candidate(target, None)

        missing_dof = evaluation.missing_dof
        if missing_dof is not None:
            heapq.heappush(
                dof_heap,
                (
                    missing_dof.requested_member_index,
                    missing_dof.rule_index,
                    atom_sort_key(mol, target),
                    next_serial(),
                    target,
                    version,
                    missing_dof,
                ),
            )

    def refresh_full(targets: Iterable[AtomID]) -> None:
        for target in sorted(set(targets), key=lambda atom: atom_sort_key(mol, atom)):
            if target not in pending_atoms:
                continue
            evaluation = evaluate_atom(
                mol,
                template,
                target,
                available_atoms,
                available_dofs,
            )
            planner_stats.full_atom_evaluations += 1
            full_versions[target] += 1
            availability[target] = evaluation
            push_full_availability(target, evaluation)

    def refresh_rigid(targets: Iterable[AtomID]) -> None:
        for target in sorted(set(targets), key=lambda atom: atom_sort_key(mol, atom)):
            if target not in pending_atoms:
                continue
            residue = get_residue(mol, target)
            rt = residue_template(template, residue)
            candidate: Optional[BuildableCandidate] = None
            for rule_id in dependency_index.target_to_rigid_rules.get(target, ()):
                rule = rt.get("build_rules", {})[target.atom_name][rule_id.rule_index]
                status, buildable, _missing_dof, _missing_refs = evaluate_rule(
                    mol,
                    template,
                    target,
                    rule,
                    rule_id.rule_index,
                    available_atoms,
                    available_dofs,
                )
                planner_stats.rigid_rule_evaluations += 1
                if status == "buildable" and buildable is not None:
                    candidate = buildable
                    break
            set_rigid_candidate(target, candidate)

    def pop_rigid() -> Optional[BuildableCandidate]:
        while rigid_heap:
            _rule_index, _order, _serial, target, version, candidate = (
                heapq.heappop(rigid_heap)
            )
            if target not in pending_atoms:
                continue
            if rigid_versions[target] != version:
                continue
            if current_rigid.get(target) != candidate:
                continue
            return candidate
        return None

    def pop_buildable(
        heap: List[Tuple[Any, ...]],
        expected_class: str,
    ) -> Optional[BuildableCandidate]:
        while heap:
            _rule_index, _order, _serial, target, version, candidate = (
                heapq.heappop(heap)
            )
            if target not in pending_atoms:
                continue
            if full_versions[target] != version:
                continue
            evaluation = availability.get(target)
            if evaluation is None or evaluation.buildable != candidate:
                continue
            if candidate.torsion_class != expected_class:
                continue
            return candidate
        return None

    def pop_missing_dof() -> Optional[MissingDOFCandidate]:
        while dof_heap:
            (
                _member_index,
                _rule_index,
                _order,
                _serial,
                target,
                version,
                candidate,
            ) = heapq.heappop(dof_heap)
            if target not in pending_atoms:
                continue
            if full_versions[target] != version:
                continue
            evaluation = availability.get(target)
            if evaluation is None or evaluation.missing_dof != candidate:
                continue
            return candidate
        return None

    def append_atom_step(candidate: BuildableCandidate) -> AtomID:
        step = PlannedAtomBuildStep(
            atom_key=candidate.atom_key,
            rule_index=candidate.rule_index,
            torsion_class=candidate.torsion_class,
            torsion_source=candidate.torsion_source,
            observed_member_index=candidate.observed_member_index,
            dof_key=candidate.dof_key,
            torsion_group_index=candidate.torsion_group_index,
            target_member_index=candidate.target_member_index,
        )
        steps.append(step)
        pending_atoms.remove(candidate.atom_key)
        available_atoms.add(candidate.atom_key)
        return candidate.atom_key

    refresh_full(pending_atoms)
    dirty_atoms: Set[AtomID] = set()

    while pending_atoms:
        # Rigid closure: update only rigid rules while collecting every affected
        # target for one complete refresh after the closure.
        while True:
            rigid = pop_rigid()
            if rigid is None:
                if dirty_atoms:
                    refresh_targets = set(dirty_atoms)
                    dirty_atoms.clear()
                    refresh_full(refresh_targets)
                    # The complete refresh can expose a higher-priority rigid
                    # rule, so continue until a true rigid fixed point.
                    continue
                break

            new_atom = append_atom_step(rigid)
            affected = dependency_index.affected_targets(new_atom)
            rigid_affected = dependency_index.rigid_affected_targets(new_atom)
            planner_stats.dependency_target_updates += len(affected)
            dirty_atoms.update(affected)
            refresh_rigid(rigid_affected)

        if not pending_atoms:
            break

        # One derived step can unlock a complete rigid branch.
        derived = pop_buildable(derived_heap, "derived_rotatable")
        if derived is not None:
            new_atom = append_atom_step(derived)
            affected = dependency_index.affected_targets(new_atom)
            planner_stats.dependency_target_updates += len(affected)
            dirty_atoms.update(affected)
            refresh_rigid(
                dependency_index.rigid_affected_targets(new_atom)
            )
            continue

        # Free-rotor closure: plan every currently available group without
        # intermediate rule refreshes, then update the affected union once.
        free = pop_buildable(free_heap, "free_rotor_hydrogen")
        if free is not None:
            while free is not None:
                group_steps = plan_free_rotor_group_steps(
                    mol,
                    template,
                    free,
                    pending_atoms,
                    available_atoms,
                    available_dofs,
                )
                for step in group_steps:
                    if step.atom_key not in pending_atoms:
                        continue
                    steps.append(step)
                    pending_atoms.remove(step.atom_key)
                    available_atoms.add(step.atom_key)
                    affected = dependency_index.affected_targets(step.atom_key)
                    planner_stats.dependency_target_updates += len(affected)
                    dirty_atoms.update(affected)
                free = pop_buildable(
                    free_heap,
                    "free_rotor_hydrogen",
                )
            if dirty_atoms:
                refresh_targets = set(dirty_atoms)
                dirty_atoms.clear()
                refresh_full(refresh_targets)
            continue

        # No buildable atom remains before the next missing-DOF boundary.
        missing_dof = pop_missing_dof()
        if missing_dof is not None:
            dof_step = PlannedMissingDOFStep(
                dof_key=missing_dof.dof_key,
                central_bond=missing_dof.central_bond,
                requested_dihedral_atoms=missing_dof.requested_dihedral_atoms,
                torsion_group_index=missing_dof.torsion_group_index,
                requested_member_index=missing_dof.requested_member_index,
                reason_atom=missing_dof.atom_key,
                reason_rule_index=missing_dof.rule_index,
            )
            steps.append(dof_step)
            available_dofs.add(missing_dof.dof_key)
            rule_id = dependency_index.dof_to_rule.get(
                missing_dof.dof_key
            )
            refresh_full(
                {
                    rule_id.target_atom
                    if rule_id is not None
                    else missing_dof.atom_key
                }
            )
            continue

        # Remaining atoms are not reachable by available rules or a missing DOF.
        examples = []
        for atom_key in sorted(
            pending_atoms,
            key=lambda atom: atom_sort_key(mol, atom),
        )[:20]:
            evaluation = availability.get(atom_key)
            if evaluation is None:
                examples.append(str(atom_key))
            elif evaluation.no_rules:
                examples.append(f"{atom_key}: no_rules")
            elif evaluation.missing_refs:
                refs = ", ".join(
                    f"{ref.chain_id}:{ref.residue_index}:{ref.atom_name}"
                    for ref in evaluation.missing_refs[:5]
                )
                suffix = "..." if len(evaluation.missing_refs) > 5 else ""
                examples.append(
                    f"{atom_key}: missing_refs [{refs}{suffix}]"
                )
            else:
                examples.append(f"{atom_key}: unresolved")
        raise PlanningError(
            "No buildable atom and no missing-DOF boundary remain. "
            "This usually means disconnected missing structure, unsupported "
            "template coverage, or bad input. "
            f"Pending atoms: {len(pending_atoms)}. Examples: "
            + "; ".join(examples)
        )

    raw_planning_finished = time.perf_counter()
    reordered_steps, local_completions = analyze_and_reorder_missing_dof_branches(
        mol,
        template,
        steps,
    )
    planning_finished = time.perf_counter()
    planner_stats.raw_planning_seconds = raw_planning_finished - planning_started
    planner_stats.local_completion_analysis_seconds = (
        planning_finished - raw_planning_finished
    )
    planner_stats.total_planning_seconds = planning_finished - planning_started
    reordered_requirements = [
        step
        for step in reordered_steps
        if isinstance(step, PlannedMissingDOFStep)
    ]
    return BuildPlan(
        steps=reordered_steps,
        unresolved_atoms=set(),
        requirements=reordered_requirements,
        local_completions=local_completions,
    )



# -----------------------------------------------------------------------------
# Coordinate execution helpers
# -----------------------------------------------------------------------------

import copy
import math
import numpy as np


def _coord_array(coord: Any) -> np.ndarray:
    return np.asarray(coord, dtype=float)


def _coord_tuple(coord: np.ndarray) -> Tuple[float, float, float]:
    return (float(coord[0]), float(coord[1]), float(coord[2]))


def _vector_norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def _unit(v: np.ndarray, context: str = "vector") -> np.ndarray:
    n = _vector_norm(v)
    if n < 1.0e-12:
        raise RuntimeError(f"Cannot normalize near-zero {context}")
    return v / n


def compute_dihedral_deg(p0: Any, p1: Any, p2: Any, p3: Any) -> float:
    """Return the signed dihedral angle p0-p1-p2-p3 in degrees."""
    p0 = _coord_array(p0)
    p1 = _coord_array(p1)
    p2 = _coord_array(p2)
    p3 = _coord_array(p3)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = _unit(b1, "central dihedral bond")
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return normalize_angle(math.degrees(math.atan2(y, x)))


def place_atom_from_internal(
    ref0_coord: Any,
    ref1_coord: Any,
    ref2_coord: Any,
    r: float,
    angle_deg: float,
    dihedral_deg: float,
) -> Tuple[float, float, float]:
    """Place target atom from refs ref0-ref1-ref2-target.

    r is distance ref2-target, angle is ref1-ref2-target, and dihedral is
    ref0-ref1-ref2-target, all matching the template internal definition.
    """
    a = _coord_array(ref0_coord)
    b = _coord_array(ref1_coord)
    c = _coord_array(ref2_coord)
    theta = math.radians(float(angle_deg))
    phi = math.radians(float(dihedral_deg))

    bc = _unit(c - b, "ref1-ref2 vector")
    normal = _unit(np.cross(b - a, c - b), "reference plane normal")
    m = np.cross(normal, bc)
    d = c + float(r) * (
        -math.cos(theta) * bc
        + math.sin(theta) * (math.cos(phi) * m + math.sin(phi) * normal)
    )
    return _coord_tuple(d)


def _atom_from_id(mol: Molecule, atom_id: AtomID):
    return get_residue(mol, atom_id).atoms[atom_id.atom_name]


def _coord_from_resolved_ref(mol: Molecule, ref: ResolvedRef):
    atom = _atom_from_id(mol, ref.atom_id())
    if atom.coord is None:
        raise RuntimeError(f"Reference atom has no coordinates: {ref}")
    return atom.coord


def _rule_for_step(mol: Molecule, template: Mapping[str, Any], step: PlannedAtomBuildStep) -> Mapping[str, Any]:
    residue = get_residue(mol, step.atom_key)
    rt = residue_template(template, residue)
    rules = rt.get("build_rules", {}).get(step.atom_key.atom_name, [])
    try:
        return rules[step.rule_index]
    except IndexError as exc:
        raise RuntimeError(f"Build step rule index out of range for {step.atom_key}: {step.rule_index}") from exc


def _resolved_rule_refs(mol: Molecule, residue: Residue, rule: Mapping[str, Any], context: str) -> Tuple[ResolvedRef, ResolvedRef, ResolvedRef]:
    refs = [require_resolved_ref(mol, residue, r, context) for r in rule.get("refs", [])]
    if len(refs) != 3:
        raise RuntimeError(f"{context}: rule does not have exactly three refs")
    return refs[0], refs[1], refs[2]


def _internal_values(rule: Mapping[str, Any]) -> Tuple[float, float, float]:
    internal = rule.get("internal", {})
    try:
        return float(internal["r"]), float(internal["angle"]), float(internal["dihedral"])
    except KeyError as exc:
        raise RuntimeError(f"Rule internal coordinates missing key {exc}") from exc



# -----------------------------------------------------------------------------
# Minimal MM parameter layer and free-rotor scoring
# -----------------------------------------------------------------------------

@dataclass
class FreeRotorMMCache:
    """Topology data shared by all free-rotor searches in one build execution."""

    bond_graph: Dict[AtomID, Set[AtomID]]
    topological_distances: Dict[Tuple[AtomID, AtomID, int], Optional[int]] = field(
        default_factory=dict
    )
    free_rotor_hydrogens: Optional[Set[AtomID]] = None
    collapsed_charge_by_parent: Optional[Dict[AtomID, float]] = None
    # ``None`` retains the historical molecule-wide symmetric behavior
    # controlled by FreeRotorSearchSettings.collapse_other_free_rotors.  The
    # side-chain executors provide an explicit phase-local set; already fixed
    # hydrogens outside those branches then remain explicit MM environment.
    free_rotor_hydrogens_to_collapse: Optional[Set[AtomID]] = None
    spatial_cells: Optional[Dict[Tuple[int, int, int], List["_EnvAtom"]]] = None
    spatial_cell_size_nm: Optional[float] = None


@dataclass(frozen=True)
class _EnvAtom:
    atom_id: AtomID
    coord_nm: np.ndarray
    params: MMAtomParams
    topological_distance: Optional[int]


@dataclass(frozen=True)
class _TrialGroupAtom:
    atom_id: AtomID
    coord_angstrom: Tuple[float, float, float]
    coord_nm: np.ndarray
    params: MMAtomParams
    dihedral_deg: float


@dataclass(frozen=True)
class _VectorizedNonbondedPairs:
    """Contiguous MM arrays for all nonbonded partners of one trial H."""

    coords_nm: np.ndarray
    charges: np.ndarray
    sigmas: np.ndarray
    epsilons: np.ndarray
    scale_lj: np.ndarray
    scale_qq: np.ndarray


@dataclass(frozen=True)
class _FreeRotorPairLists:
    """Preclassified and vectorized partners for one free-rotor H group."""

    pairs_by_hydrogen: Mapping[AtomID, _VectorizedNonbondedPairs]


def _atom_coord_nm(atom: Any) -> np.ndarray:
    return 0.1 * _coord_array(atom.coord)


def _resolve_rtp_atom_token(mol: Molecule, residue: Residue, token: str) -> Optional[AtomID]:
    offset = 0
    atom_name = token
    if token.startswith("-"):
        offset = -1
        atom_name = token[1:]
    elif token.startswith("+"):
        offset = 1
        atom_name = token[1:]
    chain = mol.chains[residue.chain_id]
    ridx = residue.index_in_chain + offset
    if ridx < 0 or ridx >= len(chain.residues):
        return None
    target_res = chain.residues[ridx]
    if atom_name not in target_res.atoms:
        return None
    return AtomID(residue.chain_id, ridx, atom_name)


def build_mm_bond_graph(mol: Molecule, mm_params: MMParameterProvider) -> Dict[AtomID, Set[AtomID]]:
    graph: Dict[AtomID, Set[AtomID]] = {}
    for cid, chain in mol.chains.items():
        for residue in chain.residues:
            for atom_name in residue.atoms:
                graph.setdefault(AtomID(cid, residue.index_in_chain, atom_name), set())
            for a_tok, b_tok in mm_params.residue_bonds.get(residue.ff_resname, []):
                a = _resolve_rtp_atom_token(mol, residue, a_tok)
                b = _resolve_rtp_atom_token(mol, residue, b_tok)
                if a is None or b is None:
                    continue
                graph.setdefault(a, set()).add(b)
                graph.setdefault(b, set()).add(a)
    return graph


def _mm_library_neighbors(
    mol: Molecule,
    mm_params: MMParameterProvider,
    atom_id: AtomID,
) -> Set[AtomID]:
    """Resolve force-field neighbours of one atom without scanning the molecule."""
    residue = get_residue(mol, atom_id)
    neighbors: Set[AtomID] = set()
    for a_tok, b_tok in mm_params.residue_bonds.get(residue.ff_resname, []):
        a = _resolve_rtp_atom_token(mol, residue, a_tok)
        b = _resolve_rtp_atom_token(mol, residue, b_tok)
        if a == atom_id and b is not None:
            neighbors.add(b)
        elif b == atom_id and a is not None:
            neighbors.add(a)
    return neighbors


def build_local_mm_bond_graph(
    mol: Molecule,
    mm_params: MMParameterProvider,
    seed_atom_ids: Iterable[AtomID],
    *,
    max_depth: int = 3,
) -> Dict[AtomID, Set[AtomID]]:
    """Build only the topology needed for local exclusions and proper torsions.

    Nonbonded classification needs at most the 1-2/1-3/1-4 shells around a
    moving atom.  A residue-local side-chain search therefore does not need a
    graph node for every atom in a ribosome-sized molecule.
    """
    if max_depth < 0:
        raise ValueError("max_depth cannot be negative")
    graph: Dict[AtomID, Set[AtomID]] = {}
    best_depth: Dict[AtomID, int] = {}
    queue: deque[Tuple[AtomID, int]] = deque()
    for atom_id in set(seed_atom_ids):
        try:
            _atom_from_id(mol, atom_id)
        except (KeyError, IndexError):
            continue
        best_depth[atom_id] = 0
        queue.append((atom_id, 0))
        graph.setdefault(atom_id, set())

    while queue:
        atom_id, depth = queue.popleft()
        for neighbor in _mm_library_neighbors(mol, mm_params, atom_id):
            graph.setdefault(atom_id, set()).add(neighbor)
            graph.setdefault(neighbor, set()).add(atom_id)
            next_depth = depth + 1
            if next_depth > max_depth:
                continue
            previous_depth = best_depth.get(neighbor)
            if previous_depth is None or next_depth < previous_depth:
                best_depth[neighbor] = next_depth
                queue.append((neighbor, next_depth))
    return graph


def _free_rotor_group_atom_ids(mol: Molecule, anchor_step: PlannedAtomBuildStep) -> List[AtomID]:
    residue = get_residue(mol, anchor_step.atom_key)
    ids: List[AtomID] = []
    for atom_name in anchor_step.free_rotor_group_atoms or (anchor_step.atom_key.atom_name,):
        if atom_name in residue.atoms:
            ids.append(AtomID(anchor_step.atom_key.chain_id, anchor_step.atom_key.residue_index, atom_name))
    if anchor_step.atom_key not in ids:
        ids.insert(0, anchor_step.atom_key)
    return ids


def _free_rotor_hydrogen_ids_from_steps(
    steps: Iterable[PlanStep],
) -> Set[AtomID]:
    """Return only free-rotor H atoms actively built by these plan steps."""
    return {
        step.atom_key
        for step in steps
        if isinstance(step, PlannedAtomBuildStep)
        and step.torsion_class == "free_rotor_hydrogen"
    }


def _heavy_then_free_rotor_steps(
    steps: Iterable[PlannedAtomBuildStep],
) -> Tuple[PlannedAtomBuildStep, ...]:
    """Defer terminal free-rotor H blocks until all heavy geometry is ready."""
    ordered = tuple(steps)
    return tuple(
        step for step in ordered if step.torsion_class != "free_rotor_hydrogen"
    ) + tuple(
        step for step in ordered if step.torsion_class == "free_rotor_hydrogen"
    )


def _free_rotor_phase_offset_for_atom(anchor_step: PlannedAtomBuildStep, atom_id: AtomID) -> float:
    group = tuple(anchor_step.free_rotor_group_atoms or (anchor_step.atom_key.atom_name,))
    if atom_id.atom_name not in group or anchor_step.atom_key.atom_name not in group:
        return 0.0
    n = len(group)
    if n <= 1:
        return 0.0
    anchor_i = group.index(anchor_step.atom_key.atom_name)
    atom_i = group.index(atom_id.atom_name)
    return normalize_angle((atom_i - anchor_i) * 360.0 / n)


def _trial_group_coords(
    mol: Molecule,
    template: Mapping[str, Any],
    anchor_step: PlannedAtomBuildStep,
    anchor_rule: Mapping[str, Any],
    anchor_dihedral: float,
    mm_params: MMParameterProvider,
) -> List[_TrialGroupAtom]:
    trial_atoms: List[_TrialGroupAtom] = []
    for atom_id in _free_rotor_group_atom_ids(mol, anchor_step):
        atom = _atom_from_id(mol, atom_id)
        if atom.coord is not None and atom_id != anchor_step.atom_key:
            # Existing group atoms are not moved by the builder. They are part of
            # the environment already, not trial atoms.
            continue
        res = get_residue(mol, atom_id)
        if atom_id == anchor_step.atom_key:
            rule = anchor_rule
        else:
            rt = residue_template(template, res)
            rules = rt.get("build_rules", {}).get(atom_id.atom_name, [])
            selected_rule = None
            for r in rules:
                td = torsion_data(r)
                if td.get("class") == "free_rotor_hydrogen" and tuple(td.get("group_atoms", ())) == tuple(anchor_step.free_rotor_group_atoms):
                    selected_rule = r
                    break
            if selected_rule is None:
                # Fall back to first free-rotor rule in the same atom if group metadata differs.
                for r in rules:
                    if torsion_class(r) == "free_rotor_hydrogen":
                        selected_rule = r
                        break
            if selected_rule is None:
                continue
            rule = selected_rule
        ref0, ref1, ref2 = _resolved_rule_refs(mol, res, rule, f"free-rotor trial refs for {atom_id}")
        r, angle, _template_dihedral = _internal_values(rule)
        dihedral = normalize_angle(anchor_dihedral + _free_rotor_phase_offset_for_atom(anchor_step, atom_id))
        coord = place_atom_from_internal(
            _coord_from_resolved_ref(mol, ref0),
            _coord_from_resolved_ref(mol, ref1),
            _coord_from_resolved_ref(mol, ref2),
            r,
            angle,
            dihedral,
        )
        params = mm_params.atom_params(res, atom_id.atom_name)
        trial_atoms.append(_TrialGroupAtom(atom_id, coord, 0.1 * _coord_array(coord), params, dihedral))
    return trial_atoms


def _parent_heavy_atom_id(mol: Molecule, step: PlannedAtomBuildStep, rule: Mapping[str, Any]) -> AtomID:
    residue = get_residue(mol, step.atom_key)
    _ref0, _ref1, ref2 = _resolved_rule_refs(mol, residue, rule, f"free-rotor parent refs for {step.atom_key}")
    return ref2.atom_id()


def _template_free_rotor_hydrogen_ids(
    mol: Molecule,
    template: Mapping[str, Any],
) -> Set[AtomID]:
    """Return every atom classified as a member of a free-rotor H group."""
    atom_ids: Set[AtomID] = set()
    for cid, chain in mol.chains.items():
        for residue in chain.residues:
            if (
                not residue.group
                or residue.group not in template
                or residue.ff_resname not in template[residue.group]
            ):
                continue
            rt = residue_template(template, residue)
            build_rules = rt.get("build_rules", {})
            for atom_name, rules in build_rules.items():
                if atom_name not in residue.atoms:
                    continue
                for rule in rules:
                    if torsion_class(rule) != "free_rotor_hydrogen":
                        continue
                    group_names = _free_rotor_group_atoms(rule) or (atom_name,)
                    for group_name in group_names:
                        if group_name in residue.atoms:
                            atom_ids.add(AtomID(cid, residue.index_in_chain, group_name))
                    break
    return atom_ids


def _collapsed_free_rotor_charges(
    mol: Molecule,
    template: Mapping[str, Any],
    mm_params: MMParameterProvider,
    graph: Mapping[AtomID, Set[AtomID]],
    hydrogen_ids: Optional[Iterable[AtomID]] = None,
) -> Tuple[Set[AtomID], Dict[AtomID, float]]:
    """Collapse selected free-rotor H charges onto their parent atoms.

    With ``hydrogen_ids=None`` the historical molecule-wide set is discovered.
    Build executors pass an explicit phase-local set so fixed hydrogens outside
    the active build remain explicit. Collapsed hydrogens are omitted from the
    spatial environment; their LJ parameters are deliberately discarded.
    """
    collapsed_hydrogens = (
        _template_free_rotor_hydrogen_ids(mol, template)
        if hydrogen_ids is None
        else set(hydrogen_ids)
    )
    charge_increments: Dict[AtomID, float] = {}
    for hydrogen_id in collapsed_hydrogens:
        residue = get_residue(mol, hydrogen_id)
        params = mm_params.atom_params(residue, hydrogen_id.atom_name)
        parents = [
            atom_id
            for atom_id in graph.get(hydrogen_id, set())
            if not atom_id.atom_name.upper().startswith("H")
        ]
        if len(parents) != 1:
            raise RuntimeError(
                f"Free-rotor hydrogen {hydrogen_id} must have exactly one parent "
                f"heavy atom in the MM graph; found {len(parents)}"
            )
        parent_id = parents[0]
        charge_increments[parent_id] = charge_increments.get(parent_id, 0.0) + params.charge
    return collapsed_hydrogens, charge_increments


def _spatial_cell_key(coord_nm: np.ndarray, cell_size_nm: float) -> Tuple[int, int, int]:
    return tuple(int(math.floor(float(value) / cell_size_nm)) for value in coord_nm)  # type: ignore[return-value]


def _initialize_free_rotor_spatial_environment(
    mol: Molecule,
    template: Mapping[str, Any],
    mm_params: MMParameterProvider,
    mm_cache: FreeRotorMMCache,
    settings: FreeRotorSearchSettings,
) -> None:
    """Build the order-independent MM environment and its spatial cell list once."""
    if mm_cache.spatial_cells is not None:
        if mm_cache.spatial_cell_size_nm != settings.preselection_radius_nm:
            raise RuntimeError(
                "Free-rotor spatial cache was built with a different "
                "preselection radius"
            )
        return
    if settings.preselection_radius_nm <= 0.0:
        raise ValueError("preselection_radius_nm must be positive")

    free_hydrogens: Set[AtomID] = set()
    charge_increments: Dict[AtomID, float] = {}
    selected_hydrogens = mm_cache.free_rotor_hydrogens_to_collapse
    if selected_hydrogens is not None:
        if selected_hydrogens:
            free_hydrogens, charge_increments = _collapsed_free_rotor_charges(
                mol,
                template,
                mm_params,
                mm_cache.bond_graph,
                selected_hydrogens,
            )
    elif settings.collapse_other_free_rotors:
        free_hydrogens, charge_increments = _collapsed_free_rotor_charges(
            mol,
            template,
            mm_params,
            mm_cache.bond_graph,
        )

    cells: Dict[Tuple[int, int, int], List[_EnvAtom]] = defaultdict(list)
    for cid, chain in mol.chains.items():
        for residue in chain.residues:
            for atom_name, atom in residue.atoms.items():
                atom_id = AtomID(cid, residue.index_in_chain, atom_name)
                if atom_id in free_hydrogens or atom.coord is None:
                    continue
                try:
                    params = mm_params.atom_params(residue, atom_name)
                except KeyError:
                    continue
                charge_increment = charge_increments.get(atom_id, 0.0)
                if charge_increment:
                    params = MMAtomParams(
                        atom_type=params.atom_type,
                        charge=params.charge + charge_increment,
                        sigma=params.sigma,
                        epsilon=params.epsilon,
                    )
                coord_nm = _atom_coord_nm(atom)
                cells[_spatial_cell_key(coord_nm, settings.preselection_radius_nm)].append(
                    _EnvAtom(atom_id, coord_nm, params, None)
                )

    mm_cache.free_rotor_hydrogens = free_hydrogens
    mm_cache.collapsed_charge_by_parent = charge_increments
    mm_cache.spatial_cells = dict(cells)
    mm_cache.spatial_cell_size_nm = settings.preselection_radius_nm


def _free_rotor_environment_atoms(
    mol: Molecule,
    template: Mapping[str, Any],
    mm_params: MMParameterProvider,
    parent_atom_id: AtomID,
    group_atom_ids: Sequence[AtomID],
    mm_cache: FreeRotorMMCache,
    settings: FreeRotorSearchSettings,
) -> List[_EnvAtom]:
    """Select a spherical local environment from the prebuilt spatial cells."""
    _initialize_free_rotor_spatial_environment(
        mol,
        template,
        mm_params,
        mm_cache,
        settings,
    )
    parent_atom = _atom_from_id(mol, parent_atom_id)
    if parent_atom.coord is None:
        raise RuntimeError(f"Free-rotor parent atom has no coordinates: {parent_atom_id}")
    parent_nm = _atom_coord_nm(parent_atom)
    group_set = set(group_atom_ids)
    cells = mm_cache.spatial_cells or {}
    center = _spatial_cell_key(parent_nm, settings.preselection_radius_nm)
    cutoff_squared = settings.preselection_radius_nm ** 2
    env: List[_EnvAtom] = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for candidate in cells.get(
                    (center[0] + dx, center[1] + dy, center[2] + dz), ()
                ):
                    if candidate.atom_id in group_set:
                        continue
                    delta = candidate.coord_nm - parent_nm
                    if float(np.dot(delta, delta)) <= cutoff_squared:
                        env.append(candidate)
    return env


def _topological_shells_upto_three(
    graph: Mapping[AtomID, Set[AtomID]],
    start: AtomID,
) -> Tuple[Set[AtomID], Set[AtomID], Set[AtomID]]:
    """Backward-compatible typed wrapper around the shared MM helper."""
    shell_12, shell_13, shell_14 = topological_shells_upto_three(graph, start)
    return set(shell_12), set(shell_13), set(shell_14)


def _prepare_free_rotor_pair_lists(
    graph: Mapping[AtomID, Set[AtomID]],
    group_atom_ids: Sequence[AtomID],
    env_atoms: Sequence[_EnvAtom],
    mm_params: MMParameterProvider,
) -> _FreeRotorPairLists:
    """Classify exclusions and 1-4 scaling once, outside the angle scan."""
    pair_arrays: Dict[AtomID, _VectorizedNonbondedPairs] = {}
    for hydrogen_id in group_atom_ids:
        shell_12, shell_13, shell_14 = _topological_shells_upto_three(
            graph, hydrogen_id
        )
        excluded = shell_12 | shell_13 | {hydrogen_id}
        partners = [atom for atom in env_atoms if atom.atom_id not in excluded]
        pair_arrays[hydrogen_id] = _VectorizedNonbondedPairs(
            coords_nm=np.asarray([atom.coord_nm for atom in partners], dtype=float).reshape(-1, 3),
            charges=np.asarray([atom.params.charge for atom in partners], dtype=float),
            sigmas=np.asarray([atom.params.sigma for atom in partners], dtype=float),
            epsilons=np.asarray([atom.params.epsilon for atom in partners], dtype=float),
            scale_lj=np.asarray(
                [
                    mm_params.fudge_lj if atom.atom_id in shell_14 else 1.0
                    for atom in partners
                ],
                dtype=float,
            ),
            scale_qq=np.asarray(
                [
                    mm_params.fudge_qq if atom.atom_id in shell_14 else 1.0
                    for atom in partners
                ],
                dtype=float,
            ),
        )
    return _FreeRotorPairLists(pair_arrays)


def _torsion_atom_ids_with_h_terminal(
    graph: Mapping[AtomID, Set[AtomID]],
    h_atom_id: AtomID,
) -> List[Tuple[AtomID, AtomID, AtomID, AtomID]]:
    out: List[Tuple[AtomID, AtomID, AtomID, AtomID]] = []
    for c in graph.get(h_atom_id, set()):
        for b in graph.get(c, set()):
            if b == h_atom_id:
                continue
            for a in graph.get(b, set()):
                if a == c:
                    continue
                out.append((a, b, c, h_atom_id))
    # deterministic order
    return sorted(set(out), key=lambda t: (t[0], t[1], t[2], t[3]))


def _trial_coord_lookup(trial_atoms: Sequence[_TrialGroupAtom]) -> Dict[AtomID, Tuple[float, float, float]]:
    return {ta.atom_id: ta.coord_angstrom for ta in trial_atoms}


def _coord_for_atom_id_angstrom(mol: Molecule, atom_id: AtomID, trial_lookup: Mapping[AtomID, Tuple[float, float, float]]) -> Tuple[float, float, float]:
    if atom_id in trial_lookup:
        return trial_lookup[atom_id]
    atom = _atom_from_id(mol, atom_id)
    if atom.coord is None:
        raise RuntimeError(f"Atom lacks coordinates for torsion evaluation: {atom_id}")
    return atom.coord


def _torsion_energy_for_trial_group(
    mol: Molecule,
    mm_params: MMParameterProvider,
    graph: Mapping[AtomID, Set[AtomID]],
    trial_atoms: Sequence[_TrialGroupAtom],
) -> float:
    trial_lookup = _trial_coord_lookup(trial_atoms)
    energy = 0.0
    for trial in trial_atoms:
        for a, b, c, d in _torsion_atom_ids_with_h_terminal(graph, trial.atom_id):
            try:
                ra, rb, rc, rd = (get_residue(mol, x) for x in (a, b, c, d))
                atypes = (
                    mm_params.atom_params(ra, a.atom_name).atom_type,
                    mm_params.atom_params(rb, b.atom_name).atom_type,
                    mm_params.atom_params(rc, c.atom_name).atom_type,
                    mm_params.atom_params(rd, d.atom_name).atom_type,
                )
            except KeyError:
                continue
            terms = mm_params.proper_dihedral_terms(atypes)
            if not terms:
                continue
            coords = [_coord_for_atom_id_angstrom(mol, x, trial_lookup) for x in (a, b, c, d)]
            phi = compute_dihedral_deg(*coords)
            energy += proper_torsion_energy(phi, terms)
    return energy


def _nonbonded_energy_for_trial_group(
    mm_params: MMParameterProvider,
    trial_atoms: Sequence[_TrialGroupAtom],
    pair_lists: _FreeRotorPairLists,
    settings: FreeRotorSearchSettings,
) -> float:
    energy = 0.0
    for trial in trial_atoms:
        pairs = pair_lists.pairs_by_hydrogen.get(trial.atom_id)
        if pairs is None:
            continue
        energy += vectorized_nonbonded_energy(
            trial.coord_nm,
            trial.params,
            pairs.coords_nm,
            pairs.charges,
            pairs.sigmas,
            pairs.epsilons,
            pairs.scale_lj,
            pairs.scale_qq,
            switch_radius_nm=settings.switch_radius_nm,
            cutoff_radius_nm=settings.cutoff_radius_nm,
            include_lj=settings.include_lj,
            include_electrostatics=settings.include_electrostatics,
        )
    return energy


def score_free_rotor_group(
    mol: Molecule,
    template: Mapping[str, Any],
    anchor_step: PlannedAtomBuildStep,
    anchor_rule: Mapping[str, Any],
    anchor_dihedral: float,
    mm_params: MMParameterProvider,
    mm_cache: FreeRotorMMCache,
    pair_lists: _FreeRotorPairLists,
    settings: FreeRotorSearchSettings,
) -> float:
    trial_atoms = _trial_group_coords(mol, template, anchor_step, anchor_rule, anchor_dihedral, mm_params)
    energy = 0.0
    if settings.include_torsions:
        energy += _torsion_energy_for_trial_group(
            mol,
            mm_params,
            mm_cache.bond_graph,
            trial_atoms,
        )
    energy += _nonbonded_energy_for_trial_group(
        mm_params,
        trial_atoms,
        pair_lists,
        settings,
    )
    return energy


def choose_free_rotor_dihedral_mm_search(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    rule: Mapping[str, Any],
    mm_params: MMParameterProvider,
    settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
) -> float:
    """Choose the anchor H dihedral by local MM scan of the whole free-rotor H group."""
    settings = settings or FreeRotorSearchSettings()
    _r, _angle, template_dihedral = _internal_values(rule)

    if mm_cache is None:
        mm_cache = FreeRotorMMCache(build_mm_bond_graph(mol, mm_params))
    graph = mm_cache.bond_graph
    parent_id = _parent_heavy_atom_id(mol, step, rule)
    group_ids = _free_rotor_group_atom_ids(mol, step)
    env_atoms = _free_rotor_environment_atoms(
        mol,
        template,
        mm_params,
        parent_id,
        group_ids,
        mm_cache,
        settings,
    )
    pair_lists = _prepare_free_rotor_pair_lists(
        graph, group_ids, env_atoms, mm_params
    )

    group_size = max(
        1,
        len(step.free_rotor_group_atoms or (step.atom_key.atom_name,)),
    )
    return optimize_symmetric_periodic_1d(
        template_dihedral,
        group_size,
        settings.grid_step_deg,
        lambda phi: score_free_rotor_group(
            mol,
            template,
            step,
            rule,
            phi,
            mm_params,
            mm_cache,
            pair_lists,
            settings,
        ),
    )


def choose_free_rotor_dihedral_fallback(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    rule: Mapping[str, Any],
) -> float:
    """Fallback free-rotor policy: use the template/internal dihedral."""
    _r, _angle, dihedral = _internal_values(rule)
    return dihedral


def choose_free_rotor_dihedral(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    rule: Mapping[str, Any],
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
) -> float:
    if mm_params is None:
        return choose_free_rotor_dihedral_fallback(mol, template, step, rule)
    return choose_free_rotor_dihedral_mm_search(
        mol,
        template,
        step,
        rule,
        mm_params,
        free_rotor_settings,
        mm_cache,
    )


def _target_member_phase(group: Mapping[str, Any], member_index: int) -> float:
    return float(group.get("members", [])[member_index].get("phase", 0.0))


def _validated_dof_value(
    dof_key: DOFKey,
    dof_values: Optional[Mapping[DOFKey, float]],
) -> float:
    if dof_values is None or dof_key not in dof_values:
        raise RuntimeError(f"No external value was supplied for missing DOF {dof_key}")
    try:
        value = float(dof_values[dof_key])
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Missing DOF {dof_key} has a non-numeric value: "
            f"{dof_values[dof_key]!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(f"Missing DOF {dof_key} must be finite, got {value!r}")
    return normalize_angle(value)


def _build_dihedral_for_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    rule: Mapping[str, Any],
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
    dof_values: Optional[Mapping[DOFKey, float]] = None,
) -> float:
    residue = get_residue(mol, step.atom_key)

    if step.torsion_source == "internal":
        _r, _angle, dihedral = _internal_values(rule)
        return dihedral

    if step.torsion_source == "free_rotor_search":
        return choose_free_rotor_dihedral(
            mol,
            template,
            step,
            rule,
            mm_params,
            free_rotor_settings,
            mm_cache,
        )

    if step.torsion_source == "observed_member":
        if step.observed_member_index is None:
            raise RuntimeError(f"Observed-member build step lacks observed_member_index: {step}")
        rt = residue_template(template, residue)
        gi, target_mi, group, _rev = find_matching_torsion_group(mol, residue, rt, rule, step.atom_key)
        observed_atoms = member_dihedral_atoms(mol, residue, group, step.observed_member_index)
        if observed_atoms is None:
            raise RuntimeError(
                f"Planned observed torsion member is no longer resolvable: {step}"
            )
        coords = [_coord_from_resolved_ref(mol, r) for r in observed_atoms]
        observed_phi = compute_dihedral_deg(*coords)
        observed_phase = _target_member_phase(group, step.observed_member_index)
        target_phase = _target_member_phase(group, target_mi)
        return normalize_angle(observed_phi + target_phase - observed_phase)

    if step.torsion_source == "free_rotor_group_phase":
        if step.free_rotor_anchor is None or step.phase_offset is None:
            raise RuntimeError(f"Free-rotor group phase step lacks anchor/phase_offset: {step}")
        residue = get_residue(mol, step.atom_key)
        ref0, ref1, ref2 = _resolved_rule_refs(mol, residue, rule, f"free-rotor group phase refs for {step.atom_key}")
        atom0 = _coord_from_resolved_ref(mol, ref0)
        atom1 = _coord_from_resolved_ref(mol, ref1)
        atom2 = _coord_from_resolved_ref(mol, ref2)
        anchor_atom = _atom_from_id(mol, step.free_rotor_anchor)
        if anchor_atom.coord is None:
            raise RuntimeError(f"Free-rotor anchor has no coordinates: {step.free_rotor_anchor}")
        anchor_phi = compute_dihedral_deg(atom0, atom1, atom2, anchor_atom.coord)
        return normalize_angle(anchor_phi + float(step.phase_offset))

    if step.torsion_source == "supplied_dof":
        if step.dof_key is None:
            raise RuntimeError(
                f"Supplied-DOF build step lacks a DOF key: {step.atom_key}"
            )
        return _validated_dof_value(step.dof_key, dof_values)

    raise RuntimeError(f"Unsupported torsion_source {step.torsion_source!r} for {step.atom_key}")


def _assign_atom_from_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    *,
    require_existing: bool,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
    dof_values: Optional[Mapping[DOFKey, float]] = None,
) -> None:
    residue = get_residue(mol, step.atom_key)
    atom = residue.atoms[step.atom_key.atom_name]
    if require_existing and atom.coord is None:
        raise RuntimeError(f"Rebuild target has no coordinates: {step.atom_key}")
    if not require_existing and atom.coord is not None:
        # A plan generated for an older molecule was reused. Treat already-built
        # atoms as an error rather than silently overwriting coordinates.
        raise RuntimeError(f"Target atom already has coordinates: {step.atom_key}")

    rule = _rule_for_step(mol, template, step)
    ref0, ref1, ref2 = _resolved_rule_refs(mol, residue, rule, f"build refs for {step.atom_key}")
    r, angle, _template_dihedral = _internal_values(rule)
    dihedral = _build_dihedral_for_step(
        mol,
        template,
        step,
        rule,
        mm_params,
        free_rotor_settings,
        mm_cache,
        dof_values,
    )

    coord = place_atom_from_internal(
        _coord_from_resolved_ref(mol, ref0),
        _coord_from_resolved_ref(mol, ref1),
        _coord_from_resolved_ref(mol, ref2),
        r,
        angle,
        dihedral,
    )
    atom.coord = coord
    atom.built = True
    atom.build_source = step.torsion_source
    atom.build_rule_index = step.rule_index
    atom.occupancy = 1.0
    atom.bfactor = 0.0


def build_atom_from_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
    dof_values: Optional[Mapping[DOFKey, float]] = None,
) -> None:
    """Build one currently missing target atom in-place."""
    _assign_atom_from_step(
        mol,
        template,
        step,
        require_existing=False,
        mm_params=mm_params,
        free_rotor_settings=free_rotor_settings,
        mm_cache=mm_cache,
        dof_values=dof_values,
    )


def rebuild_atom_from_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
    dof_values: Optional[Mapping[DOFKey, float]] = None,
) -> None:
    """Explicitly overwrite one already built target atom in-place."""
    _assign_atom_from_step(
        mol,
        template,
        step,
        require_existing=True,
        mm_params=mm_params,
        free_rotor_settings=free_rotor_settings,
        mm_cache=mm_cache,
        dof_values=dof_values,
    )


def _copy_molecule_for_coordinate_build(molecule: Molecule) -> Molecule:
    """Copy the mutable molecule model without generic recursive deepcopy.

    Coordinate building mutates Atom records and, downstream, callers may move
    or otherwise edit the returned molecule. Every mutable molecule/chain/
    residue/atom record is therefore independent, while immutable strings,
    numbers and coordinate tuples can be shared safely. This is materially
    faster than `copy.deepcopy` for ribosome-sized GUI previews.
    """
    result = copy.copy(molecule)
    result.chains = {}
    for chain_id, chain in molecule.chains.items():
        copied_chain = copy.copy(chain)
        copied_chain.residues = []
        for residue in chain.residues:
            copied_residue = copy.copy(residue)
            copied_residue.atoms = {
                name: copy.copy(atom)
                for name, atom in residue.atoms.items()
            }
            copied_residue.observed_extra_atoms = {
                name: copy.copy(atom)
                for name, atom in residue.observed_extra_atoms.items()
            }
            copied_residue.connectivity_parts = [
                list(part) for part in residue.connectivity_parts
            ]
            copied_chain.residues.append(copied_residue)
        result.chains[chain_id] = copied_chain

    result.covalent_bonds = list(molecule.covalent_bonds)
    result.passthrough_atoms = [
        copy.copy(atom) for atom in molecule.passthrough_atoms
    ]
    result.unassigned_records = [
        copy.copy(atom) for atom in molecule.unassigned_records
    ]
    result.warnings = list(molecule.warnings)
    result.periodic_box = copy.copy(molecule.periodic_box)
    return result


def _copy_molecule_for_local_coordinate_build(
    molecule: Molecule,
    mutable_residue_ids: Iterable[ResidueID],
) -> Molecule:
    """Copy only residues whose trial atom coordinates will be mutated.

    Side-chain MM optimization returns DOF values, not a molecule.  All atoms
    outside the selected residue-local branches are therefore read-only and
    can safely be shared with the input molecule.  This avoids copying every
    Atom object in a ribosome merely to rotate a handful of local atoms.
    """
    mutable = set(mutable_residue_ids)
    result = copy.copy(molecule)
    result.chains = {}
    for chain_id, chain in molecule.chains.items():
        copied_chain = copy.copy(chain)
        copied_chain.residues = list(chain.residues)
        for residue_id in mutable:
            if residue_id.chain_id != chain_id:
                continue
            residue = chain.residues[residue_id.residue_index]
            copied_residue = copy.copy(residue)
            copied_residue.atoms = {
                name: copy.copy(atom) for name, atom in residue.atoms.items()
            }
            copied_residue.observed_extra_atoms = {
                name: copy.copy(atom)
                for name, atom in residue.observed_extra_atoms.items()
            }
            copied_residue.connectivity_parts = [
                list(part) for part in residue.connectivity_parts
            ]
            copied_chain.residues[residue_id.residue_index] = copied_residue
        result.chains[chain_id] = copied_chain
    return result


def execute_build_plan_until_missing_dof(
    molecule: Molecule,
    template: Mapping[str, Any],
    build_plan: BuildPlan,
    *,
    dof_values: Optional[Mapping[DOFKey, float]] = None,
    modify_myself: bool = False,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
) -> Tuple[Molecule, BuildPlan]:
    """Execute until the first unresolved missing-DOF step or plan end.

    If modify_myself is False, both molecule and build_plan are deep-copied and
    the returned objects are independent trial outputs. If True, the input
    objects are modified in-place and also returned. Missing-DOF steps present
    in `dof_values` are consumed and their values are used by the corresponding
    `supplied_dof` atom steps. With no supplied values, the historical behavior
    of stopping at the first missing DOF is unchanged.
    """
    mol = molecule if modify_myself else _copy_molecule_for_coordinate_build(molecule)
    plan = build_plan if modify_myself else copy.deepcopy(build_plan)

    effective_free_rotor_settings = (
        free_rotor_settings or FreeRotorSearchSettings()
    )
    mm_cache = (
        FreeRotorMMCache(
            build_mm_bond_graph(mol, mm_params),
            # ``None`` activates the established molecule-wide symmetric
            # monopole environment.  This is intentional for the initial
            # build node so free-rotor H placement is independent of plan
            # order.  ``set()`` explicitly disables it when requested.
            free_rotor_hydrogens_to_collapse=(
                None
                if effective_free_rotor_settings.collapse_other_free_rotors
                else set()
            ),
        )
        if mm_params is not None
        else None
    )
    consumed_count = 0
    for step in plan.steps:
        if isinstance(step, PlannedMissingDOFStep):
            if dof_values is None or step.dof_key not in dof_values:
                break
            # Validate at the boundary so a malformed value cannot consume the
            # request and leave the copied plan in a misleading state.
            _validated_dof_value(step.dof_key, dof_values)
            consumed_count += 1
            continue
        if not isinstance(step, PlannedAtomBuildStep):
            raise RuntimeError(f"Unknown plan step type: {step!r}")
        build_atom_from_step(
            mol,
            template,
            step,
            mm_params,
            effective_free_rotor_settings,
            mm_cache,
            dof_values,
        )
        consumed_count += 1
    if consumed_count:
        plan.steps = plan.steps[consumed_count:]

    # Remaining requirement list should describe only missing-DOF steps still in
    # the unexecuted plan prefix/suffix.
    plan.requirements = [s for s in plan.steps if isinstance(s, PlannedMissingDOFStep)]
    remaining_dofs = {step.dof_key for step in plan.requirements}
    plan.local_completions = tuple(
        group
        for group in plan.local_completions
        if any(dof_key in remaining_dofs for dof_key in group.dof_keys)
    )
    plan.unresolved_atoms = collect_pending_atoms(mol)
    return mol, plan


def _open_branch_groups(
    plan: BuildPlan,
) -> Dict[ResidueID, ResidueLocalCompletionGroup]:
    return {
        group.residue_id: group
        for group in plan.local_completions
        if group.classification == "residue_local_open_branch"
    }


def prepare_sidechain_execution_index(
    remaining_plan: BuildPlan,
) -> SidechainExecutionIndex:
    """Partition safe local branches and index their transitive DOF effects.

    The source plan is never modified. Step objects are intentionally shared as
    immutable execution recipes; only the containing tuples and lookup mapping
    are new. A DOF affects exactly those atom steps whose planner-propagated
    `required_dofs` contains its key.
    """
    groups = _open_branch_groups(remaining_plan)
    safe_group_ids = set(groups)
    sidechain_steps: List[PlanStep] = []
    residual_steps: List[PlanStep] = []
    dof_steps: List[PlannedMissingDOFStep] = []
    atom_steps: List[PlannedAtomBuildStep] = []
    residue_groups: List[ResidueID] = []
    seen_groups: Set[ResidueID] = set()

    for step in remaining_plan.steps:
        residue_id = getattr(step, "local_completion_group", None)
        if residue_id not in safe_group_ids:
            residual_steps.append(step)
            continue
        sidechain_steps.append(step)
        if residue_id not in seen_groups:
            seen_groups.add(residue_id)
            residue_groups.append(residue_id)
        if isinstance(step, PlannedMissingDOFStep):
            if step.dof_key not in groups[residue_id].dof_keys:
                raise PlanningError(
                    f"Safe local group {residue_id} contains an unexpected DOF "
                    f"step {step.dof_key}"
                )
            dof_steps.append(step)
        elif isinstance(step, PlannedAtomBuildStep):
            atom_steps.append(step)
        else:
            raise PlanningError(f"Unknown side-chain plan step: {step!r}")

    remaining_dof_keys = {step.dof_key for step in dof_steps}
    expected_groups = {
        residue_id
        for residue_id, group in groups.items()
        if any(dof_key in remaining_dof_keys for dof_key in group.dof_keys)
    }
    if seen_groups != expected_groups:
        raise PlanningError(
            "Safe local completion groups and remaining plan steps disagree; "
            f"groups_without_steps={sorted(expected_groups - seen_groups)}, "
            f"steps_without_groups={sorted(seen_groups - expected_groups)}"
        )

    affected: Dict[DOFKey, Tuple[PlannedAtomBuildStep, ...]] = {}
    for dof_step in dof_steps:
        affected_steps = tuple(
            step
            for step in atom_steps
            if dof_step.dof_key in step.required_dofs
        )
        if not affected_steps:
            raise PlanningError(
                f"Safe local missing DOF affects no planned atom: {dof_step.dof_key}"
            )
        affected[dof_step.dof_key] = affected_steps

    return SidechainExecutionIndex(
        residue_groups=tuple(residue_groups),
        dof_steps=tuple(dof_steps),
        atom_steps=tuple(atom_steps),
        sidechain_steps=tuple(sidechain_steps),
        residual_steps=tuple(residual_steps),
        affected_atom_steps_by_dof=affected,
    )


def _sidechain_execution_index(
    remaining_plan: BuildPlan,
    execution_index: Optional[SidechainExecutionIndex],
) -> SidechainExecutionIndex:
    return (
        execution_index
        if execution_index is not None
        else prepare_sidechain_execution_index(remaining_plan)
    )


def _open_branch_dof_steps(
    plan: BuildPlan,
    execution_index: Optional[SidechainExecutionIndex] = None,
) -> List[PlannedMissingDOFStep]:
    return list(_sidechain_execution_index(plan, execution_index).dof_steps)


def get_template_sidechain_dof_defaults(
    molecule: Molecule,
    template: Mapping[str, Any],
    remaining_plan: BuildPlan,
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
) -> Dict[DOFKey, float]:
    """Return cheap template fallbacks for safe local completion DOFs.

    Production preview preparation should normally use
    ``optimize_missing_sidechain_dofs``. This fallback deliberately has the
    same result shape and remains useful without force-field parameters.
    """
    defaults: Dict[DOFKey, float] = {}
    for step in _open_branch_dof_steps(remaining_plan, execution_index):
        _residue, _residue_template, rule = _selected_rule(
            molecule,
            template,
            step.reason_atom,
            step.reason_rule_index,
        )
        _r, _angle, dihedral = _internal_values(rule)
        defaults[step.dof_key] = normalize_angle(dihedral)
    return defaults


@dataclass(frozen=True)
class _SidechainBranchModel:
    residue_id: ResidueID
    dof_steps: Tuple[PlannedMissingDOFStep, ...]
    atom_steps: Tuple[PlannedAtomBuildStep, ...]
    energy_atom_ids: Tuple[AtomID, ...]


def _sidechain_branch_models(
    index: SidechainExecutionIndex,
    free_rotor_hydrogens: Set[AtomID],
) -> Tuple[_SidechainBranchModel, ...]:
    models: List[_SidechainBranchModel] = []
    for residue_id in index.residue_groups:
        dof_steps = tuple(
            step
            for step in index.dof_steps
            if step.local_completion_group == residue_id
        )
        atom_steps = tuple(
            step
            for step in index.atom_steps
            if step.local_completion_group == residue_id
        )
        models.append(
            _SidechainBranchModel(
                residue_id=residue_id,
                dof_steps=dof_steps,
                atom_steps=atom_steps,
                energy_atom_ids=tuple(
                    step.atom_key
                    for step in atom_steps
                    if step.atom_key not in free_rotor_hydrogens
                ),
            )
        )
    return tuple(models)


def _sidechain_coordinate_generator(
    geometry_molecule: Molecule,
    template: Mapping[str, Any],
    atom_steps: Sequence[PlannedAtomBuildStep],
    dof_order: Sequence[DOFKey],
):
    """Return a history-independent local coordinate callback with memoization."""
    ordered_dofs = tuple(dof_order)
    cache: Dict[Tuple[float, ...], Dict[AtomID, Tuple[float, float, float]]] = {}

    def generate(
        values: Mapping[DOFKey, float],
    ) -> Mapping[AtomID, Tuple[float, float, float]]:
        key = tuple(round(normalize_angle(values[dof]), 10) for dof in ordered_dofs)
        if key in cache:
            return cache[key]
        for step in atom_steps:
            atom = _atom_from_id(geometry_molecule, step.atom_key)
            atom.coord = None
            atom.built = False
            atom.build_source = None
            atom.build_rule_index = None
        for step in atom_steps:
            # MM free-rotor placement is intentionally not nested inside the
            # heavy-branch search.  Its charge is collapsed onto the parent in
            # this optimization and the H geometry is optimized afterwards.
            build_atom_from_step(
                geometry_molecule,
                template,
                step,
                mm_params=None,
                dof_values=values,
            )
        result = {
            step.atom_key: tuple(_atom_from_id(geometry_molecule, step.atom_key).coord)
            for step in atom_steps
        }
        cache[key] = result
        return result

    return generate


def _coarse_sidechain_conformers(
    dof_steps: Sequence[PlannedMissingDOFStep],
    default_values: Mapping[DOFKey, float],
    coordinates,
    settings: SidechainOptimizationSettings,
) -> List[Tuple[Dict[DOFKey, float], Mapping[AtomID, Tuple[float, float, float]]]]:
    """Sample every DOF independently to obtain a deterministic swept envelope."""
    dof_order = tuple(step.dof_key for step in dof_steps)
    base = {dof: float(default_values[dof]) for dof in dof_order}
    n_grid = max(3, int(round(360.0 / settings.grid_step_deg)))
    step_deg = 360.0 / n_grid
    states: Dict[Tuple[float, ...], Dict[DOFKey, float]] = {}

    def remember(values: Mapping[DOFKey, float]) -> None:
        clean = {dof: normalize_angle(values[dof]) for dof in dof_order}
        key = tuple(round(clean[dof], 10) for dof in dof_order)
        states[key] = clean

    remember(base)
    for dof in dof_order:
        for grid_index in range(n_grid):
            trial = dict(base)
            trial[dof] = normalize_angle(base[dof] + grid_index * step_deg)
            remember(trial)
    return [(state, coordinates(state)) for state in states.values()]


@dataclass(frozen=True)
class _PreparedSidechainMMEnvironment:
    fixed_coords: Mapping[AtomID, Tuple[float, float, float]]
    params_by_atom: Mapping[AtomID, MMAtomParams]
    spatial_index: MMSpatialIndex


def _prepare_sidechain_mm_environment(
    molecule: Molecule,
    mm_params: MMParameterProvider,
    collapsed_hydrogens: Set[AtomID],
    collapsed_charge_by_parent: Mapping[AtomID, float],
    settings: SidechainOptimizationSettings,
) -> _PreparedSidechainMMEnvironment:
    """Create coordinates, effective parameters and the cell list in one pass."""
    fixed_coords: Dict[AtomID, Tuple[float, float, float]] = {}
    params_by_atom: Dict[AtomID, MMAtomParams] = {}
    spatial_atoms: List[SpatialAtom] = []
    for chain_id, chain in molecule.chains.items():
        for residue in chain.residues:
            for atom_name, atom in residue.atoms.items():
                if atom.coord is None:
                    continue
                atom_id = AtomID(chain_id, residue.index_in_chain, atom_name)
                coord = (
                    float(atom.coord[0]),
                    float(atom.coord[1]),
                    float(atom.coord[2]),
                )
                fixed_coords[atom_id] = coord
                if atom_id in collapsed_hydrogens:
                    continue
                try:
                    params = mm_params.atom_params(residue, atom_name)
                except KeyError:
                    continue
                increment = float(collapsed_charge_by_parent.get(atom_id, 0.0))
                if increment:
                    params = MMAtomParams(
                        atom_type=params.atom_type,
                        charge=params.charge + increment,
                        sigma=params.sigma,
                        epsilon=params.epsilon,
                    )
                params_by_atom[atom_id] = params
                spatial_atoms.append(
                    SpatialAtom(atom_id, 0.1 * _coord_array(coord), params)
                )
    return _PreparedSidechainMMEnvironment(
        fixed_coords=fixed_coords,
        params_by_atom=params_by_atom,
        spatial_index=MMSpatialIndex(
            spatial_atoms,
            settings.preselection_radius_nm,
        ),
    )


def _pair_scaling(
    graph: Mapping[AtomID, Set[AtomID]],
    atom_a: AtomID,
    atom_b: AtomID,
    mm_params: MMParameterProvider,
    shell_cache: Dict[AtomID, Tuple[Set[AtomID], Set[AtomID], Set[AtomID]]],
) -> Optional[Tuple[float, float]]:
    if atom_a == atom_b:
        return None
    if atom_a not in shell_cache:
        shell_cache[atom_a] = _topological_shells_upto_three(graph, atom_a)
    shell_12, shell_13, shell_14 = shell_cache[atom_a]
    if atom_b in shell_12 or atom_b in shell_13:
        return None
    if atom_b in shell_14:
        return mm_params.fudge_lj, mm_params.fudge_qq
    return 1.0, 1.0


def _proper_torsions_for_missing_dofs(
    molecule: Molecule,
    dof_steps: Sequence[PlannedMissingDOFStep],
    graph: Mapping[AtomID, Set[AtomID]],
    mm_params: MMParameterProvider,
    excluded_atom_ids: Iterable[AtomID] = (),
) -> Tuple[Tuple[Tuple[AtomID, AtomID, AtomID, AtomID], Tuple[MMTorsionTerm, ...]], ...]:
    excluded = set(excluded_atom_ids)
    torsions: Dict[
        Tuple[AtomID, AtomID, AtomID, AtomID],
        Tuple[MMTorsionTerm, ...],
    ] = {}
    for dof_step in dof_steps:
        central_a = dof_step.central_bond[0].atom_id()
        central_b = dof_step.central_bond[1].atom_id()
        for terminal_a in graph.get(central_a, set()):
            if terminal_a == central_b:
                continue
            for terminal_b in graph.get(central_b, set()):
                if terminal_b == central_a:
                    continue
                atoms = (terminal_a, central_a, central_b, terminal_b)
                reverse = tuple(reversed(atoms))
                canonical = min(atoms, reverse)
                if any(atom_id in excluded for atom_id in canonical):
                    continue
                try:
                    atom_types = tuple(
                        mm_params.atom_params(
                            get_residue(molecule, atom_id), atom_id.atom_name
                        ).atom_type
                        for atom_id in canonical
                    )
                except KeyError:
                    continue
                terms = tuple(mm_params.proper_dihedral_terms(atom_types))
                if terms:
                    torsions[canonical] = terms
    return tuple(sorted(torsions.items(), key=lambda item: item[0]))


class _SidechainLocalEnergy:
    """Cached local MM evaluator for one cluster of moving branches."""

    def __init__(
        self,
        molecule: Molecule,
        mm_params: MMParameterProvider,
        graph: Mapping[AtomID, Set[AtomID]],
        moving_atom_ids: Sequence[AtomID],
        dof_steps: Sequence[PlannedMissingDOFStep],
        fixed_coords: Mapping[AtomID, Tuple[float, float, float]],
        params_by_atom: Mapping[AtomID, MMAtomParams],
        coordinate_samples: Sequence[Mapping[AtomID, Tuple[float, float, float]]],
        coordinates,
        settings: SidechainOptimizationSettings,
        spatial_index: Optional[MMSpatialIndex] = None,
        collapsed_hydrogens: Iterable[AtomID] = (),
    ) -> None:
        self.molecule = molecule
        self.mm_params = mm_params
        self.graph = graph
        self.moving_atom_ids = tuple(
            atom_id for atom_id in moving_atom_ids if atom_id in params_by_atom
        )
        self.moving_set = set(self.moving_atom_ids)
        self.fixed_coords = fixed_coords
        self.params_by_atom = params_by_atom
        self.coordinates = coordinates
        self.settings = settings
        self.shell_cache: Dict[
            AtomID, Tuple[Set[AtomID], Set[AtomID], Set[AtomID]]
        ] = {}
        if spatial_index is None:
            spatial_atoms = [
                SpatialAtom(
                    atom_id,
                    0.1 * _coord_array(coord),
                    params_by_atom[atom_id],
                )
                for atom_id, coord in fixed_coords.items()
                if atom_id in params_by_atom
            ]
            spatial_index = MMSpatialIndex(
                spatial_atoms,
                settings.preselection_radius_nm,
            )
        self.spatial_index = spatial_index
        self.torsions = _proper_torsions_for_missing_dofs(
            molecule,
            dof_steps,
            graph,
            mm_params,
            collapsed_hydrogens,
        )
        self.moving_pairs: List[Tuple[AtomID, AtomID, float, float]] = []
        for index, atom_a in enumerate(self.moving_atom_ids):
            for atom_b in self.moving_atom_ids[index + 1 :]:
                scaling = _pair_scaling(
                    graph,
                    atom_a,
                    atom_b,
                    mm_params,
                    self.shell_cache,
                )
                if scaling is not None:
                    self.moving_pairs.append((atom_a, atom_b, *scaling))
        self.pair_lists: Dict[AtomID, _VectorizedNonbondedPairs] = {}
        self.candidate_ids: Dict[AtomID, Set[AtomID]] = {}
        self.covered_coords_nm: Dict[AtomID, List[np.ndarray]] = defaultdict(list)
        self.neighbor_expansions = 0
        self.add_coordinate_samples(coordinate_samples)

    def add_coordinate_samples(
        self,
        coordinate_samples: Sequence[Mapping[AtomID, Tuple[float, float, float]]],
    ) -> None:
        for moving_id in self.moving_atom_ids:
            for sample in coordinate_samples:
                if moving_id in sample:
                    self.covered_coords_nm[moving_id].append(
                        0.1 * _coord_array(sample[moving_id])
                    )
            selected: Dict[AtomID, SpatialAtom] = {
                atom.atom_id: atom
                for atom in self._candidate_spatial_atoms(moving_id, coordinate_samples)
            }
            existing = self.candidate_ids.setdefault(moving_id, set())
            for atom_id in existing:
                coord = self.fixed_coords[atom_id]
                selected.setdefault(
                    atom_id,
                    SpatialAtom(
                        atom_id,
                        0.1 * _coord_array(coord),
                        self.params_by_atom[atom_id],
                    ),
                )
            existing.update(selected)
            partners = [selected[atom_id] for atom_id in sorted(selected)]
            self.pair_lists[moving_id] = _VectorizedNonbondedPairs(
                coords_nm=np.asarray(
                    [atom.coord_nm for atom in partners], dtype=float
                ).reshape(-1, 3),
                charges=np.asarray(
                    [atom.params.charge for atom in partners], dtype=float
                ),
                sigmas=np.asarray(
                    [atom.params.sigma for atom in partners], dtype=float
                ),
                epsilons=np.asarray(
                    [atom.params.epsilon for atom in partners], dtype=float
                ),
                scale_lj=np.asarray(
                    [
                        _pair_scaling(
                            self.graph,
                            moving_id,
                            atom.atom_id,
                            self.mm_params,
                            self.shell_cache,
                        )[0]
                        for atom in partners
                    ],
                    dtype=float,
                ),
                scale_qq=np.asarray(
                    [
                        _pair_scaling(
                            self.graph,
                            moving_id,
                            atom.atom_id,
                            self.mm_params,
                            self.shell_cache,
                        )[1]
                        for atom in partners
                    ],
                    dtype=float,
                ),
            )

    def _candidate_spatial_atoms(
        self,
        moving_id: AtomID,
        coordinate_samples: Sequence[Mapping[AtomID, Tuple[float, float, float]]],
    ) -> List[SpatialAtom]:
        selected: Dict[AtomID, SpatialAtom] = {}
        for sample in coordinate_samples:
            if moving_id not in sample:
                continue
            for atom in self.spatial_index.query(
                0.1 * _coord_array(sample[moving_id]),
                self.settings.preselection_radius_nm,
            ):
                if atom.atom_id in self.moving_set:
                    continue
                scaling = _pair_scaling(
                    self.graph,
                    moving_id,
                    atom.atom_id,
                    self.mm_params,
                    self.shell_cache,
                )
                if scaling is not None:
                    selected[atom.atom_id] = atom
        return list(selected.values())

    def ensure_neighbor_coverage(
        self,
        coords: Mapping[AtomID, Tuple[float, float, float]],
    ) -> None:
        """Extend the pair lists only outside the proven swept-volume skin.

        If a new coordinate is within ``preselection-cutoff`` of a previously
        queried coordinate, every partner inside the actual cutoff is
        guaranteed to have been inside the earlier preselection sphere.
        """
        skin_nm = self.settings.preselection_radius_nm - self.settings.cutoff_radius_nm
        needs_extension = False
        for moving_id in self.moving_atom_ids:
            current = 0.1 * _coord_array(coords[moving_id])
            covered = self.covered_coords_nm.get(moving_id, [])
            if not covered:
                needs_extension = True
                break
            covered_array = np.asarray(covered, dtype=float).reshape(-1, 3)
            delta = covered_array - current
            if float(np.min(np.linalg.norm(delta, axis=1))) > skin_nm:
                needs_extension = True
                break
        if needs_extension:
            self.add_coordinate_samples([coords])
            self.neighbor_expansions += 1

    def score(self, values: Mapping[DOFKey, float]) -> float:
        coords = self.coordinates(values)
        # The swept coarse envelope normally contains every relevant partner.
        # Multi-DOF combinations can nevertheless reach positions absent from
        # independent scans. Grow the candidate set lazily before evaluating
        # that state. A newly discovered partner was farther than the 0.9 nm
        # skin in every earlier checked state and therefore contributed exactly
        # zero there under the 0.8 nm cutoff; cached earlier energies stay valid.
        self.ensure_neighbor_coverage(coords)
        energy = 0.0
        if self.settings.include_torsions:
            for atom_ids, terms in self.torsions:
                try:
                    torsion_coords = [
                        coords[atom_id]
                        if atom_id in coords
                        else self.fixed_coords[atom_id]
                        for atom_id in atom_ids
                    ]
                except KeyError:
                    continue
                energy += proper_torsion_energy(
                    compute_dihedral_deg(*torsion_coords), terms
                )

        for moving_id in self.moving_atom_ids:
            pair_list = self.pair_lists[moving_id]
            energy += vectorized_nonbonded_energy(
                0.1 * _coord_array(coords[moving_id]),
                self.params_by_atom[moving_id],
                pair_list.coords_nm,
                pair_list.charges,
                pair_list.sigmas,
                pair_list.epsilons,
                pair_list.scale_lj,
                pair_list.scale_qq,
                switch_radius_nm=self.settings.switch_radius_nm,
                cutoff_radius_nm=self.settings.cutoff_radius_nm,
                include_lj=self.settings.include_lj,
                include_electrostatics=self.settings.include_electrostatics,
            )

        for atom_a, atom_b, scale_lj, scale_qq in self.moving_pairs:
            distance_nm = 0.1 * math.dist(coords[atom_a], coords[atom_b])
            if distance_nm <= 1.0e-12:
                return math.inf
            energy += pair_nonbonded_energy(
                0.1 * _coord_array(coords[atom_a]),
                self.params_by_atom[atom_a],
                0.1 * _coord_array(coords[atom_b]),
                self.params_by_atom[atom_b],
                scale_lj=scale_lj,
                scale_qq=scale_qq,
                switch_radius_nm=self.settings.switch_radius_nm,
                cutoff_radius_nm=self.settings.cutoff_radius_nm,
                include_lj=self.settings.include_lj,
                include_electrostatics=self.settings.include_electrostatics,
            )
        return energy


def _branch_pair_energy_matrix(
    branch_a: _SidechainBranchModel,
    conformers_a: Sequence[
        Tuple[Mapping[DOFKey, float], Mapping[AtomID, Tuple[float, float, float]]]
    ],
    branch_b: _SidechainBranchModel,
    conformers_b: Sequence[
        Tuple[Mapping[DOFKey, float], Mapping[AtomID, Tuple[float, float, float]]]
    ],
    params_by_atom: Mapping[AtomID, MMAtomParams],
    graph: Mapping[AtomID, Set[AtomID]],
    mm_params: MMParameterProvider,
    settings: SidechainOptimizationSettings,
) -> Tuple[float, bool, float]:
    """Return minimum distance, possible-contact flag, and mutual-energy range."""
    energy = np.zeros((len(conformers_a), len(conformers_b)), dtype=float)
    min_distance_nm = math.inf
    possible_contact = False
    shell_cache: Dict[AtomID, Tuple[Set[AtomID], Set[AtomID], Set[AtomID]]] = {}
    considered_pair = False
    for atom_a in branch_a.energy_atom_ids:
        if atom_a not in params_by_atom:
            continue
        coords_a = 0.1 * np.asarray(
            [coords[atom_a] for _values, coords in conformers_a], dtype=float
        )
        for atom_b in branch_b.energy_atom_ids:
            if atom_b not in params_by_atom:
                continue
            scaling = _pair_scaling(
                graph, atom_a, atom_b, mm_params, shell_cache
            )
            if scaling is None:
                continue
            considered_pair = True
            coords_b = 0.1 * np.asarray(
                [coords[atom_b] for _values, coords in conformers_b], dtype=float
            )
            delta = coords_a[:, None, :] - coords_b[None, :, :]
            distances = np.linalg.norm(delta, axis=2)
            pair_min = float(np.min(distances))
            min_distance_nm = min(min_distance_nm, pair_min)
            params_a = params_by_atom[atom_a]
            params_b = params_by_atom[atom_b]
            sigma = 0.5 * (params_a.sigma + params_b.sigma)
            if sigma > 0.0 and pair_min < settings.coupling_sigma_factor * sigma:
                possible_contact = True

            valid = distances > 1.0e-12
            if not np.all(valid):
                possible_contact = True
            safe_distances = np.where(valid, distances, 1.0e-12)
            weights = switch_weights(
                safe_distances,
                settings.switch_radius_nm,
                settings.cutoff_radius_nm,
            )
            pair_energy = np.zeros_like(safe_distances)
            if settings.include_lj:
                epsilon = math.sqrt(max(0.0, params_a.epsilon * params_b.epsilon))
                if sigma > 0.0 and epsilon > 0.0:
                    sr6 = (sigma / safe_distances) ** 6
                    pair_energy += scaling[0] * 4.0 * epsilon * (sr6 * sr6 - sr6)
            if (
                settings.include_electrostatics
                and params_a.charge != 0.0
                and params_b.charge != 0.0
            ):
                pair_energy += (
                    scaling[1]
                    * COULOMB_KJ_MOL_NM
                    * params_a.charge
                    * params_b.charge
                    / safe_distances
                )
            energy += weights * pair_energy
    if not considered_pair:
        return math.inf, False, 0.0
    finite = energy[np.isfinite(energy)]
    energy_range = (
        float(np.max(finite) - np.min(finite)) if finite.size else math.inf
    )
    return min_distance_nm, possible_contact, energy_range


def _sidechain_branch_clusters(
    branches: Sequence[_SidechainBranchModel],
    conformers: Mapping[
        ResidueID,
        Sequence[
            Tuple[
                Mapping[DOFKey, float],
                Mapping[AtomID, Tuple[float, float, float]],
            ]
        ],
    ],
    params_by_atom: Mapping[AtomID, MMAtomParams],
    graph: Mapping[AtomID, Set[AtomID]],
    mm_params: MMParameterProvider,
    settings: SidechainOptimizationSettings,
) -> Tuple[Tuple[Tuple[_SidechainBranchModel, ...], ...], int]:
    adjacency: Dict[ResidueID, Set[ResidueID]] = {
        branch.residue_id: set() for branch in branches
    }
    branch_by_id = {branch.residue_id: branch for branch in branches}
    coupled_pairs = 0
    for index, branch_a in enumerate(branches):
        for branch_b in branches[index + 1 :]:
            min_distance, possible_contact, energy_range = _branch_pair_energy_matrix(
                branch_a,
                conformers[branch_a.residue_id],
                branch_b,
                conformers[branch_b.residue_id],
                params_by_atom,
                graph,
                mm_params,
                settings,
            )
            if min_distance >= settings.cutoff_radius_nm:
                continue
            if not possible_contact and (
                energy_range <= settings.coupling_energy_threshold_kj_mol
            ):
                continue
            adjacency[branch_a.residue_id].add(branch_b.residue_id)
            adjacency[branch_b.residue_id].add(branch_a.residue_id)
            coupled_pairs += 1

    order = {branch.residue_id: i for i, branch in enumerate(branches)}
    visited: Set[ResidueID] = set()
    clusters: List[Tuple[_SidechainBranchModel, ...]] = []
    for branch in branches:
        if branch.residue_id in visited:
            continue
        stack = [branch.residue_id]
        component: List[ResidueID] = []
        visited.add(branch.residue_id)
        while stack:
            residue_id = stack.pop()
            component.append(residue_id)
            for neighbor in sorted(adjacency[residue_id], key=order.get, reverse=True):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        component.sort(key=order.get)
        clusters.append(tuple(branch_by_id[residue_id] for residue_id in component))
    return tuple(clusters), coupled_pairs


def optimize_missing_sidechain_dofs(
    base_molecule: Molecule,
    template: Mapping[str, Any],
    base_remaining_plan: BuildPlan,
    mm_params: MMParameterProvider,
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
    settings: Optional[SidechainOptimizationSettings] = None,
    stats: Optional[SidechainMMOptimizationStats] = None,
) -> Dict[DOFKey, float]:
    """Find deterministic local-MM defaults for all safe open branches.

    The input molecule and plan are never changed.  The returned mapping has
    the same shape as ``get_template_sidechain_dof_defaults`` and can be passed
    directly to ``build_missing_sidechain_gui_payload``.
    """
    started = time.perf_counter()
    settings = settings or SidechainOptimizationSettings()
    settings.validate()
    index = _sidechain_execution_index(base_remaining_plan, execution_index)
    if not index.dof_steps:
        if stats is not None:
            stats.total_seconds = time.perf_counter() - started
        return {}

    defaults = get_template_sidechain_dof_defaults(
        base_molecule,
        template,
        base_remaining_plan,
        execution_index=index,
    )
    geometry_molecule = _copy_molecule_for_local_coordinate_build(
        base_molecule,
        index.residue_groups,
    )
    ordered_atom_steps = _heavy_then_free_rotor_steps(index.atom_steps)
    for step in ordered_atom_steps:
        build_atom_from_step(
            geometry_molecule,
            template,
            step,
            mm_params=None,
            dof_values=defaults,
        )
    # Heavy side-chain geometry is searched with only its own unresolved
    # terminal free-rotor groups represented by parent-centred monopoles.
    # Already fixed free-rotor hydrogens elsewhere in the molecule remain
    # explicit atoms with their actual coordinates, charges, and LJ terms.
    free_hydrogens = _free_rotor_hydrogen_ids_from_steps(index.atom_steps)
    topology_seeds: Set[AtomID] = {
        step.atom_key for step in index.atom_steps
    }
    for dof_step in index.dof_steps:
        topology_seeds.update(ref.atom_id() for ref in dof_step.central_bond)
        topology_seeds.update(
            ref.atom_id() for ref in dof_step.requested_dihedral_atoms
        )
    graph = build_local_mm_bond_graph(
        geometry_molecule,
        mm_params,
        topology_seeds,
    )
    free_hydrogens, collapsed_charges = _collapsed_free_rotor_charges(
        geometry_molecule,
        template,
        mm_params,
        graph,
        free_hydrogens,
    )
    environment = _prepare_sidechain_mm_environment(
        geometry_molecule,
        mm_params,
        free_hydrogens,
        collapsed_charges,
        settings,
    )
    fixed_coords = environment.fixed_coords
    params_by_atom = environment.params_by_atom
    branches = _sidechain_branch_models(index, free_hydrogens)
    preparation_finished = time.perf_counter()

    conformers: Dict[
        ResidueID,
        List[
            Tuple[
                Dict[DOFKey, float],
                Mapping[AtomID, Tuple[float, float, float]],
            ]
        ],
    ] = {}
    branch_generators: Dict[ResidueID, Any] = {}
    for branch in branches:
        dof_order = tuple(step.dof_key for step in branch.dof_steps)
        generator = _sidechain_coordinate_generator(
            geometry_molecule,
            template,
            branch.atom_steps,
            dof_order,
        )
        branch_generators[branch.residue_id] = generator
        conformers[branch.residue_id] = _coarse_sidechain_conformers(
            branch.dof_steps,
            defaults,
            generator,
            settings,
        )
    coarse_sampling_finished = time.perf_counter()

    clusters, coupled_pairs = _sidechain_branch_clusters(
        branches,
        conformers,
        params_by_atom,
        graph,
        mm_params,
        settings,
    )
    clustering_finished = time.perf_counter()
    optimized = dict(defaults)
    total_evaluations = 0
    total_sweeps = 0
    neighbor_rebuilds = 0

    for cluster in clusters:
        dof_steps = tuple(step for branch in cluster for step in branch.dof_steps)
        atom_steps = tuple(step for branch in cluster for step in branch.atom_steps)
        dof_order = tuple(step.dof_key for step in dof_steps)
        moving_ids = tuple(
            atom_id
            for branch in cluster
            for atom_id in branch.energy_atom_ids
        )
        generator = _sidechain_coordinate_generator(
            geometry_molecule,
            template,
            atom_steps,
            dof_order,
        )
        samples = [
            coords
            for branch in cluster
            for _values, coords in conformers[branch.residue_id]
        ]
        energy_model = _SidechainLocalEnergy(
            geometry_molecule,
            mm_params,
            graph,
            moving_ids,
            dof_steps,
            fixed_coords,
            params_by_atom,
            samples,
            generator,
            settings,
            spatial_index=environment.spatial_index,
            collapsed_hydrogens=free_hydrogens,
        )
        affected = {
            dof: tuple(
                step.atom_key
                for step in index.affected_atom_steps_by_dof[dof]
                if step.local_completion_group
                in {branch.residue_id for branch in cluster}
            )
            for dof in dof_order
        }

        initial = {dof: optimized[dof] for dof in dof_order}
        result = optimize_periodic_dofs(
            dof_order,
            initial,
            energy_model.score,
            generator,
            affected,
            settings,
        )
        total_evaluations += result.evaluations
        total_sweeps += result.refinement_sweeps
        neighbor_rebuilds += energy_model.neighbor_expansions
        optimized.update(
            {dof: normalize_angle(value) for dof, value in result.values.items()}
        )

    if stats is not None:
        stats.branch_count = len(branches)
        stats.cluster_count = len(clusters)
        stats.coupled_branch_pairs = coupled_pairs
        stats.coarse_conformations = sum(len(items) for items in conformers.values())
        stats.energy_evaluations = total_evaluations
        stats.refinement_sweeps = total_sweeps
        stats.neighbor_rebuilds = neighbor_rebuilds
        finished = time.perf_counter()
        stats.preparation_seconds = preparation_finished - started
        stats.coarse_sampling_seconds = (
            coarse_sampling_finished - preparation_finished
        )
        stats.clustering_seconds = (
            clustering_finished - coarse_sampling_finished
        )
        stats.optimization_seconds = finished - clustering_finished
        stats.total_seconds = finished - started
    return optimized


def _sidechain_residue_id(
    value: Union[ResidueID, Mapping[str, Any]],
) -> ResidueID:
    if isinstance(value, ResidueID):
        return value
    if not isinstance(value, Mapping):
        raise ValueError(
            "Side-chain residue must be a ResidueID or a structured mapping"
        )
    nested = value.get("residue")
    data = nested if isinstance(nested, Mapping) else value
    try:
        return ResidueID(
            chain_id=str(data["chain_id"]),
            residue_index=int(data["residue_index"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Malformed side-chain residue identity: {value!r}") from exc


def prepare_missing_sidechain_local_optimization(
    working_molecule: Molecule,
    template: Mapping[str, Any],
    base_remaining_plan: BuildPlan,
    residue: Union[ResidueID, Mapping[str, Any]],
    mm_params: MMParameterProvider,
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
    settings: Optional[SidechainOptimizationSettings] = None,
) -> SidechainLocalOptimizationContext:
    """Prepare local MM scoring for one GUI-visible side-chain residue.

    No branch clustering is consulted: all atoms outside the selected residue
    are represented as the fixed, explicit environment found in
    ``working_molecule``.  The working molecule must already contain the first
    complete side-chain preview.
    """

    effective_settings = settings or SidechainOptimizationSettings()
    effective_settings.validate()
    index = _sidechain_execution_index(base_remaining_plan, execution_index)
    residue_id = _sidechain_residue_id(residue)
    if residue_id not in set(index.residue_groups):
        raise ValueError(
            f"Residue {residue_id} is not a safe residue-local open branch"
        )

    missing_preview_atoms = [
        step.atom_key
        for step in index.atom_steps
        if _atom_from_id(working_molecule, step.atom_key).coord is None
    ]
    if missing_preview_atoms:
        raise ValueError(
            "Local side-chain optimization requires the initial complete "
            "side-chain preview; atoms still missing="
            f"{missing_preview_atoms}"
        )

    dof_steps = tuple(
        step
        for step in index.dof_steps
        if step.local_completion_group == residue_id
    )
    atom_steps = tuple(
        step
        for step in index.atom_steps
        if step.local_completion_group == residue_id
    )
    if not dof_steps or not atom_steps:
        raise ValueError(f"Residue {residue_id} has no optimizable local branch")

    geometry_molecule = _copy_molecule_for_local_coordinate_build(
        working_molecule,
        (residue_id,),
    )
    free_hydrogens = _free_rotor_hydrogen_ids_from_steps(atom_steps)
    topology_seeds: Set[AtomID] = {step.atom_key for step in atom_steps}
    for dof_step in dof_steps:
        topology_seeds.update(ref.atom_id() for ref in dof_step.central_bond)
        topology_seeds.update(
            ref.atom_id() for ref in dof_step.requested_dihedral_atoms
        )
    graph = build_local_mm_bond_graph(
        geometry_molecule,
        mm_params,
        topology_seeds,
    )
    free_hydrogens, collapsed_charges = _collapsed_free_rotor_charges(
        geometry_molecule,
        template,
        mm_params,
        graph,
        free_hydrogens,
    )
    environment = _prepare_sidechain_mm_environment(
        geometry_molecule,
        mm_params,
        free_hydrogens,
        collapsed_charges,
        effective_settings,
    )
    branch = _SidechainBranchModel(
        residue_id=residue_id,
        dof_steps=dof_steps,
        atom_steps=atom_steps,
        energy_atom_ids=tuple(
            step.atom_key
            for step in atom_steps
            if step.atom_key not in free_hydrogens
        ),
    )
    dof_order = tuple(step.dof_key for step in dof_steps)
    generator = _sidechain_coordinate_generator(
        geometry_molecule,
        template,
        atom_steps,
        dof_order,
    )
    current_sample = {
        step.atom_key: tuple(
            _atom_from_id(working_molecule, step.atom_key).coord
        )
        for step in atom_steps
    }
    energy_model = _SidechainLocalEnergy(
        geometry_molecule,
        mm_params,
        graph,
        branch.energy_atom_ids,
        dof_steps,
        environment.fixed_coords,
        environment.params_by_atom,
        (current_sample,),
        generator,
        effective_settings,
        spatial_index=environment.spatial_index,
        collapsed_hydrogens=free_hydrogens,
    )
    affected = {
        dof: tuple(
            step.atom_key
            for step in index.affected_atom_steps_by_dof[dof]
            if step.local_completion_group == residue_id
        )
        for dof in dof_order
    }
    return SidechainLocalOptimizationContext(
        residue_id=residue_id,
        dof_order=dof_order,
        affected_atoms_by_dof=affected,
        coordinates=generator,
        energy_model=energy_model,
        settings=effective_settings,
    )


def refine_missing_sidechain_dofs(
    context: SidechainLocalOptimizationContext,
    current_dof_values: Mapping[DOFKey, float],
    *,
    initial_step_degrees: Optional[float] = None,
) -> PeriodicOptimizationResult:
    """Relax one residue from its current GUI DOFs into its local MM basin."""

    expected = set(context.dof_order)
    provided = set(current_dof_values)
    if provided != expected:
        missing = sorted(expected - provided, key=_dof_sort_key)
        extra = sorted(provided - expected, key=_dof_sort_key)
        raise ValueError(
            "Local side-chain DOF values do not match the selected residue; "
            f"missing={missing}, extra={extra}"
        )
    initial = {
        dof: _validated_dof_value(dof, current_dof_values)
        for dof in context.dof_order
    }
    return refine_periodic_dofs(
        context.dof_order,
        initial,
        context.energy_model.score,
        context.coordinates,
        context.affected_atoms_by_dof,
        context.settings,
        initial_step_degrees=initial_step_degrees,
    )


def _dof_key_payload(key: DOFKey) -> Dict[str, Any]:
    return {
        "chain_id": key.chain_id,
        "residue_index": key.residue_index,
        "atom": key.atom_name,
        "rule_index": key.rule_index,
    }


def _resolved_ref_payload(molecule: Molecule, ref: ResolvedRef) -> Dict[str, Any]:
    residue = molecule.chains[ref.chain_id].residues[ref.residue_index]
    return {
        "chain_id": ref.chain_id,
        "residue_index": ref.residue_index,
        "resseq": residue.resseq,
        "icode": residue.icode,
        "ff_resname": residue.ff_resname,
        "atom": ref.atom_name,
    }


def _residue_payload(molecule: Molecule, residue_id: ResidueID) -> Dict[str, Any]:
    residue = molecule.chains[residue_id.chain_id].residues[
        residue_id.residue_index
    ]
    return {
        "chain_id": residue_id.chain_id,
        "residue_index": residue_id.residue_index,
        "resseq": residue.resseq,
        "icode": residue.icode,
        "ff_resname": residue.ff_resname,
    }


def build_missing_sidechain_gui_payload(
    molecule: Molecule,
    remaining_plan: BuildPlan,
    default_dof_values: Mapping[DOFKey, float],
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
) -> Dict[str, Any]:
    """Distill safe local branches into a minimal JSON-compatible GUI payload.

    DOFs and residue groups follow their actual order in `remaining_plan.steps`,
    not the set-like order of `ResidueLocalCompletionGroup.dof_keys`.
    `value_degrees` is initialized from `default_degrees`; the GUI may update
    only the former while retaining the latter for reset-to-default behavior.
    """
    groups = _open_branch_groups(remaining_plan)
    index = _sidechain_execution_index(remaining_plan, execution_index)
    ordered_dof_steps = list(index.dof_steps)
    expected_dofs = {step.dof_key for step in ordered_dof_steps}
    provided_dofs = set(default_dof_values)
    if provided_dofs != expected_dofs:
        missing = sorted(expected_dofs - provided_dofs, key=_dof_sort_key)
        extra = sorted(provided_dofs - expected_dofs, key=_dof_sort_key)
        raise ValueError(
            "Default side-chain DOF values do not match the remaining safe "
            f"branches; missing={missing}, extra={extra}"
        )

    steps_by_residue: Dict[ResidueID, List[PlannedMissingDOFStep]] = defaultdict(list)
    residue_order: List[ResidueID] = []
    for step in ordered_dof_steps:
        residue_id = step.local_completion_group
        if residue_id is None:
            raise RuntimeError(f"Safe local DOF lacks a completion group: {step.dof_key}")
        if residue_id not in steps_by_residue:
            residue_order.append(residue_id)
        steps_by_residue[residue_id].append(step)

    atom_names_by_residue: Dict[ResidueID, List[str]] = defaultdict(list)
    for step in index.atom_steps:
        residue_id = step.local_completion_group
        if residue_id in groups:
            atom_names_by_residue[residue_id].append(step.atom_key.atom_name)

    sidechains: List[Dict[str, Any]] = []
    for residue_id in residue_order:
        dofs: List[Dict[str, Any]] = []
        for step in steps_by_residue[residue_id]:
            default = _validated_dof_value(step.dof_key, default_dof_values)
            dofs.append(
                {
                    "dof_key": _dof_key_payload(step.dof_key),
                    "dihedral_atoms": [
                        _resolved_ref_payload(molecule, ref)
                        for ref in step.requested_dihedral_atoms
                    ],
                    "default_degrees": default,
                    "value_degrees": default,
                }
            )
        sidechains.append(
            {
                "residue": _residue_payload(molecule, residue_id),
                "atoms": atom_names_by_residue[residue_id],
                "dofs": dofs,
            }
        )
    return {"sidechains": sidechains}


def _dof_key_from_payload(data: Mapping[str, Any]) -> DOFKey:
    try:
        return DOFKey(
            chain_id=str(data["chain_id"]),
            residue_index=int(data["residue_index"]),
            atom_name=str(data["atom"]),
            rule_index=int(data["rule_index"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Malformed GUI DOF key: {data!r}") from exc


def _sidechain_dof_values_from_payload(
    remaining_plan: BuildPlan,
    sidechain_data: Mapping[str, Any],
    execution_index: Optional[SidechainExecutionIndex] = None,
) -> Dict[DOFKey, float]:
    groups = _open_branch_groups(remaining_plan)
    index = _sidechain_execution_index(remaining_plan, execution_index)
    expected_by_residue: Dict[ResidueID, Set[DOFKey]] = defaultdict(set)
    for step in index.dof_steps:
        if step.local_completion_group is None:
            raise RuntimeError(f"Safe local DOF lacks a completion group: {step.dof_key}")
        expected_by_residue[step.local_completion_group].add(step.dof_key)

    raw_sidechains = sidechain_data.get("sidechains")
    if not isinstance(raw_sidechains, Sequence) or isinstance(
        raw_sidechains, (str, bytes)
    ):
        raise ValueError("Side-chain GUI payload must contain a 'sidechains' list")

    values: Dict[DOFKey, float] = {}
    seen_residues: Set[ResidueID] = set()
    for raw_sidechain in raw_sidechains:
        if not isinstance(raw_sidechain, Mapping):
            raise ValueError(f"Malformed side-chain entry: {raw_sidechain!r}")
        raw_residue = raw_sidechain.get("residue")
        if not isinstance(raw_residue, Mapping):
            raise ValueError(f"Side-chain entry lacks residue identity: {raw_sidechain!r}")
        try:
            residue_id = ResidueID(
                chain_id=str(raw_residue["chain_id"]),
                residue_index=int(raw_residue["residue_index"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Malformed side-chain residue identity: {raw_residue!r}") from exc
        if residue_id not in groups:
            raise ValueError(f"GUI payload references a non-open side-chain group: {residue_id}")
        if residue_id in seen_residues:
            raise ValueError(f"GUI payload repeats side-chain residue: {residue_id}")
        seen_residues.add(residue_id)

        raw_dofs = raw_sidechain.get("dofs")
        if not isinstance(raw_dofs, Sequence) or isinstance(raw_dofs, (str, bytes)):
            raise ValueError(f"Side-chain {residue_id} must contain a 'dofs' list")
        for raw_dof in raw_dofs:
            if not isinstance(raw_dof, Mapping):
                raise ValueError(f"Malformed DOF entry for {residue_id}: {raw_dof!r}")
            raw_key = raw_dof.get("dof_key")
            if not isinstance(raw_key, Mapping):
                raise ValueError(f"DOF entry lacks a structured key: {raw_dof!r}")
            dof_key = _dof_key_from_payload(raw_key)
            if dof_key not in expected_by_residue[residue_id]:
                raise ValueError(
                    f"DOF {dof_key} does not belong to safe side-chain {residue_id}"
                )
            if dof_key in values:
                raise ValueError(f"GUI payload repeats missing DOF {dof_key}")
            if "value_degrees" not in raw_dof:
                raise ValueError(f"DOF {dof_key} lacks value_degrees")
            values[dof_key] = _validated_dof_value(
                dof_key,
                {dof_key: raw_dof["value_degrees"]},
            )

    expected_dofs = set().union(*expected_by_residue.values()) if expected_by_residue else set()
    if seen_residues != set(expected_by_residue) or set(values) != expected_dofs:
        missing_residues = sorted(set(expected_by_residue) - seen_residues)
        missing_dofs = sorted(expected_dofs - set(values), key=_dof_sort_key)
        raise ValueError(
            "GUI payload does not provide every safe local branch; "
            f"missing_residues={missing_residues}, missing_dofs={missing_dofs}"
        )
    return values


def complete_missing_sidechains(
    base_molecule: Molecule,
    template: Mapping[str, Any],
    base_remaining_plan: BuildPlan,
    sidechain_data: Mapping[str, Any],
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
) -> Molecule:
    """Build every safe residue-local open branch from an immutable baseline.

    Selection, default generation, MM optimization, session state and GUI
    communication intentionally live outside this light executor. The molecule
    is copied once; the source plan and execution index remain reusable.
    """
    if base_remaining_plan.steps and not isinstance(
        base_remaining_plan.steps[0], PlannedMissingDOFStep
    ):
        raise ValueError(
            "complete_missing_sidechains requires a plan already executed up "
            "to its first missing DOF"
        )
    index = _sidechain_execution_index(base_remaining_plan, execution_index)
    dof_values = _sidechain_dof_values_from_payload(
        base_remaining_plan,
        sidechain_data,
        index,
    )
    molecule = _copy_molecule_for_coordinate_build(base_molecule)
    selected_steps = _heavy_then_free_rotor_steps(index.atom_steps)
    active_free_rotors = _free_rotor_hydrogen_ids_from_steps(selected_steps)
    sidechain_free_rotor_settings = replace(
        free_rotor_settings or FreeRotorSearchSettings(),
        collapse_other_free_rotors=False,
    )
    mm_cache = (
        FreeRotorMMCache(
            build_local_mm_bond_graph(
                molecule,
                mm_params,
                (step.atom_key for step in selected_steps),
            ),
            free_rotor_hydrogens_to_collapse=active_free_rotors,
        )
        if mm_params is not None
        else None
    )
    for step in selected_steps:
        build_atom_from_step(
            molecule,
            template,
            step,
            mm_params=mm_params,
            free_rotor_settings=sidechain_free_rotor_settings,
            mm_cache=mm_cache,
            dof_values=dof_values,
        )
    return molecule


def _changed_dof_key(value: Union[DOFKey, Mapping[str, Any]]) -> DOFKey:
    if isinstance(value, DOFKey):
        return value
    if not isinstance(value, Mapping):
        raise ValueError(f"Changed DOF must be a DOFKey or mapping, got {value!r}")
    nested = value.get("dof_key")
    return _dof_key_from_payload(nested if isinstance(nested, Mapping) else value)


def _atom_coordinate_payload(
    molecule: Molecule,
    atom_id: AtomID,
) -> Dict[str, Any]:
    residue = get_residue(molecule, atom_id)
    atom = residue.atoms[atom_id.atom_name]
    if atom.coord is None:
        raise RuntimeError(f"Updated atom has no coordinates: {atom_id}")
    return {
        "atom": {
            "chain_id": atom_id.chain_id,
            "residue_index": atom_id.residue_index,
            "resseq": residue.resseq,
            "icode": residue.icode,
            "ff_resname": residue.ff_resname,
            "atom_name": atom_id.atom_name,
        },
        "coord": [float(value) for value in atom.coord],
    }


def update_missing_sidechains(
    working_molecule: Molecule,
    template: Mapping[str, Any],
    base_remaining_plan: BuildPlan,
    sidechain_data: Mapping[str, Any],
    changed_dof_keys: Iterable[Union[DOFKey, Mapping[str, Any]]],
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
) -> Tuple[Molecule, Dict[str, Any]]:
    """Transactionally rebuild only atoms affected by changed side-chain DOFs.

    The working molecule is updated in-place and returned for convenient API
    symmetry. If any coordinate step fails, all atoms selected for this update
    are restored to their exact pre-call state and no coordinate patch is
    returned. The canonical remaining plan and execution index are never
    modified.
    """
    index = _sidechain_execution_index(base_remaining_plan, execution_index)
    dof_values = _sidechain_dof_values_from_payload(
        base_remaining_plan,
        sidechain_data,
        index,
    )
    changed: Set[DOFKey] = {
        _changed_dof_key(value)
        for value in changed_dof_keys
    }
    known_dofs = set(index.affected_atom_steps_by_dof)
    unknown = sorted(changed - known_dofs, key=_dof_sort_key)
    if unknown:
        raise ValueError(f"Incremental update references non-side-chain DOFs: {unknown}")
    if not changed:
        return working_molecule, {"updated_atoms": []}

    affected_atoms: Set[AtomID] = set()
    for dof_key in changed:
        affected_atoms.update(
            step.atom_key
            for step in index.affected_atom_steps_by_dof[dof_key]
        )
    selected_steps = _heavy_then_free_rotor_steps(
        step
        for step in index.atom_steps
        if step.atom_key in affected_atoms
    )

    snapshots: Dict[AtomID, Tuple[Any, ...]] = {}
    for step in selected_steps:
        atom = _atom_from_id(working_molecule, step.atom_key)
        snapshots[step.atom_key] = (
            atom.coord,
            atom.built,
            atom.build_source,
            atom.build_rule_index,
            atom.occupancy,
            atom.bfactor,
        )

    active_free_rotors = _free_rotor_hydrogen_ids_from_steps(selected_steps)
    sidechain_free_rotor_settings = replace(
        free_rotor_settings or FreeRotorSearchSettings(),
        collapse_other_free_rotors=False,
    )
    mm_cache = (
        FreeRotorMMCache(
            build_local_mm_bond_graph(
                working_molecule,
                mm_params,
                (step.atom_key for step in selected_steps),
            ),
            free_rotor_hydrogens_to_collapse=active_free_rotors,
        )
        if mm_params is not None
        else None
    )
    try:
        for step in selected_steps:
            rebuild_atom_from_step(
                working_molecule,
                template,
                step,
                mm_params=mm_params,
                free_rotor_settings=sidechain_free_rotor_settings,
                mm_cache=mm_cache,
                dof_values=dof_values,
            )
        patch = {
            "updated_atoms": [
                _atom_coordinate_payload(working_molecule, step.atom_key)
                for step in selected_steps
            ]
        }
    except Exception:
        for atom_id, snapshot in snapshots.items():
            atom = _atom_from_id(working_molecule, atom_id)
            (
                atom.coord,
                atom.built,
                atom.build_source,
                atom.build_rule_index,
                atom.occupancy,
                atom.bfactor,
            ) = snapshot
        raise
    return working_molecule, patch


def optimize_missing_sidechain_preview(
    working_molecule: Molecule,
    template: Mapping[str, Any],
    base_remaining_plan: BuildPlan,
    sidechain_data: Mapping[str, Any],
    residue: Union[ResidueID, Mapping[str, Any]],
    mm_params: MMParameterProvider,
    *,
    execution_index: Optional[SidechainExecutionIndex] = None,
    settings: Optional[SidechainOptimizationSettings] = None,
    optimization_context: Optional[SidechainLocalOptimizationContext] = None,
    initial_step_degrees: Optional[float] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
) -> Tuple[Molecule, Dict[str, Any]]:
    """Locally optimize one GUI side-chain and return its coordinate patch.

    The current slider values are read from ``sidechain_data``.  Only the DOFs
    and atoms of the selected residue are refined/rebuilt; every other branch
    remains fixed in its current working-molecule geometry.  ``sidechain_data``
    is not modified, so its ``default_degrees`` remain a stable reset target.

    A prepared context may be reused while only this residue changes.  It must
    be discarded after any coordinate outside this residue is modified.
    """

    index = _sidechain_execution_index(base_remaining_plan, execution_index)
    residue_id = _sidechain_residue_id(residue)
    all_values = _sidechain_dof_values_from_payload(
        base_remaining_plan,
        sidechain_data,
        index,
    )
    selected_keys = tuple(
        step.dof_key
        for step in index.dof_steps
        if step.local_completion_group == residue_id
    )
    if not selected_keys:
        raise ValueError(
            f"Residue {residue_id} is not a safe optimizable side-chain branch"
        )
    current_values = {key: all_values[key] for key in selected_keys}
    context = optimization_context or prepare_missing_sidechain_local_optimization(
        working_molecule,
        template,
        base_remaining_plan,
        residue_id,
        mm_params,
        execution_index=index,
        settings=settings,
    )
    if context.residue_id != residue_id:
        raise ValueError(
            "The supplied local-optimization context belongs to "
            f"{context.residue_id}, not {residue_id}"
        )
    initial_energy = float(context.energy_model.score(current_values))
    result = refine_missing_sidechain_dofs(
        context,
        current_values,
        initial_step_degrees=initial_step_degrees,
    )
    optimized_values = {
        key: normalize_angle(value) for key, value in result.values.items()
    }

    updated_data = copy.deepcopy(sidechain_data)
    selected_entry: Optional[Mapping[str, Any]] = None
    for sidechain in updated_data.get("sidechains", []):
        if not isinstance(sidechain, Mapping):
            continue
        raw_residue = sidechain.get("residue")
        if isinstance(raw_residue, Mapping) and _sidechain_residue_id(
            raw_residue
        ) == residue_id:
            selected_entry = sidechain
            break
    if selected_entry is None:
        raise ValueError(f"GUI payload omits selected side-chain {residue_id}")
    entries_by_key: Dict[DOFKey, Any] = {}
    for entry in selected_entry.get("dofs", []):
        if isinstance(entry, Mapping) and isinstance(entry.get("dof_key"), Mapping):
            entries_by_key[_dof_key_from_payload(entry["dof_key"])] = entry
    for dof_key, value in optimized_values.items():
        if dof_key not in entries_by_key:
            raise ValueError(
                f"GUI payload omits selected side-chain DOF {dof_key}"
            )
        entries_by_key[dof_key]["value_degrees"] = value

    changed = tuple(
        dof_key
        for dof_key in selected_keys
        if abs(normalize_angle(optimized_values[dof_key] - current_values[dof_key]))
        > 1.0e-10
    )
    _molecule, patch = update_missing_sidechains(
        working_molecule,
        template,
        base_remaining_plan,
        updated_data,
        changed,
        execution_index=index,
        mm_params=mm_params,
        free_rotor_settings=free_rotor_settings,
    )
    response = {
        "residue": _residue_payload(working_molecule, residue_id),
        "dofs": [
            {
                "dof_key": _dof_key_payload(dof_key),
                "value_degrees": float(optimized_values[dof_key]),
            }
            for dof_key in selected_keys
        ],
        "initial_energy_kj_mol": initial_energy,
        "optimized_energy_kj_mol": float(result.energy_kj_mol),
        "energy_change_kj_mol": float(result.energy_kj_mol - initial_energy),
        "evaluations": int(result.evaluations),
        "refinement_sweeps": int(result.refinement_sweeps),
        "updated_atoms": patch["updated_atoms"],
    }
    return working_molecule, response


def execution_report(original_plan: BuildPlan, remaining_plan: BuildPlan) -> str:
    before = summarize_build_plan(original_plan)
    after = summarize_build_plan(remaining_plan)
    executed_atom_steps = before["atom_build_steps"] - after["atom_build_steps"]
    lines = []
    lines.append("Molecule build execution report")
    lines.append("===============================\n")
    lines.append(f"Executed atom build steps: {executed_atom_steps}")
    lines.append(f"Remaining plan steps: {after['total_steps']}")
    lines.append(f"Remaining atom build steps: {after['atom_build_steps']}")
    lines.append(f"Remaining missing DOF steps: {after['missing_dof_steps']}")
    lines.append(f"Remaining unresolved atoms: {after['unresolved_atoms']}")
    if remaining_plan.steps and isinstance(remaining_plan.steps[0], PlannedMissingDOFStep):
        lines.append("Stopped at next missing DOF: yes")
    else:
        lines.append("Stopped at next missing DOF: no")
    return "\n".join(lines) + "\n"

# -----------------------------------------------------------------------------
# Serialization/reporting helpers
# -----------------------------------------------------------------------------

def atom_label(mol: Molecule, atom_id: AtomID) -> str:
    res = get_residue(mol, atom_id)
    return f"{atom_id.chain_id}:{res.resseq}{res.icode or ''} {res.ff_resname} {atom_id.atom_name}"


def ref_label(mol: Molecule, ref: ResolvedRef) -> str:
    res = mol.chains[ref.chain_id].residues[ref.residue_index]
    return f"{ref.chain_id}:{res.resseq}{res.icode or ''} {res.ff_resname} {ref.atom_name}"


def summarize_build_plan(plan: BuildPlan) -> Dict[str, Any]:
    atom_steps = [s for s in plan.steps if isinstance(s, PlannedAtomBuildStep)]
    dof_steps = [s for s in plan.steps if isinstance(s, PlannedMissingDOFStep)]
    by_class: Dict[str, int] = {}
    by_source: Dict[str, int] = {}
    for s in atom_steps:
        by_class[s.torsion_class] = by_class.get(s.torsion_class, 0) + 1
        by_source[s.torsion_source] = by_source.get(s.torsion_source, 0) + 1
    return {
        "total_steps": len(plan.steps),
        "atom_build_steps": len(atom_steps),
        "missing_dof_steps": len(dof_steps),
        "requirements": len(plan.requirements),
        "unresolved_atoms": len(plan.unresolved_atoms),
        "atom_steps_by_torsion_class": by_class,
        "atom_steps_by_torsion_source": by_source,
    }



def plan_summary(plan: BuildPlan) -> Dict[str, Any]:
    """Backward-compatible alias for summarize_build_plan."""
    return summarize_build_plan(plan)

def build_plan_report(mol: Molecule, plan: BuildPlan, *, include_steps: bool = False) -> str:
    s = summarize_build_plan(plan)
    lines: List[str] = []
    lines.append("Molecule build planner report")
    lines.append("=============================\n")
    lines.append(f"Total plan steps: {s['total_steps']}")
    lines.append(f"Atom build steps: {s['atom_build_steps']}")
    lines.append(f"Missing DOF steps: {s['missing_dof_steps']}")
    lines.append(f"Unresolved atoms: {s['unresolved_atoms']}")
    lines.append("")
    lines.append("Atom steps by torsion class:")
    for k in ("rigid", "derived_rotatable", "free_rotor_hydrogen"):
        lines.append(f"  {k}: {s['atom_steps_by_torsion_class'].get(k, 0)}")
    lines.append("")
    lines.append("Atom steps by torsion source:")
    for k in sorted(s["atom_steps_by_torsion_source"]):
        lines.append(f"  {k}: {s['atom_steps_by_torsion_source'][k]}")
    lines.append("")
    if plan.requirements:
        lines.append("Required missing DOFs:")
        for i, req in enumerate(plan.requirements, 1):
            atoms = " - ".join(ref_label(mol, r) for r in req.requested_dihedral_atoms)
            cb = " - ".join(ref_label(mol, r) for r in req.central_bond)
            lines.append(f"  {i}. DOF {req.dof_key.chain_id}:{req.dof_key.residue_index}:{req.dof_key.atom_name}[rule {req.dof_key.rule_index}]")
            lines.append(f"     central bond: {cb}")
            lines.append(f"     requested dihedral: {atoms}")
            lines.append(f"     unlocks first: {atom_label(mol, req.reason_atom)} via rule {req.reason_rule_index}")
            if req.prerequisite_dofs:
                lines.append(
                    "     prerequisite DOFs: "
                    + ", ".join(
                        f"{key.chain_id}:{key.residue_index}:{key.atom_name}[{key.rule_index}]"
                        for key in req.prerequisite_dofs
                    )
                )
            lines.append(
                f"     completion classification: {req.completion_classification}"
            )
    else:
        lines.append("Required missing DOFs: none")
    if plan.local_completions:
        lines.append("")
        lines.append("Residue-local completion groups:")
        for group in plan.local_completions:
            residue = mol.chains[group.residue_id.chain_id].residues[
                group.residue_id.residue_index
            ]
            lines.append(
                f"  {group.residue_id.chain_id}:{residue.resseq}{residue.icode or ''} "
                f"{residue.ff_resname}: {group.classification}; "
                f"DOFs={len(group.dof_keys)}, atoms={len(group.atom_keys)}, "
                f"anchors={len(group.anchor_atoms)}"
            )
    if include_steps:
        lines.append("")
        lines.append("Plan steps:")
        for i, step in enumerate(plan.steps, 1):
            if isinstance(step, PlannedMissingDOFStep):
                atoms = " - ".join(ref_label(mol, r) for r in step.requested_dihedral_atoms)
                lines.append(f"  {i:5d}. REQUEST_DOF {atoms}")
            else:
                extra = ""
                if step.torsion_source == "observed_member":
                    extra = f" observed_member={step.observed_member_index}"
                elif step.torsion_source == "supplied_dof":
                    extra = f" dof={step.dof_key}"
                elif step.torsion_source == "free_rotor_group_phase":
                    extra = f" anchor={atom_label(mol, step.free_rotor_anchor) if step.free_rotor_anchor else None} phase_offset={step.phase_offset}"
                lines.append(
                    f"  {i:5d}. BUILD {atom_label(mol, step.atom_key)} rule={step.rule_index} "
                    f"class={step.torsion_class} source={step.torsion_source}{extra}"
                )
    return "\n".join(lines) + "\n"


def plan_to_dict(mol: Molecule, plan: BuildPlan) -> Dict[str, Any]:
    def atom_id_dict(a: AtomID) -> Dict[str, Any]:
        res = get_residue(mol, a)
        return {
            "chain_id": a.chain_id,
            "residue_index": a.residue_index,
            "resseq": res.resseq,
            "icode": res.icode,
            "ff_resname": res.ff_resname,
            "atom": a.atom_name,
        }

    def ref_dict(r: ResolvedRef) -> Dict[str, Any]:
        res = mol.chains[r.chain_id].residues[r.residue_index]
        return {
            "chain_id": r.chain_id,
            "residue_index": r.residue_index,
            "resseq": res.resseq,
            "icode": res.icode,
            "ff_resname": res.ff_resname,
            "atom": r.atom_name,
        }

    def dof_key_dict(key: DOFKey) -> Dict[str, Any]:
        return {
            "chain_id": key.chain_id,
            "residue_index": key.residue_index,
            "atom": key.atom_name,
            "rule_index": key.rule_index,
        }

    def residue_id_dict(key: ResidueID) -> Dict[str, Any]:
        residue = mol.chains[key.chain_id].residues[key.residue_index]
        return {
            "chain_id": key.chain_id,
            "residue_index": key.residue_index,
            "resseq": residue.resseq,
            "icode": residue.icode,
            "ff_resname": residue.ff_resname,
        }

    out_steps: List[Dict[str, Any]] = []
    for step in plan.steps:
        if isinstance(step, PlannedMissingDOFStep):
            item = {
                "type": "missing_dof",
                "dof_key": dof_key_dict(step.dof_key),
                "central_bond": [ref_dict(r) for r in step.central_bond],
                "requested_dihedral_atoms": [ref_dict(r) for r in step.requested_dihedral_atoms],
                "torsion_group_index": step.torsion_group_index,
                "requested_member_index": step.requested_member_index,
                "reason_atom": atom_id_dict(step.reason_atom),
                "reason_rule_index": step.reason_rule_index,
            }
            if step.prerequisite_dofs:
                item["prerequisite_dofs"] = [
                    dof_key_dict(key) for key in step.prerequisite_dofs
                ]
            if step.local_completion_group is not None:
                item["local_completion_group"] = residue_id_dict(
                    step.local_completion_group
                )
            if step.completion_classification != "unclassified":
                item["completion_classification"] = (
                    step.completion_classification
                )
            out_steps.append(item)
        else:
            item = {
                "type": "build_atom",
                "atom": atom_id_dict(step.atom_key),
                "rule_index": step.rule_index,
                "torsion_class": step.torsion_class,
                "torsion_source": step.torsion_source,
                "observed_member_index": step.observed_member_index,
                "dof_key": None if step.dof_key is None else dof_key_dict(step.dof_key),
                "free_rotor_anchor": None if step.free_rotor_anchor is None else atom_id_dict(step.free_rotor_anchor),
                "phase_offset": step.phase_offset,
                "free_rotor_group_atoms": list(step.free_rotor_group_atoms),
            }
            if step.required_dofs:
                item["required_dofs"] = [
                    dof_key_dict(key) for key in step.required_dofs
                ]
            if step.local_completion_group is not None:
                item["local_completion_group"] = residue_id_dict(
                    step.local_completion_group
                )
            out_steps.append(item)
    output = {"summary": summarize_build_plan(plan), "steps": out_steps}
    if plan.local_completions:
        output["local_completions"] = [
            {
                "residue": residue_id_dict(group.residue_id),
                "classification": group.classification,
                "dof_keys": [dof_key_dict(key) for key in group.dof_keys],
                "atom_keys": [atom_id_dict(key) for key in group.atom_keys],
                "anchor_atoms": [atom_id_dict(key) for key in group.anchor_atoms],
                "non_bridge_dofs": [
                    dof_key_dict(key) for key in group.non_bridge_dofs
                ],
                "external_pending_atoms": [
                    atom_id_dict(key) for key in group.external_pending_atoms
                ],
            }
            for group in plan.local_completions
        ]
    return output
