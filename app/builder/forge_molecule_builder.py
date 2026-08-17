#!/usr/bin/env python3
"""Planning and coordinate execution for deterministic missing-atom building.

The module contains:
  * a symbolic planner that converts missing atoms into build steps, including
    missing-DOF boundaries;
  * a coordinate executor that applies build steps until the next missing DOF;
  * a narrow MM parameter/scoring layer used only for free-rotor hydrogen search.
"""

from __future__ import annotations

import heapq
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

from forge_molecule_parser import Molecule, Residue, Chain


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


@dataclass
class PlannedMissingDOFStep:
    dof_key: DOFKey
    central_bond: Tuple[ResolvedRef, ResolvedRef]
    requested_dihedral_atoms: Tuple[ResolvedRef, ResolvedRef, ResolvedRef, ResolvedRef]
    torsion_group_index: int
    requested_member_index: int
    reason_atom: AtomID
    reason_rule_index: int


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


PlanStep = Union[PlannedAtomBuildStep, PlannedMissingDOFStep]


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
# Planner
# -----------------------------------------------------------------------------

def plan_build_steps(
    mol: Molecule,
    template: Mapping[str, Any],
    *,
    stats: Optional[PlannerStats] = None,
) -> BuildPlan:
    pending_atoms = collect_pending_atoms(mol)
    available_atoms = collect_available_atoms(mol)
    available_dofs: Set[DOFKey] = set()
    steps: List[PlanStep] = []
    requirements: List[PlannedMissingDOFStep] = []
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
            requirements.append(dof_step)
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

    return BuildPlan(
        steps=steps,
        unresolved_atoms=set(),
        requirements=requirements,
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

COULOMB_KJ_MOL_NM = 138.935458  # kJ mol^-1 nm e^-2


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
class FreeRotorMMCache:
    """Topology data shared by all free-rotor searches in one build execution."""

    bond_graph: Dict[AtomID, Set[AtomID]]
    topological_distances: Dict[Tuple[AtomID, AtomID, int], Optional[int]] = field(
        default_factory=dict
    )
    free_rotor_hydrogens: Optional[Set[AtomID]] = None
    collapsed_charge_by_parent: Optional[Dict[AtomID, float]] = None
    spatial_cells: Optional[Dict[Tuple[int, int, int], List["_EnvAtom"]]] = None
    spatial_cell_size_nm: Optional[float] = None


@dataclass
class MMParameterProvider:
    """Thin MM parameter layer for builder-local free-rotor searches.

    The provider stores only what the builder needs: per-residue atom type/charge,
    atom-type LJ parameters, proper dihedral types, residue bonds, and default
    1-4 scaling factors. It is intentionally not a general force-field engine.
    """

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
        fudge_lj, fudge_qq = _read_umfff_defaults(force_field_file) if force_field_file is not None else (0.5, 0.8333)

        residues: Dict[str, Dict[str, MMAtomParams]] = {}
        for resname, atoms in residues_raw.items():
            residues[resname] = {}
            for atom_name, atom_type, charge in atoms:
                if atom_type not in atomtypes:
                    raise ValueError(f"Atom type {atom_type!r} for {resname}:{atom_name} missing in nonbonded parameters")
                sigma, epsilon = atomtypes[atom_type]
                residues[resname][atom_name] = MMAtomParams(atom_type, charge, sigma, epsilon)
        return cls(residues=residues, residue_bonds=bonds, dihedraltypes=dihedraltypes, fudge_lj=fudge_lj, fudge_qq=fudge_qq)

    @classmethod
    def from_ff_ida_object(cls, ff_obj: Any, force_field_file: Optional[Any] = None) -> "MMParameterProvider":
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
                raise ValueError("FF_IDA object does not contain nonbonded sigma/R and eps arrays")
            for name, atype, q, sig, eps in zip(names, types, charges, sigmas, epsilons):
                residues[resname][name] = MMAtomParams(str(atype), float(q), float(sig), float(eps))
            bonds[resname] = [tuple(b[:2]) for b in unit.get("bonds", [])]

        dihedraltypes: Dict[Tuple[str, str, str, str], List[MMTorsionTerm]] = {}
        for key, terms in ff_obj.b.get("dihedraltypes", {}).items():
            dihedraltypes[tuple(key)] = [MMTorsionTerm(float(t[0]), float(t[1]), int(float(t[2]))) for t in terms]
        fudge_lj, fudge_qq = _read_umfff_defaults(force_field_file) if force_field_file is not None else (0.5, 0.8333)
        return cls(residues=residues, residue_bonds=bonds, dihedraltypes=dihedraltypes, fudge_lj=fudge_lj, fudge_qq=fudge_qq)

    def atom_params(self, residue: Residue, atom_name: str) -> MMAtomParams:
        try:
            return self.residues[residue.ff_resname][atom_name]
        except KeyError as exc:
            raise KeyError(f"MM parameters missing for {residue.ff_resname}:{atom_name}") from exc

    def proper_dihedral_terms(self, atom_types: Tuple[str, str, str, str]) -> List[MMTorsionTerm]:
        key = _canonical_proper_dihedral_key(atom_types)
        return self.dihedraltypes.get(key, [])


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


def _iter_text_lines(source: Any) -> Iterable[str]:
    if hasattr(source, "read"):
        try:
            source.seek(0)
        except Exception:
            pass
        for line in source:
            if isinstance(line, bytes):
                yield line.decode("utf-8")
            else:
                yield line
        return
    with open(source, "r", encoding="utf-8") as handle:
        for line in handle:
            yield line


def _strip_umfff_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def _read_umfff_residue_lib(path: Any) -> Tuple[Dict[str, List[Tuple[str, str, float]]], Dict[str, List[Tuple[str, str]]]]:
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
            residues[current_residue].append((fields[0], fields[1], float(fields[2])))
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


def _read_umfff_bonded_dihedrals(path: Any) -> Dict[Tuple[str, str, str, str], List[MMTorsionTerm]]:
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
        atypes = (fields[0], fields[1], fields[2], fields[3])
        key = _canonical_proper_dihedral_key(atypes)
        dihedrals.setdefault(key, []).append(MMTorsionTerm(float(fields[6]), float(fields[5]), int(float(fields[7]))))
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


def _canonical_proper_dihedral_key(atom_types: Tuple[str, str, str, str]) -> Tuple[str, str, str, str]:
    rev = (atom_types[3], atom_types[2], atom_types[1], atom_types[0])
    return min(tuple(atom_types), rev)


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


def _topological_distance_upto(graph: Mapping[AtomID, Set[AtomID]], start: AtomID, target: AtomID, max_depth: int = 3) -> Optional[int]:
    if start == target:
        return 0
    visited = {start}
    frontier = {start}
    for depth in range(1, max_depth + 1):
        nxt: Set[AtomID] = set()
        for node in frontier:
            for nb in graph.get(node, set()):
                if nb == target:
                    return depth
                if nb not in visited:
                    visited.add(nb)
                    nxt.add(nb)
        frontier = nxt
        if not frontier:
            break
    return None


def _cached_topological_distance_upto(
    cache: FreeRotorMMCache,
    start: AtomID,
    target: AtomID,
    max_depth: int = 3,
) -> Optional[int]:
    first, second = min((start, target), (target, start))
    key = (first, second, max_depth)
    if key not in cache.topological_distances:
        cache.topological_distances[key] = _topological_distance_upto(
            cache.bond_graph,
            start,
            target,
            max_depth,
        )
    return cache.topological_distances[key]


def _switch_weight(r_nm: float, settings: FreeRotorSearchSettings) -> float:
    if r_nm >= settings.cutoff_radius_nm:
        return 0.0
    if r_nm <= settings.switch_radius_nm:
        return 1.0
    x = (r_nm - settings.switch_radius_nm) / (settings.cutoff_radius_nm - settings.switch_radius_nm)
    return 1.0 - (3.0 * x * x - 2.0 * x * x * x)


def _combined_lj_params(p1: MMAtomParams, p2: MMAtomParams) -> Tuple[float, float]:
    sigma = 0.5 * (p1.sigma + p2.sigma)
    epsilon = math.sqrt(max(0.0, p1.epsilon * p2.epsilon))
    return sigma, epsilon


def _lj_energy_kj(r_nm: float, p1: MMAtomParams, p2: MMAtomParams) -> float:
    sigma, epsilon = _combined_lj_params(p1, p2)
    if sigma <= 0.0 or epsilon <= 0.0 or r_nm <= 1.0e-12:
        return 0.0
    sr6 = (sigma / r_nm) ** 6
    return 4.0 * epsilon * (sr6 * sr6 - sr6)


def _coulomb_energy_kj(r_nm: float, p1: MMAtomParams, p2: MMAtomParams) -> float:
    if r_nm <= 1.0e-12 or p1.charge == 0.0 or p2.charge == 0.0:
        return 0.0
    return COULOMB_KJ_MOL_NM * p1.charge * p2.charge / r_nm


def _free_rotor_group_atom_ids(mol: Molecule, anchor_step: PlannedAtomBuildStep) -> List[AtomID]:
    residue = get_residue(mol, anchor_step.atom_key)
    ids: List[AtomID] = []
    for atom_name in anchor_step.free_rotor_group_atoms or (anchor_step.atom_key.atom_name,):
        if atom_name in residue.atoms:
            ids.append(AtomID(anchor_step.atom_key.chain_id, anchor_step.atom_key.residue_index, atom_name))
    if anchor_step.atom_key not in ids:
        ids.insert(0, anchor_step.atom_key)
    return ids


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
) -> Tuple[Set[AtomID], Dict[AtomID, float]]:
    """Precompute all free-rotor H charges collapsed onto parent atoms.

    Free-rotor hydrogens are omitted from the spatial environment.  Their
    charges are represented only when the corresponding parent heavy atom is
    selected into a local environment; hydrogen LJ parameters are discarded.
    """
    collapsed_hydrogens = _template_free_rotor_hydrogen_ids(mol, template)
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
    if settings.collapse_other_free_rotors:
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
    """Return exact 1-2, 1-3, and 1-4 atom shells around start."""
    visited = {start}
    frontier = {start}
    shells: List[Set[AtomID]] = []
    for _depth in range(3):
        next_frontier: Set[AtomID] = set()
        for atom_id in frontier:
            next_frontier.update(graph.get(atom_id, set()))
        next_frontier.difference_update(visited)
        shells.append(next_frontier)
        visited.update(next_frontier)
        frontier = next_frontier
    return shells[0], shells[1], shells[2]


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
            for term in terms:
                energy += term.k * (1.0 + math.cos(math.radians(term.periodicity * phi - term.phase)))
    return energy


def _nonbonded_energy_for_trial_group(
    mm_params: MMParameterProvider,
    trial_atoms: Sequence[_TrialGroupAtom],
    pair_lists: _FreeRotorPairLists,
    settings: FreeRotorSearchSettings,
) -> float:
    energy = 0.0
    cutoff_squared = settings.cutoff_radius_nm ** 2
    switch_radius = settings.switch_radius_nm
    switch_width = settings.cutoff_radius_nm - switch_radius
    for trial in trial_atoms:
        pairs = pair_lists.pairs_by_hydrogen.get(trial.atom_id)
        if pairs is None or pairs.coords_nm.shape[0] == 0:
            continue

        delta = pairs.coords_nm - trial.coord_nm
        r_squared = np.einsum("ij,ij->i", delta, delta)
        active = (r_squared < cutoff_squared) & (r_squared > 1.0e-24)
        if not np.any(active):
            continue

        r_nm = np.sqrt(r_squared[active])
        weights = np.ones_like(r_nm)
        switching = r_nm > switch_radius
        if np.any(switching):
            x = (r_nm[switching] - switch_radius) / switch_width
            weights[switching] = 1.0 - (3.0 * x * x - 2.0 * x * x * x)

        pair_energy = np.zeros_like(r_nm)
        if settings.include_lj:
            sigma = 0.5 * (trial.params.sigma + pairs.sigmas[active])
            epsilon = np.sqrt(
                np.maximum(0.0, trial.params.epsilon * pairs.epsilons[active])
            )
            valid_lj = (sigma > 0.0) & (epsilon > 0.0)
            if np.any(valid_lj):
                sr6 = (sigma[valid_lj] / r_nm[valid_lj]) ** 6
                pair_energy[valid_lj] += (
                    pairs.scale_lj[active][valid_lj]
                    * 4.0
                    * epsilon[valid_lj]
                    * (sr6 * sr6 - sr6)
                )
        if settings.include_electrostatics and trial.params.charge != 0.0:
            pair_energy += (
                pairs.scale_qq[active]
                * COULOMB_KJ_MOL_NM
                * trial.params.charge
                * pairs.charges[active]
                / r_nm
            )
        energy += float(np.sum(weights * pair_energy))
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


def _parabolic_grid_minimum(center_phi: float, e_prev: float, e_center: float, e_next: float, step: float) -> float:
    denom = e_prev - 2.0 * e_center + e_next
    if denom <= 1.0e-12 or not math.isfinite(denom):
        return center_phi
    delta = 0.5 * step * (e_prev - e_next) / denom
    if not math.isfinite(delta) or abs(delta) > step:
        return center_phi
    return center_phi + delta


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

    group_size = max(1, len(step.free_rotor_group_atoms or (step.atom_key.atom_name,)))
    search_period = 360.0 / group_size
    n_grid = max(3, int(round(search_period / settings.grid_step_deg)))
    step_deg = search_period / n_grid

    phis_unwrapped = [template_dihedral + i * step_deg for i in range(n_grid)]
    energies = [
        score_free_rotor_group(
            mol,
            template,
            step,
            rule,
            normalize_angle(phi),
            mm_params,
            mm_cache,
            pair_lists,
            settings,
        )
        for phi in phis_unwrapped
    ]

    min_idx = min(range(n_grid), key=lambda i: (energies[i], abs(normalize_angle(phis_unwrapped[i] - template_dihedral))))
    prev_idx = (min_idx - 1) % n_grid
    next_idx = (min_idx + 1) % n_grid
    center = phis_unwrapped[min_idx]
    refined = _parabolic_grid_minimum(center, energies[prev_idx], energies[min_idx], energies[next_idx], step_deg)
    return normalize_angle(refined)


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


def _build_dihedral_for_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    rule: Mapping[str, Any],
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
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
        raise RuntimeError(
            "Cannot execute supplied_dof step without an external DOF value. "
            "execute_build_plan_until_missing_dof should stop before such steps."
        )

    raise RuntimeError(f"Unsupported torsion_source {step.torsion_source!r} for {step.atom_key}")


def build_atom_from_step(
    mol: Molecule,
    template: Mapping[str, Any],
    step: PlannedAtomBuildStep,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
    mm_cache: Optional[FreeRotorMMCache] = None,
) -> None:
    """Execute one PlannedAtomBuildStep in-place on `mol`."""
    residue = get_residue(mol, step.atom_key)
    atom = residue.atoms[step.atom_key.atom_name]
    if atom.coord is not None:
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


def execute_build_plan_until_missing_dof(
    molecule: Molecule,
    template: Mapping[str, Any],
    build_plan: BuildPlan,
    *,
    modify_myself: bool = False,
    mm_params: Optional[MMParameterProvider] = None,
    free_rotor_settings: Optional[FreeRotorSearchSettings] = None,
) -> Tuple[Molecule, BuildPlan]:
    """Execute atom-build steps until the first missing-DOF step or plan end.

    If modify_myself is False, both molecule and build_plan are deep-copied and
    the returned objects are independent trial outputs. If True, the input
    objects are modified in-place and also returned.
    """
    mol = molecule if modify_myself else copy.deepcopy(molecule)
    plan = build_plan if modify_myself else copy.deepcopy(build_plan)

    mm_cache = (
        FreeRotorMMCache(build_mm_bond_graph(mol, mm_params))
        if mm_params is not None
        else None
    )
    executed_count = 0
    for step in plan.steps:
        if isinstance(step, PlannedMissingDOFStep):
            break
        if not isinstance(step, PlannedAtomBuildStep):
            raise RuntimeError(f"Unknown plan step type: {step!r}")
        build_atom_from_step(
            mol,
            template,
            step,
            mm_params,
            free_rotor_settings,
            mm_cache,
        )
        executed_count += 1
    if executed_count:
        plan.steps = plan.steps[executed_count:]

    # Remaining requirement list should describe only missing-DOF steps still in
    # the unexecuted plan prefix/suffix.
    plan.requirements = [s for s in plan.steps if isinstance(s, PlannedMissingDOFStep)]
    plan.unresolved_atoms = collect_pending_atoms(mol)
    return mol, plan


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
    else:
        lines.append("Required missing DOFs: none")
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

    out_steps: List[Dict[str, Any]] = []
    for step in plan.steps:
        if isinstance(step, PlannedMissingDOFStep):
            out_steps.append({
                "type": "missing_dof",
                "dof_key": {
                    "chain_id": step.dof_key.chain_id,
                    "residue_index": step.dof_key.residue_index,
                    "atom": step.dof_key.atom_name,
                    "rule_index": step.dof_key.rule_index,
                },
                "central_bond": [ref_dict(r) for r in step.central_bond],
                "requested_dihedral_atoms": [ref_dict(r) for r in step.requested_dihedral_atoms],
                "torsion_group_index": step.torsion_group_index,
                "requested_member_index": step.requested_member_index,
                "reason_atom": atom_id_dict(step.reason_atom),
                "reason_rule_index": step.reason_rule_index,
            })
        else:
            out_steps.append({
                "type": "build_atom",
                "atom": atom_id_dict(step.atom_key),
                "rule_index": step.rule_index,
                "torsion_class": step.torsion_class,
                "torsion_source": step.torsion_source,
                "observed_member_index": step.observed_member_index,
                "dof_key": None if step.dof_key is None else {
                    "chain_id": step.dof_key.chain_id,
                    "residue_index": step.dof_key.residue_index,
                    "atom": step.dof_key.atom_name,
                    "rule_index": step.dof_key.rule_index,
                },
                "free_rotor_anchor": None if step.free_rotor_anchor is None else atom_id_dict(step.free_rotor_anchor),
                "phase_offset": step.phase_offset,
                "free_rotor_group_atoms": list(step.free_rotor_group_atoms),
            })
    return {"summary": summarize_build_plan(plan), "steps": out_steps}
