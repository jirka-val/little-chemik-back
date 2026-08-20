"""Public end-to-end integration API for the FORGE structure-treatment workflow.

The caller is responsible for the upstream cleaning contract documented in
``INTEGRATION_CONTRACT.md``.  In particular, this module deliberately does not
apply the development-only HIS/alias/residue filtering shim.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from forge_molecule_builder import (
    BuildPlan,
    execute_build_plan_until_missing_dof,
    plan_build_steps,
)
from forge_molecule_mm import FreeRotorSearchSettings
from forge_molecule_ions import (
    CrystalIonCleanupReport,
    IonAdditionReport,
    IonPlacementSettings,
    SaltSpecification,
    add_ions_to_solvated_molecule,
    clean_crystal_ions,
)
from forge_molecule_parser import Molecule, build_molecule_from_forge_json, load_json
from forge_molecule_solvation import (
    SolvationReport,
    SolvationSettings,
    SolvationVdwParameters,
    WaterSolvationTemplate,
    load_solvation_template,
    solvate_molecule,
)
from forge_molecule_state_assignment import (
    HydrogenBondGeometrySettings,
    StateAssignmentReport,
    assign_molecule_states,
)


@dataclass(frozen=True)
class WorkflowResources:
    """Validated runtime libraries and force-field parameters."""

    converting_dictionary: Mapping[str, Any]
    building_template: Mapping[str, Any]
    state_definitions: Mapping[str, Any]
    water_template: WaterSolvationTemplate
    force_field_parameters: SolvationVdwParameters

    @classmethod
    def from_paths(
        cls,
        *,
        converting_dictionary: str | Path,
        building_template: str | Path,
        state_definitions: str | Path,
        solvation_template: str | Path,
        force_field_root: str | Path,
    ) -> "WorkflowResources":
        """Load the complete runtime resource set once at process startup."""

        return cls(
            converting_dictionary=load_json(converting_dictionary),
            building_template=load_json(building_template),
            state_definitions=load_json(state_definitions),
            water_template=load_solvation_template(solvation_template),
            force_field_parameters=SolvationVdwParameters.from_force_field_root(
                force_field_root
            ),
        )


@dataclass
class WorkflowSettings:
    """Runtime choices for one structure-treatment invocation."""

    pH: float = 7.0
    covalent_cutoff_angstrom: float = 2.3
    hydrogen_bond: HydrogenBondGeometrySettings = field(
        default_factory=HydrogenBondGeometrySettings
    )
    free_rotor: FreeRotorSearchSettings = field(
        default_factory=FreeRotorSearchSettings
    )
    solvation: SolvationSettings = field(default_factory=SolvationSettings)
    ions: IonPlacementSettings = field(default_factory=IonPlacementSettings)
    add_solvent_and_ions: bool = True


@dataclass
class WorkflowResult:
    """Structured result; no user-facing text parsing is required."""

    status: str
    molecule: Molecule
    remaining_plan: BuildPlan
    state_assignment: StateAssignmentReport
    crystal_ion_cleanup: Optional[CrystalIonCleanupReport] = None
    solvation: Optional[SolvationReport] = None
    ion_addition: Optional[IonAdditionReport] = None

    @property
    def completed(self) -> bool:
        return self.status == "complete"

    @property
    def stopped_at_missing_dof(self) -> bool:
        return self.status == "missing_dof"


def run_forge_workflow(
    structure_data: Mapping[str, Any],
    resources: WorkflowResources,
    *,
    salts: Sequence[SaltSpecification] = (),
    settings: Optional[WorkflowSettings] = None,
) -> WorkflowResult:
    """Parse, assign states, build, solvate and add ions.

    A missing DOF is a normal boundary, not an exception.  The function then
    returns the partially built molecule and the remaining plan without
    cleaning ions, moving coordinates, solvating, or adding ions.  Production
    code must explicitly resolve the DOF and resume through a future API; it
    must never substitute the benchmark-only template-dihedral bridge.
    """

    settings = settings or WorkflowSettings()
    molecule = build_molecule_from_forge_json(
        structure_data, resources.converting_dictionary
    )
    molecule, state_report = assign_molecule_states(
        molecule,
        resources.converting_dictionary,
        resources.building_template,
        resources.state_definitions,
        pH=settings.pH,
        covalent_cutoff_angstrom=settings.covalent_cutoff_angstrom,
        hydrogen_bond_settings=settings.hydrogen_bond,
        modify_myself=True,
    )

    plan = plan_build_steps(molecule, resources.building_template)
    molecule, remaining = execute_build_plan_until_missing_dof(
        molecule,
        resources.building_template,
        plan,
        modify_myself=True,
        mm_params=resources.force_field_parameters.merged_builder_provider(),
        free_rotor_settings=settings.free_rotor,
    )
    if remaining.steps:
        return WorkflowResult(
            status="missing_dof",
            molecule=molecule,
            remaining_plan=remaining,
            state_assignment=state_report,
        )

    if not settings.add_solvent_and_ions:
        return WorkflowResult(
            status="complete",
            molecule=molecule,
            remaining_plan=remaining,
            state_assignment=state_report,
        )

    molecule, cleanup_report = clean_crystal_ions(
        molecule,
        resources.building_template,
        resources.force_field_parameters,
        settings.ions,
        modify_myself=True,
    )
    molecule, solvation_report = solvate_molecule(
        molecule,
        resources.water_template,
        resources.force_field_parameters,
        settings.solvation,
        modify_myself=True,
    )
    molecule, ion_report = add_ions_to_solvated_molecule(
        molecule,
        salts,
        resources.force_field_parameters,
        settings.ions,
        modify_myself=True,
    )
    return WorkflowResult(
        status="complete",
        molecule=molecule,
        remaining_plan=remaining,
        state_assignment=state_report,
        crystal_ion_cleanup=cleanup_report,
        solvation=solvation_report,
        ion_addition=ion_report,
    )
