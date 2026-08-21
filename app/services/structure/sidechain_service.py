"""
Interaktivní dostavování bezpečných "residue-local open branch" side-chainů
(chybějící báze/postranní řetězce s jedním kotevním bodem) - viz
app/builder/INTEGRATION_CONTRACT.md, sekce "Interactive missing-side-chain API".

Na rozdíl od zbytku backendu (workspace_manager, editor_service, ...), který je
čistě souborový (přečti PDB, uprav, zapiš), tohle drží živé Python objekty
(Molecule/BuildPlan/MM parametry) mezi jednotlivými HTTP requesty jedné
GUI relace (start -> libovolně mnoho update/optimize -> commit nebo cancel).
Backend běží jako jeden proces (stejný předpoklad jako workspace_manager a
_cached_ff_parameters/_static_resources v forge_service.py), takže jednoduchý
in-memory dict klíčovaný workspace_id stačí - žádná perzistence napříč restarty
není potřeba, session je vždy jen jedna aktivní GUI relace na workspace.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

logger = logging.getLogger(__name__)

_BUILDER_DIR = Path(__file__).resolve().parents[2] / "builder"
if str(_BUILDER_DIR) not in sys.path:
    sys.path.insert(0, str(_BUILDER_DIR))

from forge_molecule_builder import (  # noqa: E402
    BuildPlan,
    DOFKey,
    ResidueID,
    SidechainExecutionIndex,
    SidechainLocalOptimizationContext,
    build_missing_sidechain_gui_payload,
    complete_missing_sidechains,
    execute_build_plan_until_missing_dof,
    optimize_missing_sidechain_dofs,
    optimize_missing_sidechain_preview,
    prepare_missing_sidechain_local_optimization,
    prepare_sidechain_execution_index,
    update_missing_sidechains,
)
from forge_molecule_ions import add_ions_to_solvated_molecule, clean_crystal_ions  # noqa: E402
from forge_molecule_mm import MMParameterProvider  # noqa: E402
from forge_molecule_parser import Molecule  # noqa: E402
from forge_molecule_solvation import solvate_molecule  # noqa: E402
from forge_workflow import WorkflowResources, WorkflowSettings  # noqa: E402

from app.core.exceptions import NotFoundError
from app.services.structure.forge_service import (
    ForgeMissingDOFError,
    ForgePreparationResult,
    ForgeStructureService,
    build_forge_meta,
    molecule_to_pdb,
)


class SidechainSessionNotFoundError(NotFoundError):
    """Žádná aktivní side-chain relace pro tenhle workspace (nikdy nezačala, nebo už byla commit/cancel)."""

    def __init__(self):
        super().__init__(
            "No active side-chain completion session for this workspace. Call /sidechains/start first."
        )


@dataclass
class SidechainStartResult:
    status: str  # "complete" | "missing_dof"
    prepared: Optional[ForgePreparationResult] = None
    gui_payload: Optional[Dict[str, Any]] = None
    preview_filename: Optional[str] = None
    preview_pdb_text: Optional[str] = None


@dataclass
class SidechainSession:
    resources: WorkflowResources
    settings: WorkflowSettings
    salts: List[Any]
    base_molecule: Molecule
    remaining_plan: BuildPlan
    state_assignment: Any
    execution_index: SidechainExecutionIndex
    mm_parameters: MMParameterProvider
    sidechain_data: Dict[str, Any]
    working_molecule: Molecule
    local_opt_contexts: Dict[ResidueID, SidechainLocalOptimizationContext] = field(default_factory=dict)


def _dof_key_from_payload(data: Mapping[str, Any]) -> DOFKey:
    """
    Vlastní odpovídač GUI DOF payloadu na DOFKey - zrcadlí tvar, jaký
    build_missing_sidechain_gui_payload() sám generuje (viz
    forge_molecule_builder._dof_key_payload), knihovna ale svoji verzi téhle
    funkce nevystavuje jako veřejné API.
    """
    return DOFKey(
        chain_id=str(data["chain_id"]),
        residue_index=int(data["residue_index"]),
        atom_name=str(data["atom"]),
        rule_index=int(data["rule_index"]),
    )


def _dof_values_from_sidechain_data(sidechain_data: Mapping[str, Any]) -> Dict[DOFKey, float]:
    """
    Odvodí kompletní {DOFKey: hodnota ve stupních} z aktuálního GUI JSON
    payloadu. Sidechain_data je jediný zdroj pravdy držený v session - tohle
    se počítá on-demand (hlavně pro commit()), místo aby se muselo ručně
    udržovat v synchronizaci při každém update/optimize volání.
    """
    values: Dict[DOFKey, float] = {}
    for sidechain in sidechain_data.get("sidechains", []):
        for dof in sidechain.get("dofs", []):
            key = _dof_key_from_payload(dof["dof_key"])
            values[key] = float(dof["value_degrees"])
    return values


def _residue_id_from_payload(data: Mapping[str, Any]) -> ResidueID:
    return ResidueID(chain_id=str(data["chain_id"]), residue_index=int(data["residue_index"]))


class SidechainSessionService:
    """Drží jednu aktivní interaktivní side-chain relaci na workspace_id."""

    def __init__(self, forge_service: ForgeStructureService):
        self.forge_service = forge_service
        self._sessions: Dict[str, SidechainSession] = {}

    def _require_session(self, workspace_id: str) -> SidechainSession:
        session = self._sessions.get(workspace_id)
        if session is None:
            raise SidechainSessionNotFoundError()
        return session

    def start(
        self,
        workspace_id: str,
        pdb_text: str,
        ff_selections: Dict[str, Any],
        ph: float = 7.0,
        add_solvent_and_ions: bool = True,
        salts: Optional[List[Dict[str, Any]]] = None,
        box_shape: Optional[str] = None,
        box_padding_angstrom: Optional[float] = None,
        keep_crystal_waters: Optional[bool] = None,
        crystal_water_mode: str = "remove_all",
    ) -> SidechainStartResult:
        run = self.forge_service.run_workflow(
            pdb_text,
            ff_selections,
            ph=ph,
            add_solvent_and_ions=add_solvent_and_ions,
            salts=salts,
            box_shape=box_shape,
            box_padding_angstrom=box_padding_angstrom,
            keep_crystal_waters=keep_crystal_waters,
            crystal_water_mode=crystal_water_mode,
        )
        result = run.result

        if not result.stopped_at_missing_dof:
            prepared = ForgePreparationResult(
                pdb_text=molecule_to_pdb(result.molecule),
                forge_meta=build_forge_meta(result.molecule),
                warnings=list(result.molecule.warnings),
                state_assignment=result.state_assignment,
                crystal_ion_cleanup=result.crystal_ion_cleanup,
                solvation=result.solvation,
                ion_addition=result.ion_addition,
            )
            self._sessions.pop(workspace_id, None)
            return SidechainStartResult(status="complete", prepared=prepared)

        template = run.resources.building_template
        mm_parameters = run.resources.force_field_parameters.merged_builder_provider()
        execution_index = prepare_sidechain_execution_index(result.remaining_plan)
        if not execution_index.dof_steps:
            # První missing DOF vůbec není bezpečná residue-local větev (multi-
            # anchor/cyklický/chybějící polymer span - viz "Known v1 limits" v
            # INTEGRATION_CONTRACT.md) - GUI by tu neměla co nabídnout (prázdný
            # panel), takže se chováme přesně jako dřív: stejný 409 kontrakt
            # jako prepare_structure().
            raise ForgeMissingDOFError(result.remaining_plan.steps[0], result.molecule)

        default_dofs = optimize_missing_sidechain_dofs(
            result.molecule,
            template,
            result.remaining_plan,
            mm_parameters,
            execution_index=execution_index,
        )
        gui_payload = build_missing_sidechain_gui_payload(
            result.molecule,
            result.remaining_plan,
            default_dofs,
            execution_index=execution_index,
        )
        working_molecule = complete_missing_sidechains(
            result.molecule,
            template,
            result.remaining_plan,
            gui_payload,
            execution_index=execution_index,
            mm_params=mm_parameters,
        )

        self._sessions[workspace_id] = SidechainSession(
            resources=run.resources,
            settings=run.settings,
            salts=run.salts,
            base_molecule=result.molecule,
            remaining_plan=result.remaining_plan,
            state_assignment=result.state_assignment,
            execution_index=execution_index,
            mm_parameters=mm_parameters,
            sidechain_data=gui_payload,
            working_molecule=working_molecule,
        )
        return SidechainStartResult(
            status="missing_dof",
            gui_payload=gui_payload,
            preview_filename="structure_preview.pdb",
            preview_pdb_text=molecule_to_pdb(working_molecule),
        )

    def update(
        self,
        workspace_id: str,
        sidechain_data: Dict[str, Any],
        changed_dof_keys: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        session = self._require_session(workspace_id)
        changed = [_dof_key_from_payload(entry) for entry in changed_dof_keys]

        session.working_molecule, patch = update_missing_sidechains(
            session.working_molecule,
            session.resources.building_template,
            session.remaining_plan,
            sidechain_data,
            changed,
            execution_index=session.execution_index,
            mm_params=session.mm_parameters,
        )
        session.sidechain_data = sidechain_data
        touched = {ResidueID(key.chain_id, key.residue_index) for key in changed}
        self._invalidate_contexts_except(session, touched)
        return patch

    def optimize(
        self,
        workspace_id: str,
        sidechain_data: Dict[str, Any],
        residue: Dict[str, Any],
    ) -> Dict[str, Any]:
        session = self._require_session(workspace_id)
        residue_id = _residue_id_from_payload(residue)
        context = session.local_opt_contexts.get(residue_id)
        if context is None:
            context = prepare_missing_sidechain_local_optimization(
                session.working_molecule,
                session.resources.building_template,
                session.remaining_plan,
                residue_id,
                session.mm_parameters,
                execution_index=session.execution_index,
            )
            session.local_opt_contexts[residue_id] = context

        session.working_molecule, response = optimize_missing_sidechain_preview(
            session.working_molecule,
            session.resources.building_template,
            session.remaining_plan,
            sidechain_data,
            residue_id,
            session.mm_parameters,
            execution_index=session.execution_index,
            optimization_context=context,
        )
        session.sidechain_data = sidechain_data
        self._invalidate_contexts_except(session, {residue_id})
        return response

    def _invalidate_contexts_except(
        self,
        session: SidechainSession,
        touched_residue_ids: set,
    ) -> None:
        """
        SidechainLocalOptimizationContext zůstává platný, dokud se mění jen DOF
        v jeho vlastním reziduu (viz docstring SidechainLocalOptimizationContext
        v forge_molecule_builder.py) - musí se zahodit, jakmile se změní
        souřadnice MIMO vybrané reziduum. Zahodíme tedy cache pro každé
        reziduum, které NENÍ mezi právě změněnými.
        """
        for residue_id in list(session.local_opt_contexts):
            if residue_id not in touched_residue_ids:
                session.local_opt_contexts.pop(residue_id, None)

    def commit(self, workspace_id: str) -> ForgePreparationResult:
        session = self._require_session(workspace_id)
        template = session.resources.building_template
        dof_values = _dof_values_from_sidechain_data(session.sidechain_data)

        molecule, plan = execute_build_plan_until_missing_dof(
            session.base_molecule,
            template,
            session.remaining_plan,
            dof_values=dof_values,
            modify_myself=False,
            mm_params=session.mm_parameters,
            free_rotor_settings=session.settings.free_rotor,
        )

        if plan.steps:
            # Další missing DOF, který nebyl mezi bezpečnými side-chain
            # větvemi (multi-anchor/cyklický/chybějící polymer span) - stejný
            # 409 kontrakt jako u prepare_structure(), session zůstává (uživatel
            # o svoje rozhodnutí nepřijde, může to zkusit jinak nebo zrušit).
            raise ForgeMissingDOFError(plan.steps[0], molecule)

        cleanup_report = solvation_report = ion_report = None
        if session.settings.add_solvent_and_ions:
            # Zrcadlí ocas run_forge_workflow() (viz forge_workflow.py) - ten
            # samý krok se tam přeskakuje, když settings.add_solvent_and_ions
            # je False, takže stejnou podmínku musíme mít i tady, jinak by
            # /sidechains/commit dělal solvataci/ionty i tehdy, kdy o ně
            # /sidechains/start vůbec nepožádal.
            molecule, cleanup_report = clean_crystal_ions(
                molecule,
                template,
                session.resources.force_field_parameters,
                session.settings.ions,
                modify_myself=True,
            )
            molecule, solvation_report = solvate_molecule(
                molecule,
                session.resources.water_template,
                session.resources.force_field_parameters,
                session.settings.solvation,
                modify_myself=True,
            )
            molecule, ion_report = add_ions_to_solvated_molecule(
                molecule,
                session.salts,
                session.resources.force_field_parameters,
                session.settings.ions,
                modify_myself=True,
            )

        self._sessions.pop(workspace_id, None)
        return ForgePreparationResult(
            pdb_text=molecule_to_pdb(molecule),
            forge_meta=build_forge_meta(molecule),
            warnings=list(molecule.warnings),
            state_assignment=session.state_assignment,
            crystal_ion_cleanup=cleanup_report,
            solvation=solvation_report,
            ion_addition=ion_report,
        )

    def cancel(self, workspace_id: str) -> None:
        self._sessions.pop(workspace_id, None)

    def preview_pdb_text(self, workspace_id: str) -> str:
        """
        Aktuální stav working_molecule jako PDB text - endpoint tohle po
        každém update()/optimize() volání znovu zapisuje do
        structure_preview.pdb, protože ten se (na rozdíl od working_molecule
        v paměti) sám od sebe neaktualizuje. Dokud viewer.applyCoordinatePatch()
        dělá tichý plný reload místo skutečného coordinate patche (viz jeho
        docstring), bez tohohle by GUI po prvním posunu slideru ukazovalo
        zastaralou (počáteční) geometrii.
        """
        session = self._require_session(workspace_id)
        return molecule_to_pdb(session.working_molecule)


sidechain_session_service = SidechainSessionService(ForgeStructureService())
