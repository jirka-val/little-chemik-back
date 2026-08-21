# app/api/v1/endpoints/sidechains.py
"""
Interaktivní dostavování bezpečných side-chain větví (missing_dof z
/api/validation/prepare, které builder umí bezpečně dostavět přes GUI slidery
místo tvrdého 409 dead-endu) - viz app/builder/INTEGRATION_CONTRACT.md a
app/services/structure/sidechain_service.py.

start -> libovolně mnoho update/optimize -> commit (nebo cancel). Mezi
jednotlivými voláními drží stav sidechain_session_service (in-memory,
per-workspace) - "start" ho založí, "commit"/"cancel" ho zahodí.
"""
import json
import logging
import time
from typing import Any, Dict, List

import aiofiles
from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.api.v1.endpoints.validation import PreparationRequest, _BOX_SHAPE_MAP, _build_salt_specs
from app.core.exceptions import AppBaseException, InternalError
from app.services.structure.sidechain_service import sidechain_session_service
from app.services.validation.service import ValidationService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)
router = APIRouter()
validation_service = ValidationService()


class SidechainUpdateRequest(BaseModel):
    sidechain_data: Dict[str, Any] = Field(..., description="Aktuální GUI payload (viz gui_payload z /start)")
    changed_dof_keys: List[Dict[str, Any]] = Field(
        ..., description="DOF klíče změněné tímto sliderovým tahem, tvar {chain_id, residue_index, atom, rule_index}"
    )


class SidechainOptimizeRequest(BaseModel):
    sidechain_data: Dict[str, Any] = Field(...)
    residue: Dict[str, Any] = Field(..., description="{chain_id, residue_index} vybraného rezidua pro Opt")


@router.post("/start/{workspace_id}", summary="Spustí přípravu; otevře GUI relaci, pokud narazí na bezpečný missing DOF")
async def start_sidechain_session(workspace_id: str, request: PreparationRequest):
    """
    Stejný vstup jako /api/validation/prepare. Pokud struktura žádnou
    interaktivní volbu nepotřebuje, chová se identicky (rovnou hotovo).
    Jinak založí side-chain relaci a vrátí gui_payload pro slidery + jméno
    PDB souboru s počátečním FF-optimálním náhledem (pro první plné načtení
    do Mol* - další updaty už jdou přes malé coordinate patche z /update
    a /optimize, ne přes další plné stažení).
    """
    if workspace_id != request.workspace_id:
        raise InternalError("workspace_id in path and body must match.")

    start_time = time.time()
    workspace_manager.require_workspace(workspace_id)

    try:
        pdb_path = workspace_manager.get_file_path(workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        outcome = await run_in_threadpool(
            sidechain_session_service.start,
            workspace_id=workspace_id,
            pdb_text=pdb_content,
            ff_selections=request.ff_selections,
            ph=request.ph,
            add_solvent_and_ions=request.add_solvent,
            salts=_build_salt_specs(request.positive_ion, request.negative_ion, request.ionic_strength),
            box_shape=_BOX_SHAPE_MAP.get(request.box_shape),
            box_padding_angstrom=request.box_padding_nm * 10.0,
            keep_crystal_waters=request.crystal_water_mode != "remove_all",
            crystal_water_mode=request.crystal_water_mode,
        )

        if outcome.status == "complete":
            prepared = outcome.prepared
            meta_path = pdb_path.with_name(pdb_path.name.replace(".pdb", ".forge_meta.json"))
            async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
                await f.write(prepared.pdb_text)
            async with aiofiles.open(meta_path, "w", encoding="utf-8") as f:
                await f.write(json.dumps(prepared.forge_meta))

            validation_results = await run_in_threadpool(
                validation_service.validate_pdb_content, prepared.pdb_text, label="prepared_state"
            )
            logger.info(f"Sidechain start: structure already complete for {workspace_id} in {time.time() - start_time:.2f}s.")
            return {
                "status": "complete",
                "message": "Structure successfully prepared.",
                "warnings": prepared.warnings,
                "validation": validation_results,
            }

        preview_path = workspace_manager.get_file_path(workspace_id, outcome.preview_filename)
        async with aiofiles.open(preview_path, "w", encoding="utf-8") as f:
            await f.write(outcome.preview_pdb_text)

        logger.info(
            f"Sidechain start: opened GUI session for {workspace_id} "
            f"({len(outcome.gui_payload.get('sidechains', []))} residues) in {time.time() - start_time:.2f}s."
        )
        return {
            "status": "missing_dof",
            "gui_payload": outcome.gui_payload,
            "preview_filename": outcome.preview_filename,
        }

    except AppBaseException:
        raise
    except Exception as e:
        logger.exception(f"Sidechain start failed for workspace {workspace_id}: {e}")
        raise InternalError(str(e))


async def _refresh_preview_file(workspace_id: str) -> None:
    """
    structure_preview.pdb je jen sidecar pro GUI náhled (viz /start) - working_
    molecule v paměti session se aktualizuje při každém update()/optimize(),
    ale soubor na disku ne, dokud ho sem výslovně nezapíšeme. Dokud
    viewer.applyCoordinatePatch() na frontendu dělá tichý reload místo
    skutečného coordinate patche, bez tohohle by po prvním posunu slideru
    zůstalo GUI vizuálně na počátečním stavu.
    """
    pdb_text = await run_in_threadpool(sidechain_session_service.preview_pdb_text, workspace_id)
    preview_path = workspace_manager.get_file_path(workspace_id, "structure_preview.pdb")
    async with aiofiles.open(preview_path, "w", encoding="utf-8") as f:
        await f.write(pdb_text)


@router.post("/update/{workspace_id}", summary="Přepočte atomy ovlivněné změněnými DOF a vrátí jen jejich souřadnice")
async def update_sidechains(workspace_id: str, request: SidechainUpdateRequest):
    try:
        patch = await run_in_threadpool(
            sidechain_session_service.update,
            workspace_id,
            request.sidechain_data,
            request.changed_dof_keys,
        )
        await _refresh_preview_file(workspace_id)
        return patch
    except AppBaseException:
        raise
    except Exception as e:
        logger.exception(f"Sidechain update failed for workspace {workspace_id}: {e}")
        raise InternalError(str(e))


@router.post("/optimize/{workspace_id}", summary="Lokální MM optimalizace ('Opt') jednoho vybraného rezidua")
async def optimize_sidechain(workspace_id: str, request: SidechainOptimizeRequest):
    try:
        response = await run_in_threadpool(
            sidechain_session_service.optimize,
            workspace_id,
            request.sidechain_data,
            request.residue,
        )
        await _refresh_preview_file(workspace_id)
        return response
    except AppBaseException:
        raise
    except Exception as e:
        logger.exception(f"Sidechain optimize failed for workspace {workspace_id}: {e}")
        raise InternalError(str(e))


@router.post("/commit/{workspace_id}", summary="Přijme aktuální GUI hodnoty a dokončí zbytek přípravy (solvatace, ionty)")
async def commit_sidechains(workspace_id: str):
    start_time = time.time()
    try:
        prepared = await run_in_threadpool(sidechain_session_service.commit, workspace_id)

        pdb_path = workspace_manager.get_file_path(workspace_id)
        meta_path = pdb_path.with_name(pdb_path.name.replace(".pdb", ".forge_meta.json"))
        async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
            await f.write(prepared.pdb_text)
        async with aiofiles.open(meta_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(prepared.forge_meta))

        preview_path = workspace_manager.get_file_path(workspace_id, "structure_preview.pdb")
        if preview_path.exists():
            preview_path.unlink()

        validation_results = await run_in_threadpool(
            validation_service.validate_pdb_content, prepared.pdb_text, label="prepared_state"
        )
        logger.info(f"Sidechain commit: finished for {workspace_id} in {time.time() - start_time:.2f}s.")
        return {
            "message": "Structure successfully prepared.",
            "warnings": prepared.warnings,
            "validation": validation_results,
        }
    except AppBaseException:
        # ForgeMissingDOFError (další, ne-bezpečný missing DOF) apod. - stejný
        # 409 kontrakt jako /api/validation/prepare, session zůstává zachovaná.
        raise
    except Exception as e:
        logger.exception(f"Sidechain commit failed for workspace {workspace_id}: {e}")
        raise InternalError(str(e))


@router.post("/cancel/{workspace_id}", summary="Zahodí GUI relaci beze změny structure.pdb")
async def cancel_sidechains(workspace_id: str):
    sidechain_session_service.cancel(workspace_id)
    return {"message": "Side-chain session cancelled."}
