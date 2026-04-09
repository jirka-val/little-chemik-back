import logging
import traceback
import time
import aiofiles
from typing import Dict, Literal

from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.services.validation.service import ValidationService
from app.services.structure.hydrogenation import HydrogenationService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
validation_service = ValidationService()
hydrogenation_service = HydrogenationService()


# --- Data Models ---

class ValidationRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    label: str = Field("molecule_from_front", description="Molecule state identifier")


class FixAltLocRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    selections: Dict[str, str] = Field(
        ...,
        description="Variant selection map, e.g., {'A-42': 'B'}",
        example={"A-42": "B", "A-15": "A"}
    )


class PreparationRequest(BaseModel):
    workspace_id: str = Field(...)
    ph: float = Field(7.0)
    crystal_water_mode: Literal["remove_all", "keep_water", "keep_all"] = Field("remove_all")
    add_solvent: bool = Field(False)
    box_padding_nm: float = Field(1.0)
    ionic_strength: float = Field(0.15, description="Salt concentration (M). The system is automatically neutralized.")
    positive_ion: Literal["Na+", "K+", "Li+", "Cs+", "Rb+"] = Field("Na+")
    negative_ion: Literal["Cl-", "F-", "Br-", "I-"] = Field("Cl-")


# --- Endpoints ---

@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Asynchronně načte PDB soubor a provede úvodní analýzu struktury.
    Detekuje alternativní lokace, chybějící atomy a kompatibilitu s forcefieldem.
    """
    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Molecule file not found on the server.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)

        # Non-blocking I/O pro čtení souboru
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Delegace CPU-bound validace do threadpoolu pro zamezení blokování event loopu
        return await run_in_threadpool(
            validation_service.validate_pdb_content,
            pdb_content,
            request.label
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Validation error for workspace {request.workspace_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Validation error: {str(e)}")


@router.post("/preview-selection", summary="Náhled geometrie před aplikací")
async def preview_selection(request: FixAltLocRequest):
    """
    Provede náhled geometrických úprav na základě uživatelských selekcí
    bez trvalého zápisu do souboru.
    """
    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Molecule file not found on the server.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Izolovaný výpočet kontinuity řetězce
        issues = await run_in_threadpool(
            validation_service.conf_manager.validate_continuity,
            pdb_content,
            request.selections
        )

        is_safe = len(issues) == 0

        return {
            "is_ok": is_safe,
            "issues": issues,
            "message": "Selection is geometrically valid" if is_safe else "Critical chain gaps detected"
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Selection preview error for workspace {request.workspace_id}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/apply-selections", summary="Aplikuje výběr konformací a ověří kontinuitu")
async def apply_selections(request: FixAltLocRequest):
    """
    Aplikuje vybrané konformace a asynchronně přepíše zdrojový PDB soubor.
    """
    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Molecule file not found on the server.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        result = await run_in_threadpool(
            validation_service.apply_alt_loc_selection,
            pdb_content,
            request.selections
        )

        # Zápis upravené struktury zpět na disk (non-blocking)
        if "pdb_content" in result:
            async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
                await f.write(result["pdb_content"])

        return result

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error applying selections for workspace {request.workspace_id}")
        raise HTTPException(status_code=500, detail=f"Error applying selections: {str(e)}")


@router.post("/prepare", summary="Kompletní příprava: Protonace, Solvatace, Ionty")
async def prepare_molecule(request: PreparationRequest):
    """
    Spouští výpočetně náročný proces přípravy struktury (PDBFixer/AMBER).
    Endpoint využívá threadpool pro paralelizaci a zachování odezvy serveru.
    """
    start_time = time.time()
    logger.info(f"Started molecule preparation for workspace: {request.workspace_id} (pH: {request.ph})")

    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            logger.error(f"Workspace {request.workspace_id} not found.")
            raise HTTPException(status_code=404, detail="Molecule file not found on the server.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        logger.info(f"Starting HydrogenationService (solvent: {request.add_solvent}, ions: {request.ionic_strength})")

        # Spuštění primární chemické transformace v dedikovaném vlákně
        prepared_pdb = await run_in_threadpool(
            hydrogenation_service.prepare_structure,
            pdb_content=pdb_content,
            ph=request.ph,
            crystal_water_mode=request.crystal_water_mode,
            add_solvent=request.add_solvent,
            box_padding_nm=request.box_padding_nm,
            ionic_strength=request.ionic_strength,
            positive_ion=request.positive_ion,
            negative_ion=request.negative_ion
        )

        if not prepared_pdb:
            logger.error("HydrogenationService returned empty output.")
            raise ValueError("The resulting PDB content is empty.")

        logger.info(f"Structure prepared successfully in {time.time() - start_time:.2f}s. Starting write process.")

        async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
            await f.write(prepared_pdb)

        # Finální validace integrity po modifikaci
        validation_results = await run_in_threadpool(
            validation_service.validate_pdb_content,
            prepared_pdb,
            label="prepared_state"
        )

        return {
            "message": "Structure successfully prepared.",
            "validation": validation_results
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(
            f"Critical failure during preparation of workspace {request.workspace_id}:\n{traceback.format_exc()}")
        raise HTTPException(
            status_code=500,
            detail={
                "error": str(e),
                "type": type(e).__name__,
                "msg": "Internal server error. Check logs for details."
            }
        )