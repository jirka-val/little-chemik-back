import logging
import aiofiles
from typing import Dict

from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.core.exceptions import InternalError
from app.services.analysis import list_ion_options
from app.services.ff_catalog_service import catalog_service
from app.services.validation.service import ValidationService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
validation_service = ValidationService()


class ValidationRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    label: str = Field("molecule_from_front", description="Molecule state identifier")


class FixAltLocRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    selections: Dict[str, str] = Field(
        ...,
        description="Variant selection map, e.g., {'A-42': 'B'}",
        json_schema_extra={"example": {"A-42": "B", "A-15": "A"}},
    )


# --- Endpoints ---

@router.get("/ions", summary="Vrátí vybíratelné ionty pro salt/neutralizaci po mol_type skupině")
async def get_ion_options():
    """
    Zdroj pravdy pro pos./neg. ion dropdowny na frontendu - converting_dictionary.json
    (chemická identita, ne z hardcoded seznamu - viz poznámka u
    required_ff_groups o třech nezávislých kopiích, co tenhle seznam dřív
    měl a jak to skončilo) PROTNUTÉ s FFCatalogService.get_buildable_ion_resnames()
    (jestli pro ten ion v aktuálním katalogu vůbec existuje force field).
    Bez tohohle průniku by šlo z dropdownu vybrat ion (typicky exotický
    lanthanid/aktinid z Im+, nebo nábojem-sufixovanou variantu tam, kde
    katalog zná jen starší holý symbol), který by vždy skončil neřešitelnou
    chybou - buď 409 "missing force field" (nikde není co vybrat), nebo
    KeyError hluboko v builderu. "I1"/"I1+" jsou monovalentní, "Im"/"Im+"
    dvojmocné+ (Im obsahuje Mg2+).
    """
    known = list_ion_options()
    buildable = catalog_service.get_buildable_ion_resnames()
    return {
        mol_type: sorted(set(resnames) & buildable.get(mol_type, set()))
        for mol_type, resnames in known.items()
    }


@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Asynchronně načte PDB soubor a provede úvodní analýzu struktury.
    Detekuje alternativní lokace, chybějící atomy a kompatibilitu s forcefieldem.
    """
    workspace_manager.require_workspace(request.workspace_id)
    logger.info(f"Validating structure for workspace {request.workspace_id} (label: {request.label})...")

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)

        # Non-blocking I/O pro čtení souboru
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Delegace CPU-bound validace do threadpoolu pro zamezení blokování event loopu
        result = await run_in_threadpool(
            validation_service.validate_pdb_content,
            pdb_content,
            request.label
        )

        summary = result.get("summary", {})
        analysis = result.get("analysis", {})
        logger.info(
            f"Validation finished for workspace {request.workspace_id}: "
            f"ready_for_hpc={summary.get('is_ready_for_hpc')}, "
            f"errors={len(analysis.get('errors', []))}, warnings={len(analysis.get('warnings', []))}"
        )
        return result

    except Exception as e:
        logger.exception(f"Validation error for workspace {request.workspace_id}: {str(e)}")
        raise InternalError(f"Validation error: {str(e)}")


@router.post("/preview-selection", summary="Náhled geometrie před aplikací")
async def preview_selection(request: FixAltLocRequest):
    """
    Provede náhled geometrických úprav na základě uživatelských selekcí
    bez trvalého zápisu do souboru.
    """
    workspace_manager.require_workspace(request.workspace_id)

    try:
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

    except Exception as e:
        logger.exception(f"Selection preview error for workspace {request.workspace_id}")
        raise InternalError(str(e))


@router.post("/apply-selections", summary="Aplikuje výběr konformací a ověří kontinuitu")
async def apply_selections(request: FixAltLocRequest):
    """
    Aplikuje vybrané konformace a asynchronně přepíše zdrojový PDB soubor.
    """
    workspace_manager.require_workspace(request.workspace_id)
    logger.info(
        f"Applying {len(request.selections)} AltLoc selection(s) for workspace "
        f"{request.workspace_id} (overwrites structure.pdb)..."
    )

    try:
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

        logger.info(f"AltLoc selections applied and written for workspace {request.workspace_id}.")
        return result

    except Exception as e:
        logger.exception(f"Error applying selections for workspace {request.workspace_id}")
        raise InternalError(f"Error applying selections: {str(e)}")
