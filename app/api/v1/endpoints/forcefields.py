import logging
import aiofiles
from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool

from app.services.pdb_service import PDBService
from app.services.validation.forcefield import ForceFieldValidator
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
pdb_service = PDBService()
ff_validator = ForceFieldValidator()


@router.get("/{workspace_id}", summary="Získá dostupné forcefieldy pro danou molekulu")
async def get_my_forcefields(workspace_id: str):
    """
    Asynchronně načte PDB soubor a v dedikovaném vlákně zanalyzuje
    chemické složení (typ reziduí) pro detekci kompatibilních forcefieldů.
    """
    try:
        # Existuje ten workspace?
        if not workspace_manager.workspace_exists(workspace_id):
            logger.warning(f"Workspace {workspace_id} not found.")
            raise HTTPException(status_code=404, detail="Molecule file not found on the server.")

        # Načteme PDB ze souboru ASYNCHRONNĚ
        path = workspace_manager.get_file_path(workspace_id)
        async with aiofiles.open(path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Zjistíme, co v tom je za chemii (D, R, P...)
        types = await run_in_threadpool(
            pdb_service.get_molecule_types,
            pdb_content
        )

        # Spárujeme ty správné FF
        ffs = await run_in_threadpool(
            ff_validator.get_matching_forcefields,
            types
        )

        return {
            "detected_types": types,
            "forcefields": ffs
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error fetching forcefields for workspace {workspace_id}")
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error while fetching forcefields: {str(e)}"
        )