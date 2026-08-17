import logging
import aiofiles
from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool

from app.core.exceptions import InternalError
from app.services.pdb_service import PDBService
from app.services.forcefield_service import ForceFieldService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
pdb_service = PDBService()
ff_validator = ForceFieldService()


@router.get("/{workspace_id}", summary="Získá dostupné forcefieldy pro danou molekulu")
async def get_my_forcefields(workspace_id: str):
    """
    Asynchronně načte PDB soubor a v dedikovaném vlákně zanalyzuje
    chemické složení (typ reziduí) pro detekci kompatibilních forcefieldů.
    """
    workspace_manager.require_workspace(workspace_id)

    try:
        # Načteme PDB ze souboru ASYNCHRONNĚ
        path = workspace_manager.get_file_path(workspace_id)
        async with aiofiles.open(path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Zjistíme, co v tom je za chemii (D, R, P...)
        types = await run_in_threadpool(
            pdb_service.get_molecule_types,
            pdb_content
        )

        # "Ions" nabízíme vždy, i když struktura sama žádné krystalové ionty
        # neobsahuje - "Add Water Box & Ions" (Hydrogens tab) umí přidat
        # neutralizační/solné ionty (Na+/Cl-...) bez ohledu na to, jestli
        # tam nějaké původně byly, a builder pro ně potřebuje LJ parametry
        # vybrané předem, ne až jako reakci na pád v /prepare.
        if "I" not in types:
            types.append("I")

        # Spárujeme ty správné FF
        ffs = await run_in_threadpool(
            ff_validator.get_matching_forcefields,
            types
        )

        return {
            "detected_types": types,
            "forcefields": ffs
        }

    except Exception as e:
        logger.exception(f"Error fetching forcefields for workspace {workspace_id}")
        raise InternalError(f"Internal server error while fetching forcefields: {str(e)}")