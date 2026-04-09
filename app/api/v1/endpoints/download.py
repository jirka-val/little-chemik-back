import logging
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()

@router.get("/{workspace_id}", summary="Stažení PDB souboru z workspace")
async def download_workspace(workspace_id: str):
    """
    Slouží primárně pro zobrazení struktury v Molstar na frontendu.
    Využívá asynchronní FileResponse z FastAPI, který streamuje soubor
    po částech a garantuje, že nedojde k blokování event loopu.
    """
    try:
        # Zkontrolujeme, zda soubor ještě existuje
        if not workspace_manager.workspace_exists(workspace_id):
            logger.warning(f"Download attempt for non-existent or expired workspace: {workspace_id}")
            raise HTTPException(status_code=404, detail="Workspace not found or expired.")

        file_path = workspace_manager.get_file_path(workspace_id)
        logger.info(f"Serving PDB file for workspace: {workspace_id}")

        # FileResponse zajistí bezpečné a asynchronní odeslání souboru klientovi
        return FileResponse(
            path=file_path,
            media_type="chemical/x-pdb",
            filename=f"{workspace_id}.pdb"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error serving file for workspace {workspace_id}: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error while preparing file download: {str(e)}"
        )