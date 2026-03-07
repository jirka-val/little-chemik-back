from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from app.workspaces.manager import workspace_manager

router = APIRouter()


@router.get("/{workspace_id}")
async def download_workspace(workspace_id: str):
    """
    Slouží primárně pro Molstar na frontendu.
    Molstar si přes tuto URL bleskově stáhne 3D PDB data.
    """
    # 1. Zkontrolujeme, zda soubor ještě existuje (jestli už nevypršel jeho čas)
    if not workspace_manager.workspace_exists(workspace_id):
        raise HTTPException(status_code=404, detail="Workspace not found or expired.")

    file_path = workspace_manager.get_file_path(workspace_id)

    # 2. FileResponse zajistí přímé a bleskové odeslání souboru klientovi
    return FileResponse(
        path=file_path,
        media_type="chemical/x-pdb",
        filename=f"{workspace_id}.pdb"
    )