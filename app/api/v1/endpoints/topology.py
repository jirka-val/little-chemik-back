import logging
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse  # Přidán import pro stažení
from pydantic import BaseModel
from typing import Dict, Any
import os

from app.services.topology_service import TopologyService
from app.workspaces.manager import workspace_manager

router = APIRouter()
topology_service = TopologyService()
logger = logging.getLogger("api")

class TopologyRequest(BaseModel):
    pdb_filename: str = "structure.pdb"
    ff_selections: Dict[str, Any]

@router.post("/{workspace_id}/generate")
async def generate_topology(workspace_id: str, request: TopologyRequest):
    """
    Endpoint vygeneruje topologii a rovnou ji pošle uživateli ke stažení.
    """
    if not workspace_manager.workspace_exists(workspace_id):
        raise HTTPException(status_code=404, detail="Workspace not found")

    try:
        logger.info(f"Starting topology generation for workspace: {workspace_id}")

        # 1. Necháme servis vygenerovat soubor (vrátí název, např. 'structure.prmtop')
        output_filename = topology_service.generate_topology(
            workspace_id=workspace_id,
            pdb_filename=request.pdb_filename,
            ff_selections=request.ff_selections
        )

        # 2. Získáme absolutní cestu k novému souboru
        file_path = workspace_manager.get_workspace_dir(workspace_id) / output_filename

        if not file_path.exists():
            raise HTTPException(status_code=500, detail="Generated file not found on server.")

        # 3. Vrátíme FileResponse, která vynutí stažení
        return FileResponse(
            path=file_path,
            filename=output_filename,
            media_type='application/octet-stream'
        )

    except Exception as e:
        logger.error(f"Error in topology endpoint for {workspace_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))