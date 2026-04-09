import logging
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any

from app.services.topology_service import TopologyService
from app.workspaces.manager import workspace_manager

router = APIRouter()
topology_service = TopologyService()
logger = logging.getLogger("api")


# Definice schématu požadavku
class TopologyRequest(BaseModel):
    pdb_filename: str = "structure.pdb"
    # ff_selections obsahuje mapování např. {"R": ff_data_objekt}
    ff_selections: Dict[str, Any]


@router.post("/{workspace_id}/generate")
async def generate_topology(workspace_id: str, request: TopologyRequest):
    """
    Endpoint pro generování AMBER topologie (.prmtop).
    Vstupy: ID workspace a výběr silových polí z API.
    """
    if not workspace_manager.workspace_exists(workspace_id):
        logger.warning(f"Topology request failed: Workspace {workspace_id} not found.")
        raise HTTPException(status_code=404, detail="Workspace not found")

    try:
        logger.info(f"Starting topology generation for workspace: {workspace_id}")

        # Spuštění orchestrace v TopologyService
        output_filename = topology_service.generate_topology(
            workspace_id=workspace_id,
            pdb_filename=request.pdb_filename,
            ff_selections=request.ff_selections
        )

        return {
            "status": "success",
            "workspace_id": workspace_id,
            "prmtop_file": output_filename,
            "message": "Topology generated and saved to workspace."
        }

    except Exception as e:
        logger.error(f"Error in topology endpoint for {workspace_id}: {str(e)}")
        # Pokud šéfova knihovna vyhodí chybu (např. chybějící parametry), pošleme ji uživateli
        raise HTTPException(status_code=500, detail=str(e))