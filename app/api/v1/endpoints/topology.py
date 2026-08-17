import logging
import zipfile
import io
from fastapi import APIRouter
from fastapi.responses import StreamingResponse  # Změněno z FileResponse
from pydantic import BaseModel
from typing import Dict, Any
import os

from app.core.exceptions import InternalError
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
    workspace_manager.require_workspace(workspace_id)

    try:
        # 1. Vygenerujeme soubory na disk (vrátí dict s názvy)
        result_dict = topology_service.generate_topology(
            workspace_id=workspace_id,
            pdb_filename=request.pdb_filename,
            ff_selections=request.ff_selections
        )

        # 2. VRÁTÍME ČISTÝ JSON (Žádný StreamingResponse, žádný ZIP!)
        return {
            "status": "success",
            "files": {
                "topology": result_dict["topology_file"],
                "pdb": result_dict["coordinates_file"]
            }
        }
    except Exception as e:
        logger.error(f"Error: {e}")
        raise InternalError(str(e))