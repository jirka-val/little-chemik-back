import logging
from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel
from typing import Dict, Any

from app.core.exceptions import AppBaseException, InternalError
from app.services.topology_service import TopologyService
from app.workspaces.manager import workspace_manager
from app.core.job_queue import heavy_job
from app.core.process_pool import run_heavy
from app.services.ff_selection import resolve_ff_selections

router = APIRouter()
topology_service = TopologyService()
logger = logging.getLogger(__name__)


class TopologyRequest(BaseModel):
    pdb_filename: str = "structure.pdb"
    ff_selections: Dict[str, Any]
    # Hydrogen Mass Repartitioning - volí se v Simulation panelu (4 fs timestep).
    hmr: bool = False


def topology_job(**kwargs):
    """Generování topologie v procesu poolu (app/core/process_pool.py)."""
    return topology_service.generate_topology(**kwargs)


@router.post("/{workspace_id}/generate")
@heavy_job("topology")
async def generate_topology(workspace_id: str, request: TopologyRequest):
    workspace_manager.require_workspace(workspace_id)
    # Obsah FF souborů jen z katalogu, ne od klienta (app/services/ff_selection.py).
    ff_selections = await run_in_threadpool(resolve_ff_selections, request.ff_selections)
    logger.info(f"Generating topology for workspace {workspace_id} (ff_selections: {list(request.ff_selections.keys())}, hmr={request.hmr})...")

    try:
        # 1. Vygenerujeme soubory na disk (vrátí dict s názvy). V threadpoolu:
        # generování trvá sekundy (u 300k atomů ~6 s) a přímo v event loopu by
        # po tu dobu server neodpověděl nikomu jinému (viz tests/load/).
        result_dict = await run_heavy(
            topology_job,
            workspace_id=workspace_id,
            pdb_filename=request.pdb_filename,
            ff_selections=ff_selections,
            hmr=request.hmr,
        )

        logger.info(f"Topology generated successfully for workspace {workspace_id}.")

        # 2. VRÁTÍME ČISTÝ JSON (Žádný StreamingResponse, žádný ZIP!)
        return {
            "status": "success",
            "hmr": result_dict["hmr"],
            "total_charge": result_dict["total_charge"],
            "files": {
                "topology": result_dict["topology_file"],
                "pdb": result_dict["coordinates_file"]
            }
        }
    except AppBaseException:
        raise
    except Exception as e:
        logger.exception(f"Topology generation failed for workspace {workspace_id}: {e}")
        raise InternalError(str(e))