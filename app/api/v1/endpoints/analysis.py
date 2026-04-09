import logging
import aiofiles
from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool

from app.workspaces.manager import workspace_manager
from app.services.analysis_service import build_sequence_tokens

logger = logging.getLogger("api")
router = APIRouter()


@router.get("/sequence/{workspace_id}", summary="Analýza sekvence molekuly")
async def analyze_sequence(workspace_id: str, chain: str | None = None, fill_gaps: bool = True):
    """
    Vezme workspace_id, asynchronně přečte dočasný soubor na pozadí a
    v dedikovaném vlákně vrátí sekvenci včetně informací o chybějících atomech
    (pomocí converting_dictionary).
    """
    logger.info(f"Sequence analysis requested for workspace: {workspace_id} (chain: {chain}, fill_gaps: {fill_gaps})")

    # Zkontrolujeme, zda soubor ještě žije
    if not workspace_manager.workspace_exists(workspace_id):
        logger.error(f"Workspace {workspace_id} not found.")
        raise HTTPException(status_code=404, detail="Workspace not found. Have you uploaded a file?")

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        # Přečteme soubor do paměti
        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        # Spustíme CPU-náročnou analýzu tokenů ve vedlejším vlákně
        sequence_data = await run_in_threadpool(
            build_sequence_tokens,
            pdb_text=pdb_text,
            chain=chain,
            fill_gaps=True
        )

        logger.info(f"Sequence analysis for workspace {workspace_id} successfully completed.")

        return {
            "workspace_id": workspace_id,
            "sequence": sequence_data
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error processing sequence for workspace {workspace_id}: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Internal server error while processing the molecule sequence."
        )