import logging
import aiofiles
from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool

from app.workspaces.manager import workspace_manager
from app.services.analysis_service import build_sequence_tokens

import httpx

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

@router.get("/analyze-pdb/{pdb_code}", summary="Vzdálená analýza PDB z RCSB bez workspace")
async def analyze_remote_pdb(pdb_code: str, chain: str | None = None, fill_gaps: bool = True):
    """
    Endpoint pro externí skripty (např. pro šéfa).
    Stáhne PDB soubor přímo z RCSB podle kódu, asynchronně ho přečte
    a v dedikovaném vlákně vrátí jeho surový text a analýzu chybějících atomů.
    """
    logger.info(f"Remote PDB analysis requested for code: {pdb_code} (chain: {chain}, fill_gaps: {fill_gaps})")

    # 1. Stažení PDB souboru z RCSB přes asynchronního klienta
    # Používáme oficiální URL adresu pro download, kterou máte v projektu
    rcsb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"

    async with httpx.AsyncClient() as client:
        try:
            response = await client.get(rcsb_url)

            if response.status_code == 404:
                logger.error(f"PDB kód {pdb_code} nebyl v RCSB databázi nalezen.")
                raise HTTPException(status_code=404, detail=f"Molekula s kódem {pdb_code} neexistuje v RCSB databázi.")

            elif response.status_code != 200:
                logger.error(f"RCSB vrátilo neočekávaný stav: {response.status_code}")
                raise HTTPException(status_code=502, detail="Chyba při komunikaci s RCSB databází.")

            pdb_text = response.text

        except httpx.RequestError as e:
            logger.exception(f"Síťová chyba při stahování molekuly {pdb_code}: {str(e)}")
            raise HTTPException(status_code=503, detail="RCSB databáze je momentálně nedostupná.")

    # 2. Spuštění CPU-náročné analýzy tokenů ve vedlejším vlákně (stejně jako u workspace)
    try:
        sequence_data = await run_in_threadpool(
            build_sequence_tokens,
            pdb_text=pdb_text,
            chain=chain,
            fill_gaps=fill_gaps
        )

        logger.info(f"Remote analysis for PDB {pdb_code} successfully completed.")

        # Vracíme přesně ty klíče, které šéfův skript očekává v response.json()
        return {
            "pdb_text": pdb_text,
            "missing_atoms": sequence_data
        }

    except Exception as e:
        logger.exception(f"Chyba při zpracování sekvence pro vzdálené PDB {pdb_code}: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Interní chyba serveru při analýze molekulární sekvence."
        )