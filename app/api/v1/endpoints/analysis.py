import logging
import aiofiles
import httpx
from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool

from pydantic import BaseModel
from typing import Dict


from app.workspaces.manager import workspace_manager
from app.services.analysis_service import build_sequence_tokens, analyze_pdb_altlocs, clean_pdb_altlocs

logger = logging.getLogger("api")
router = APIRouter()


# --- PŘIDANÝ MODEL PRO PŘIJETÍ DAT Z FRONTENDU ---
class AltLocSelectionRequest(BaseModel):
    selection: Dict[str, str]
# -------------------------------------------------


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
    Endpoint pro externí skripty.
    Stáhne PDB soubor přímo z RCSB podle kódu, asynchronně ho přečte,
    AUTOMATICKY VYČISTÍ ALTERNATIVNÍ POZICE (AltLocs)
    a v dedikovaném vlákně vrátí jeho čistý text a analýzu chybějících atomů.
    """
    logger.info(f"Remote PDB analysis requested for code: {pdb_code} (chain: {chain}, fill_gaps: {fill_gaps})")

    # 1. Stažení PDB souboru z RCSB přes asynchronního klienta
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

        # ==============================================================================
        try:
            # Najdeme všechny alternativní pozice v souboru (vrací slovník!)
            result = await run_in_threadpool(analyze_pdb_altlocs, pdb_text)

            # OPRAVA: Vytáhneme seznam reziduí ze slovníku
            altlocs_data = result.get("residues", []) if isinstance(result, dict) else []

            if altlocs_data and isinstance(altlocs_data, list):
                auto_selection = {}
                for alt_item in altlocs_data:
                    # Klíč nyní poskládáme ručně, protože struktura z analyze_pdb_altlocs
                    # vrací přímo chain, resseq, resname
                    chain = alt_item.get("chain", "?")
                    resseq = str(alt_item.get("resseq", ""))
                    resname = alt_item.get("resname", "")
                    key = f"{chain}_{resseq}_{resname}"

                    # Najdeme recommended variantu (altLocs je slovník, např. {"A": {...}, "B": {...}})
                    alt_locs_dict = alt_item.get("altLocs", {})

                    # Zkusíme najít první klíč, který není mezera (v našem případě 'A' nebo 'B')
                    # Pokud bys chtěl chytřejší logiku (třeba podle occupancy), musel bys tu iterovat
                    variants = list(alt_locs_dict.keys())
                    chosen_variant = variants[0] if variants else None

                    if key and chosen_variant:
                        auto_selection[key] = chosen_variant

                if auto_selection:
                    logger.info(f"Auto-resolving AltLocs pro {pdb_code}: {auto_selection}")
                    pdb_text = await run_in_threadpool(clean_pdb_altlocs, pdb_text, auto_selection)
            else:
                logger.info(f"Žádné AltLocs k vyřešení pro {pdb_code}.")
        except Exception as e:
            logger.exception(f"Chyba při automatickém čištění: {e}")
        # ==============================================================================

    # 2. Spuštění CPU-náročné analýzy tokenů ve vedlejším vlákně (nyní nad ČISTÝM pdb_text)
    try:
        sequence_data = await run_in_threadpool(
            build_sequence_tokens,
            pdb_text=pdb_text,
            chain=chain,
            fill_gaps=fill_gaps
        )

        logger.info(f"Remote analysis for PDB {pdb_code} successfully completed.")

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


@router.get("/altlocs/{workspace_id}", summary="Analýza alternativních pozic (AltLocs)")
async def analyze_altlocs(workspace_id: str):
    """
    Zkontroluje PDB soubor ve workspace a vrátí JSON s nalezenými
    alternativními pozicemi, jejich obsazeností (occupancy) a B-faktory.
    """
    logger.info(f"AltLoc analysis requested for workspace: {workspace_id}")

    if not workspace_manager.workspace_exists(workspace_id):
        logger.error(f"Workspace {workspace_id} not found.")
        raise HTTPException(status_code=404, detail="Workspace not found. Have you uploaded a file?")

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        altloc_data = await run_in_threadpool(
            analyze_pdb_altlocs,
            pdb_text=pdb_text
        )

        logger.info(f"AltLoc analysis for workspace {workspace_id} successfully completed.")

        return altloc_data

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error processing AltLocs for workspace {workspace_id}: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Internal server error while analyzing AltLocs."
        )


@router.post("/clean-altlocs/{workspace_id}", summary="Aplikuje výběr AltLocs a vyčistí PDB")
async def apply_clean_altlocs(workspace_id: str, payload: AltLocSelectionRequest):
    """
    Přijme od frontendu zvolené konformace (např. {"A_45_TYR": "A"}),
    smaže z PDB souboru nevybrané atomy a soubor přepíše.
    """
    logger.info(f"Cleaning AltLocs for workspace: {workspace_id}")

    if not workspace_manager.workspace_exists(workspace_id):
        raise HTTPException(status_code=404, detail="Workspace not found.")

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        # Načtení starého PDB
        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        # Vyčištění přes naši novou funkci (Fáze 4)
        cleaned_pdb_text = await run_in_threadpool(
            clean_pdb_altlocs,
            pdb_text=pdb_text,
            user_selection=payload.selection
        )

        # Přepsání starého PDB tím čistým
        async with aiofiles.open(file_path, "w", encoding="utf-8") as f:
            await f.write(cleaned_pdb_text)

        logger.info(f"AltLocs cleaned successfully for workspace {workspace_id}.")
        return {"status": "success", "message": "Structure cleaned successfully."}

    except Exception as e:
        logger.exception(f"Error cleaning AltLocs for workspace {workspace_id}: {str(e)}")
        raise HTTPException(status_code=500, detail="Error cleaning structure.")