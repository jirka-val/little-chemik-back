import logging
import aiofiles
import httpx
from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool

from pydantic import BaseModel
from typing import Dict, Any

from app.core.exceptions import ExternalServiceError, InternalError, RemoteMoleculeNotFoundError
from app.workspaces.manager import workspace_manager
# ZMĚNA 1: Importujeme novou funkci process_structure místo původní clean_pdb_altlocs
from app.services.analysis_service import build_sequence_tokens, analyze_pdb_altlocs, process_structure

logger = logging.getLogger("api")
router = APIRouter()


# ZMĚNA 2: Rozšířené Pydantic schéma pro komplexní přípravu struktury
class StructurePrepRequest(BaseModel):
    model: int = 1
    apply_symmetry: bool = False
    selection: Dict[str, str] = {}


@router.get("/sequence/{workspace_id}", summary="Analýza sekvence molekuly")
async def analyze_sequence(workspace_id: str, chain: str | None = None, fill_gaps: bool = True, filename: str = "structure.pdb"):
    """
    Vezme workspace_id, asynchronně přečte dočasný soubor na pozadí a
    v dedikovaném vlákně vrátí sekvenci včetně informací o chybějících atomech
    (pomocí converting_dictionary).

    `filename` (stejný vzor jako /api/download) umožňuje analyzovat i jiný
    soubor než finální structure.pdb - typicky "structure_preview.pdb" během
    otevřené side-chain GUI relace (viz sidechains.py), kdy builder už doplnil
    vodíky/atomy pro většinu reziduí, ale ještě nebyl commit.
    """
    logger.info(f"Sequence analysis requested for workspace: {workspace_id} (chain: {chain}, fill_gaps: {fill_gaps}, filename: {filename})")

    workspace_manager.require_workspace(workspace_id)

    try:
        file_path = workspace_manager.get_file_path(workspace_id, filename)

        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

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

    except Exception as e:
        logger.exception(f"Error processing sequence for workspace {workspace_id}: {str(e)}")
        raise InternalError("Internal server error while processing the molecule sequence.")


@router.get("/analyze-pdb/{pdb_code}", summary="Vzdálená analýza PDB z RCSB bez workspace")
async def analyze_remote_pdb(pdb_code: str, chain: str | None = None, fill_gaps: bool = True):
    """
    Endpoint pro externí skripty.
    Stáhne PDB soubor přímo z RCSB podle kódu, asynchronně ho přečte,
    AUTOMATICKY VYČISTÍ A PŘIPRAVÍ STRUKTURU
    a v dedikovaném vlákně vrátí jeho čistý text a analýzu chybějících atomů.
    """
    logger.info(f"Remote PDB analysis requested for code: {pdb_code} (chain: {chain}, fill_gaps: {fill_gaps})")

    rcsb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"

    async with httpx.AsyncClient() as client:
        try:
            response = await client.get(rcsb_url)

            if response.status_code == 404:
                logger.error(f"PDB kód {pdb_code} nebyl v RCSB databázi nalezen.")
                raise RemoteMoleculeNotFoundError(pdb_code)

            elif response.status_code != 200:
                logger.error(f"RCSB vrátilo neočekávaný stav: {response.status_code}")
                raise ExternalServiceError("Chyba při komunikaci s RCSB databází.", status_code=502)

            pdb_text = response.text

        except httpx.RequestError as e:
            logger.exception(f"Síťová chyba při stahování molekuly {pdb_code}: {str(e)}")
            raise ExternalServiceError("RCSB databáze je momentálně nedostupná.", status_code=503)

        # Auto-příprava struktury pro vzdálené PDB
        try:
            result = await run_in_threadpool(analyze_pdb_altlocs, pdb_text)
            altlocs_data = result.get("residues", []) if isinstance(result, dict) else []

            auto_selection = {}
            if altlocs_data and isinstance(altlocs_data, list):
                for alt_item in altlocs_data:
                    chain_id = alt_item.get("chain", "?")
                    resseq = str(alt_item.get("resseq", ""))
                    resname = alt_item.get("resname", "")
                    key = f"{chain_id}_{resseq}_{resname}"

                    alt_locs_dict = alt_item.get("altLocs", {})
                    variants = list(alt_locs_dict.keys())
                    chosen_variant = variants[0] if variants else None

                    if key and chosen_variant:
                        auto_selection[key] = chosen_variant

            # ZMĚNA 3: Voláme novou process_structure (vzdáleně standardně bereme Model 1 a aplikujeme symetrii)
            has_symmetry = result.get("hasSymmetry", False) if isinstance(result, dict) else False
            pdb_text = await run_in_threadpool(
                process_structure,
                pdb_text=pdb_text,
                target_model=1,
                apply_symmetry=has_symmetry,
                selection=auto_selection
            )
        except Exception as e:
            logger.exception(f"Chyba při automatické přípravě struktury pro {pdb_code}: {e}")

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
        raise InternalError("Interní chyba serveru při analýze molekulární sekvence.")


@router.get("/altlocs/{workspace_id}", summary="Analýza struktury (Modely, Symetrie, AltLocs)")
async def analyze_altlocs(workspace_id: str):
    """
    Zkontroluje PDB soubor ve workspace a vrátí JSON s informacemi o:
    - Přítomnosti více modelů (NMR ensemble)
    - Přítomnosti REMARK 350 (Biological Assembly matic)
    - Alternativních pozicích (AltLocs), jejich obsazenosti a B-faktorech.
    """
    logger.info(f"Structure prep analysis requested for workspace: {workspace_id}")

    workspace_manager.require_workspace(workspace_id)

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        altloc_data = await run_in_threadpool(
            analyze_pdb_altlocs,
            pdb_text=pdb_text
        )

        logger.info(f"Structure prep analysis for workspace {workspace_id} successfully completed.")

        return altloc_data

    except Exception as e:
        logger.exception(f"Error processing structure prep analysis for workspace {workspace_id}: {str(e)}")
        raise InternalError("Internal server error while analyzing structure.")


@router.post("/clean-altlocs/{workspace_id}", summary="Aplikuje výběr Modelu, Symetrie a AltLocs")
async def apply_clean_altlocs(workspace_id: str, payload: StructurePrepRequest):
    """
    Přijme od frontenyl komplexní konfigurační požadavek:
    - Vybraný číslo modelu (`model`)
    - Příznak pro vybudování biologické symetrie (`apply_symmetry`)
    - Zvolené alternativní pozice reziduí (`selection`)

    Aplikuje úpravy a přepíše PDB soubor ve workspace.
    """
    logger.info(
        f"Preparing structure for workspace: {workspace_id} (Model: {payload.model}, Symmetry: {payload.apply_symmetry})")

    workspace_manager.require_workspace(workspace_id)

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        # Načtení stávajícího PDB
        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        # ZMĚNA 4: Vyčištění a sestavení přes novou funkci process_structure
        cleaned_pdb_text = await run_in_threadpool(
            process_structure,
            pdb_text=pdb_text,
            target_model=payload.model,
            apply_symmetry=payload.apply_symmetry,
            selection=payload.selection
        )

        # Přepsání PDB souboru vyčištěnou strukturou
        async with aiofiles.open(file_path, "w", encoding="utf-8") as f:
            await f.write(cleaned_pdb_text)

        logger.info(f"Structure prepared successfully for workspace {workspace_id}.")
        return {"status": "success", "message": "Structure prepared and cleaned successfully."}

    except Exception as e:
        logger.exception(f"Error preparing structure for workspace {workspace_id}: {str(e)}")
        raise InternalError("Error preparing structure.")
