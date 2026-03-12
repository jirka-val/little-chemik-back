import logging
from fastapi import APIRouter, HTTPException
from app.workspaces.manager import workspace_manager
from app.services.analysis_service import build_sequence_tokens

logger = logging.getLogger("api")
router = APIRouter()


@router.get("/sequence/{workspace_id}")
async def analyze_sequence(workspace_id: str, chain: str | None = None, fill_gaps: bool = True):
    """
    Vezme workspace_id, přečte dočasný soubor na pozadí a vrátí sekvenci
    včetně informací o chybějících atomech (tvůj converting_dictionary).
    """
    logger.info(f"Požadavek na analýzu sekvence pro workspace: {workspace_id}")

    # 1. Zkontrolujeme, zda soubor ještě žije (nevypršel čas)
    if not workspace_manager.workspace_exists(workspace_id):
        logger.error(f"Workspace {workspace_id} nebyl nalezen.")
        raise HTTPException(status_code=404, detail="Workspace nenalezen. Nahráli jste soubor?")

    try:
        # 2. Zjistíme si cestu k souboru na disku
        file_path = workspace_manager.get_file_path(workspace_id)

        # 3. Přečteme soubor bezpečně do paměti (pro tvou analýzu)
        with open(file_path, "r", encoding="utf-8") as f:
            pdb_text = f.read()

        # TADY JE TA ÚPRAVA - dáme fill_gaps natvrdo na True
        sequence_data = build_sequence_tokens(
            pdb_text=pdb_text,
            chain=chain,
            fill_gaps=True  # <-- Vynuceno! Backend teď MUSÍ hledat chybějící atomy
        )

        logger.info(f"Analýza pro {workspace_id} byla úspěšně dokončena.")
        return {
            "workspace_id": workspace_id,
            "sequence": sequence_data
        }

    except Exception as e:
        logger.exception(f"Chyba při zpracování sekvence pro workspace {workspace_id}: {e}")
        raise HTTPException(status_code=500, detail="Interní chyba při zpracování molekuly.")