from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from app.services.pdb_service import PDBService
import logging

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()

@router.get("/process/{pdb_code}", summary="Zpracuje PDB kód molekuly", response_class=PlainTextResponse)
async def process_molecule(pdb_code: str):

    try:
        content = await pdb_service.fetch_pdb_content(pdb_code.lower())
        return content
    except FileNotFoundError as e:
        logger.error(f"Soubor nenalezen: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception(f"Chyba při zpracování {pdb_code}")
        raise HTTPException(status_code=500, detail="Interní chyba při zpracování molekuly.")