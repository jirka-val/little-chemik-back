from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from app.services.pdb_service import PDBService
import logging

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()

@router.get("/fetch/{pdb_code}", summary="Stáhne PDB data a pošle je přímo na frontend", response_class=PlainTextResponse)
async def fetch_molecule(pdb_code: str):

    try:
        content = await pdb_service.get_remote_pdb_content(pdb_code.lower())
        return content
    except FileNotFoundError as e:
        logger.error(f"Molekula nenalezena v externí databázi: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception(f"Chyba při zprostředkování PDB dat pro {pdb_code}")
        raise HTTPException(status_code=500, detail="Interní chyba při stahování molekuly.")