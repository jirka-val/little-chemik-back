from fastapi import APIRouter, HTTPException
from app.services.pdb_service import PDBService
import logging

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()

@router.get("/process/{pdb_code}", summary="Zpracuje PDB kód molekuly")
async def process_molecule(pdb_code: str):
    """
    Stáhne a analyzuje molekulu na základě jejího RCSB kódu.
    """
    try:
        content = await pdb_service.fetch_pdb_content(pdb_code.lower())
        return {"pdb_code": pdb_code, "data": content}
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception(f"Chyba při zpracování {pdb_code}")
        raise HTTPException(status_code=500, detail="Interní chyba při zpracování molekuly.")