from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from app.services.pdb_service import PDBService
from app.services.analysis_service import build_sequence_tokens
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


@router.get("/fetch/{pdb_code}/sequence", summary="Stáhne PDB a vrátí i sequence tokens")
async def fetch_molecule_sequence(
    pdb_code: str,
    chain: str | None = None,
    fill_gaps: bool = True
):
    try:
        pdb_text = await pdb_service.get_remote_pdb_content(pdb_code.lower())
        return {
            "pdb_code": pdb_code.lower(),
            "sequence": build_sequence_tokens(
                pdb_text=pdb_text,
                chain=chain,
                fill_gaps=fill_gaps
            )
        }
    except FileNotFoundError as e:
        logger.error(f"Molekula nenalezena v externí databázi: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception:
        logger.exception(f"Chyba při zpracování sekvence pro {pdb_code}")
        raise HTTPException(status_code=500, detail="Interní chyba při zpracování molekuly.")
