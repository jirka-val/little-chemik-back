from fastapi import APIRouter, HTTPException, Body
from fastapi.responses import PlainTextResponse
from app.services.pdb_service import PDBService
from app.services.analysis_service import build_sequence_tokens
# Import nové služby
from app.services.structure.hydrogenation import HydrogenationService
import logging

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()
# Instance hydrogenační služby
hydrogen_service = HydrogenationService()


@router.get("/fetch/{pdb_code}", summary="Stáhne PDB data a pošle je přímo na frontend",
            response_class=PlainTextResponse)
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


# NOVÝ ENDPOINT
@router.post("/add-hydrogens", summary="Doplní chybějící vodíky do poskytnutého PDB obsahu")
async def add_hydrogens(
        pdb_content: str = Body(..., embed=True),
        ph: float = Body(7.0, embed=True)
):
    """
    Přijme PDB jako string, použije PDBFixer k doplnění vodíků a vrátí upravené PDB.
    """
    try:
        if not pdb_content:
            raise HTTPException(status_code=400, detail="PDB content je prázdný.")

        logger.info(f"Probíhá doplňování vodíků při pH {ph}")
        updated_pdb = hydrogen_service.add_hydrogen_atoms(pdb_content, ph=ph)

        return {
            "status": "success",
            "pdb_content": updated_pdb
        }
    except Exception as e:
        logger.exception("Chyba při doplňování vodíků")
        raise HTTPException(status_code=500, detail=f"Chyba při zpracování struktury: {str(e)}")