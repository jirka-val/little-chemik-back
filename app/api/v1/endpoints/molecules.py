import logging
import io
import aiofiles
from fastapi import APIRouter, HTTPException, File, UploadFile, Body
from fastapi.concurrency import run_in_threadpool

from app.services.pdb_service import PDBService
from app.workspaces.manager import workspace_manager
from app.services.structure.hydrogenation import HydrogenationService

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()
hydrogen_service = HydrogenationService()

@router.post("/upload")
async def upload_molecule(file: UploadFile = File(...)):
    """
    Přijímá PDB soubor přes multipart form data, vytvoří workspace a vrátí jeho ID.
    """
    if not file.filename.endswith('.pdb'):
        logger.warning(f"Unsupported format attempt: {file.filename}")
        raise HTTPException(status_code=400, detail="Only .pdb files are currently supported.")

    try:
        workspace_id = await workspace_manager.create_from_upload(file)
        logger.info(f"Successfully created workspace {workspace_id} from {file.filename}")

        return {
            "workspace_id": workspace_id,
            "filename": file.filename,
            "message": "Molecule uploaded and workspace created successfully."
        }
    except Exception as e:
        logger.exception(f"Critical error saving file {file.filename}: {e}")
        raise HTTPException(status_code=500, detail="Failed to save the uploaded file.")


@router.get("/fetch-pdb/{pdb_code}")
async def fetch_pdb_by_code(pdb_code: str):
    """
    Stáhne molekulu z PDB (Protein Data Bank) asynchronně a uloží ji do workspace.
    """
    try:
        logger.info(f"Fetching PDB code: {pdb_code}")
        # Toto neblokuje event loop, protože httpx/aiohttp běží asynchronně pod kapotou
        pdb_content = await pdb_service.get_remote_pdb_content(pdb_code.lower())

        workspace_id = workspace_manager.create_from_string(pdb_content)
        logger.info(f"Successfully fetched and created workspace {workspace_id} for {pdb_code}")

        return {
            "workspace_id": workspace_id,
            "filename": f"{pdb_code}.pdb",
            "message": f"Molecule {pdb_code} successfully fetched from PDB."
        }
    except FileNotFoundError:
        logger.error(f"PDB code {pdb_code} not found.")
        raise HTTPException(status_code=404, detail=f"Molecule {pdb_code} does not exist in the RCSB database.")
    except Exception as e:
        logger.exception(f"Error fetching molecule {pdb_code} from external database: {e}")
        raise HTTPException(status_code=500, detail=f"Error fetching from PDB: {str(e)}")


@router.post("/add-hydrogens/{workspace_id}")
async def add_hydrogens(
        workspace_id: str,
        ph: float = Body(7.0, embed=True),
        optimize: bool = Body(False, embed=True)
):
    """
    Přidá atomy vodíku do molekuly a volitelně je zoptimalizuje (např. pomocí AMBER14).
    Tato CPU-náročná operace je delegována do threadpoolu, aby neblokovala FastAPI event loop.
    """
    logger.info(f"Hydrogenation request for workspace: {workspace_id} (pH: {ph}, optimize: {optimize})")

    if not workspace_manager.workspace_exists(workspace_id):
        logger.error(f"Workspace {workspace_id} not found.")
        raise HTTPException(status_code=404, detail="Workspace not found. Please upload a file first.")

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        # Asynchronní I/O pro čtení
        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        # Delegace těžkého chemického výpočtu na vedlejší vlákno
        updated_pdb_text = await run_in_threadpool(
            hydrogen_service.add_hydrogen_atoms,
            pdb_text,
            ph=ph,
            optimize=optimize
        )

        # Asynchronní I/O pro zápis výsledku
        async with aiofiles.open(file_path, "w", encoding="utf-8") as f:
            await f.write(updated_pdb_text)

        logger.info(f"Hydrogens added to {workspace_id}. Optimization: {optimize}")

        return {
            "workspace_id": workspace_id,
            "message": "Hydrogens successfully added.",
            "ph": ph,
            "optimized": optimize
        }

    except Exception as e:
        logger.exception(f"Error during hydrogenation for {workspace_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Structure modification failed: {str(e)}")