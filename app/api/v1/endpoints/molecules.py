import logging
from fastapi import APIRouter, File, UploadFile, Request
from fastapi.concurrency import run_in_threadpool

from app.core.exceptions import BadRequestError, InternalError, NotFoundError, RemoteMoleculeNotFoundError
from app.services.pdb_service import PDBService, remove_residue_from_pdb
from app.workspaces.manager import workspace_manager
from app.services.incidents import history as incident_history

logger = logging.getLogger(__name__)
router = APIRouter()
pdb_service = PDBService()

@router.post("/upload")
async def upload_molecule(file: UploadFile = File(...)):
    """
    Přijímá PDB soubor přes multipart form data, vytvoří workspace a vrátí jeho ID.
    """
    if not file.filename.endswith('.pdb'):
        logger.warning(f"Unsupported format attempt: {file.filename}")
        raise BadRequestError("Only .pdb files are currently supported.")

    try:
        workspace_id = await workspace_manager.create_from_upload(file)
        logger.info(f"Successfully created workspace {workspace_id} from {file.filename}")
        await run_in_threadpool(incident_history.record_workspace_created, workspace_id, f"upload:{file.filename}")

        return {
            "workspace_id": workspace_id,
            "filename": file.filename,
            "message": "Molecule uploaded and workspace created successfully."
        }
    except Exception as e:
        logger.exception(f"Critical error saving file {file.filename}: {e}")
        raise InternalError("Failed to save the uploaded file.")


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
        await run_in_threadpool(incident_history.record_workspace_created, workspace_id, f"rcsb:{pdb_code}")

        return {
            "workspace_id": workspace_id,
            "filename": f"{pdb_code}.pdb",
            "message": f"Molecule {pdb_code} successfully fetched from PDB."
        }
    except FileNotFoundError:
        logger.error(f"PDB code {pdb_code} not found.")
        raise RemoteMoleculeNotFoundError(pdb_code)
    except Exception as e:
        logger.exception(f"Error fetching molecule {pdb_code} from external database: {e}")
        raise InternalError(f"Error fetching from PDB: {str(e)}")


@router.post("/remove-residue/{workspace_id}")
async def delete_residue(workspace_id: str, request: Request):
    workspace_manager.require_workspace(workspace_id)
    data = await request.json()

    FILENAME = "structure.pdb"
    pdb_path = workspace_manager.get_workspace_dir(workspace_id) / FILENAME

    logger.info(
        f"Removing residue {data.get('resseq')} from chain {data.get('chain')} "
        f"in workspace {workspace_id}..."
    )

    success = remove_residue_from_pdb(
        pdb_path=pdb_path,
        chain=data.get("chain"),
        resseq=int(data.get("resseq"))
    )

    if not success:
        logger.warning(f"Residue removal failed for in {workspace_id}")
        raise NotFoundError("Residue not found or could not be removed.")

    logger.info(f"Residue {data.get('resseq')} removed from chain {data.get('chain')} in workspace {workspace_id}")
    return {"status": "success"}