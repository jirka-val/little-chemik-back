import logging
import io
import json
from typing import Any, Dict
import aiofiles
from fastapi import APIRouter, File, UploadFile, Body, Request
from fastapi.concurrency import run_in_threadpool

from app.core.exceptions import AppBaseException, BadRequestError, InternalError, NotFoundError, RemoteMoleculeNotFoundError
from app.services.pdb_service import PDBService, remove_residue_from_pdb
from app.workspaces.manager import workspace_manager
from app.services.structure.forge_service import ForgeStructureService

logger = logging.getLogger("api")
router = APIRouter()
pdb_service = PDBService()
forge_service = ForgeStructureService()

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


@router.post("/add-hydrogens/{workspace_id}", deprecated=True,
             summary="[Deprecated] Použij POST /api/validation/prepare (add_solvent=False)")
async def add_hydrogens(
        workspace_id: str,
        ff_selections: Dict[str, Any] = Body(..., description="mol_type -> FF metadata, stejný tvar jako /prepare"),
        ph: float = Body(7.0, embed=True),
        optimize: bool = Body(False, embed=True)
):
    """
    Deprecated tenký alias nad ForgeStructureService.prepare_structure(add_solvent_and_ions=False) -
    stejná logika jako /api/validation/prepare, jen bez solvatace/iontů. `optimize` je zachováno
    v kontraktu z historických důvodů, builder žádnou samostatnou optimalizaci navíc neprovádí.
    """
    logger.info(f"Hydrogenation request for workspace: {workspace_id} (pH: {ph}, optimize: {optimize})")

    workspace_manager.require_workspace(workspace_id)

    try:
        file_path = workspace_manager.get_file_path(workspace_id, "structure.pdb")

        # Asynchronní I/O pro čtení
        async with aiofiles.open(file_path, "r", encoding="utf-8") as f:
            pdb_text = await f.read()

        # Delegace těžkého chemického výpočtu na vedlejší vlákno
        result = await run_in_threadpool(
            forge_service.prepare_structure,
            pdb_text=pdb_text,
            ff_selections=ff_selections,
            ph=ph,
            add_solvent_and_ions=False,
        )

        # Asynchronní I/O pro zápis výsledku
        async with aiofiles.open(file_path, "w", encoding="utf-8") as f:
            await f.write(result.pdb_text)

        meta_path = file_path.with_name("structure.forge_meta.json")
        async with aiofiles.open(meta_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(result.forge_meta))

        logger.info(f"Hydrogens added to {workspace_id}.")

        return {
            "workspace_id": workspace_id,
            "message": "Hydrogens successfully added.",
            "warnings": result.warnings,
            "ph": ph,
        }

    except AppBaseException:
        # ForgeMissingDOFError apod. - necháme propadnout ke globálnímu handleru
        # se svým vlastním status_code, ne zabalit do generické 500 níže.
        raise
    except Exception as e:
        logger.exception(f"Error during hydrogenation for {workspace_id}: {e}")
        raise InternalError(f"Structure modification failed: {str(e)}")


@router.post("/remove-residue/{workspace_id}")
async def delete_residue(workspace_id: str, request: Request):
    data = await request.json()

    FILENAME = "structure.pdb"
    pdb_path = workspace_manager.get_workspace_dir(workspace_id) / FILENAME

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