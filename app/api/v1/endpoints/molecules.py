import logging
import shutil
from typing import Optional

from fastapi import APIRouter, File, UploadFile, Request
from fastapi.concurrency import run_in_threadpool

from app.services import input_files
from app.services.report import records as report_records
from app.services.report.citation import parse_jrnl
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
    # .ent = PDB soubor pod jménem z archivu RCSB (pdbXXXX.ent), frontend ho nabízí
    if not file.filename.lower().endswith(('.pdb', '.ent')):
        logger.warning(f"Unsupported format attempt: {file.filename}")
        raise BadRequestError("Only .pdb and .ent files are currently supported.")

    try:
        workspace_id = await workspace_manager.create_from_upload(file)
        logger.info(f"Successfully created workspace {workspace_id} from {file.filename}")
        await run_in_threadpool(incident_history.record_workspace_created, workspace_id, f"upload:{file.filename}")
        pdb_text = workspace_manager.get_file_path(workspace_id).read_text(encoding="utf-8", errors="replace")
        report_records.save_record(workspace_id, report_records.SOURCE, {
            "kind": "upload", "filename": file.filename, "citation": parse_jrnl(pdb_text),
        })

        return {
            "workspace_id": workspace_id,
            "filename": file.filename,
            "message": "Molecule uploaded and workspace created successfully."
        }
    except Exception as e:
        logger.exception(f"Critical error saving file {file.filename}: {e}")
        raise InternalError("Failed to save the uploaded file.")


@router.post("/{workspace_id}/inputs", summary="Přidá k nahranému PDB topologii a/nebo souřadnice")
async def upload_inputs(
    workspace_id: str,
    topology: Optional[UploadFile] = File(None),
    coordinates: Optional[UploadFile] = File(None),
):
    """
    Topologie (.prmtop/.parm7/.top) a souřadnice (.crd/.inpcrd/.rst7/.rst/.ncrst)
    z upload dialogu. Zatím se jen uloží a vrátí se jejich souhrn - viz
    app/services/input_files.py.
    """
    workspace_manager.require_workspace(workspace_id)
    if topology is None and coordinates is None:
        raise BadRequestError("No topology or coordinate file was sent.", code="nothing_selected")
    for upload, extensions, kind in ((topology, input_files.TOPOLOGY_EXTENSIONS, "topology"),
                                     (coordinates, input_files.COORDINATE_EXTENSIONS, "coordinates")):
        if upload is not None and not upload.filename.lower().endswith(extensions):
            raise BadRequestError(
                f"\"{upload.filename}\" is not a supported {kind} file ({', '.join(extensions)}).",
                code=f"unsupported_{kind}_format",
            )

    workspace_dir = workspace_manager.get_workspace_dir(workspace_id)
    manifest = input_files.load_manifest(workspace_dir)
    if topology is not None:
        with open(workspace_dir / input_files.UPLOADED_TOPOLOGY, "wb") as fh:
            shutil.copyfileobj(topology.file, fh)
        manifest["topology"] = topology.filename
    if coordinates is not None:
        for old in workspace_dir.glob(input_files.UPLOADED_COORDINATES + ".*"):
            old.unlink()
        with open(workspace_dir / input_files.coordinates_filename(coordinates.filename), "wb") as fh:
            shutil.copyfileobj(coordinates.file, fh)
        manifest["coordinates"] = coordinates.filename
    input_files.save_manifest(workspace_dir, manifest)
    logger.info(f"Workspace {workspace_id}: uploaded inputs {manifest}")

    return await run_in_threadpool(input_files.summarize_inputs, workspace_dir)


@router.get("/{workspace_id}/inputs", summary="Souhrn nahraných vstupů (PDB, topologie, souřadnice)")
async def get_inputs_summary(workspace_id: str):
    workspace_manager.require_workspace(workspace_id)
    return await run_in_threadpool(input_files.summarize_inputs, workspace_manager.get_workspace_dir(workspace_id))


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
        report_records.save_record(workspace_id, report_records.SOURCE, {
            "kind": "pdb", "pdb_code": pdb_code, "citation": parse_jrnl(pdb_content),
        })

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