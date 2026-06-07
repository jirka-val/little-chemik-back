# app/api/v1/endpoints/editor.py
import logging
import aiofiles
from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.services.editor_service import StructureEditorService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)
router = APIRouter()
editor_service = StructureEditorService()


# --- Modely pro příchozí data z frontendu ---

class MutateRequest(BaseModel):
    workspace_id: str
    chain_id: str
    residue_number: int
    mutate_to: str


class RenameResidueRequest(BaseModel):
    workspace_id: str
    chain_id: str
    residue_number: int
    new_res_name: str = Field(..., max_length=3)


class RenameAtomRequest(BaseModel):
    workspace_id: str
    chain_id: str
    residue_number: int
    old_atom_name: str
    new_atom_name: str = Field(..., max_length=4)


class RemoveAtomRequest(BaseModel):
    workspace_id: str
    chain_id: str
    residue_number: int
    atom_name: str


# --- Pomocná funkce pro čtení/zápis ---

async def _process_editor_action(workspace_id: str, action_func, *args):
    """Načte PDB, provede úpravu a uloží ho zpět."""
    try:
        pdb_path = workspace_manager.get_file_path(workspace_id)

        # Načtení současného stavu
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Provedení akce v threadpoolu (aby neblokovalo asynchronní smyčku)
        new_pdb_content = await run_in_threadpool(action_func, pdb_content, *args)

        if not new_pdb_content:
            raise ValueError("Editor service returned empty content.")

        # Uložení nového stavu
        async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
            await f.write(new_pdb_content)

        return {"message": "Success"}

    except Exception as e:
        logger.error(f"Editor action failed for workspace {workspace_id}: {str(e)}")
        raise HTTPException(status_code=500, detail={"error": str(e)})


# --- Endpointy ---

@router.post("/rename-residue")
async def rename_residue(request: RenameResidueRequest):
    logger.info(f"Renaming residue {request.residue_number} ({request.chain_id}) to {request.new_res_name}")
    return await _process_editor_action(
        request.workspace_id,
        editor_service.rename_residue,
        request.chain_id,
        request.residue_number,
        request.new_res_name
    )


@router.post("/mutate")
async def mutate_residue(request: MutateRequest):
    logger.info(f"Mutating residue {request.residue_number} ({request.chain_id}) to {request.mutate_to}")
    return await _process_editor_action(
        request.workspace_id,
        editor_service.mutate_residue,
        request.chain_id,
        request.residue_number,
        request.mutate_to
    )


@router.post("/rename-atom")
async def rename_atom(request: RenameAtomRequest):
    return await _process_editor_action(
        request.workspace_id,
        editor_service.rename_atom,
        request.chain_id,
        request.residue_number,
        request.old_atom_name,
        request.new_atom_name
    )


@router.post("/remove-atom")
async def remove_atom(request: RemoveAtomRequest):
    return await _process_editor_action(
        request.workspace_id,
        editor_service.remove_atom,
        request.chain_id,
        request.residue_number,
        request.atom_name
    )