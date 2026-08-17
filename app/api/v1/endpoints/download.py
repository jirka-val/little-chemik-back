import io
import zipfile
import logging
from pathlib import Path
from fastapi import APIRouter, Query
from fastapi.responses import StreamingResponse, FileResponse
from pydantic import BaseModel

from app.core.exceptions import BadRequestError, NotFoundError
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()


# 1. Pydantic model pro přijetí dat z frontendu (checkboxy)
class DownloadRequest(BaseModel):
    wants_pdb: bool = False
    wants_top: bool = False
    wants_crd: bool = False
    as_zip: bool = True


# 2. Mapování formátů na konkrétní jména souborů
# Zkontroluj si, jestli se tvůj soubor se souřadnicemi jmenuje "structure.crd" nebo např. "structure.inpcrd"
FILE_MAPPING = {
    "pdb": {"filename": "structure.pdb", "media_type": "chemical/x-pdb"},
    "top": {"filename": "structure.prmtop", "media_type": "text/plain"},
    "crd": {"filename": "structure.crd", "media_type": "text/plain"}
}


@router.post("/{workspace_id}/export", summary="Export vybraných souborů (jednotlivě nebo jako ZIP)")
async def export_workspace_files(workspace_id: str, req: DownloadRequest):
    """
    Tento endpoint přijímá JSON z frontendu, kde uživatel vybral,
    které soubory chce stáhnout a zda je chce zabalit do ZIPu.
    """
    workspace_manager.require_workspace(workspace_id)

    files_to_pack = []

    # Získání bezpečných absolutních cest přes workspace_manager
    if req.wants_pdb:
        files_to_pack.append(workspace_manager.get_file_path(workspace_id, FILE_MAPPING["pdb"]["filename"]))
    if req.wants_top:
        files_to_pack.append(workspace_manager.get_file_path(workspace_id, FILE_MAPPING["top"]["filename"]))
    if req.wants_crd:
        files_to_pack.append(workspace_manager.get_file_path(workspace_id, FILE_MAPPING["crd"]["filename"]))

    # Filtrace pouze existujících souborů
    valid_files = [f for f in files_to_pack if f.exists()]

    if not valid_files:
        raise NotFoundError("Žádný z vybraných souborů nebyl ve workspace nalezen. Byla už vygenerována topologie?")

    # Pokud uživatel chce ZIP, NEBO vybral více souborů (přes HTTP nelze poslat více souborů najednou bez ZIPu)
    if req.as_zip or len(valid_files) > 1:
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for file_path in valid_files:
                # Zapisujeme pouze název souboru, ne celou cestu na serveru
                zip_file.write(file_path, file_path.name)

        # Vrácení ukazatele na začátek souboru, aby ho šlo přečíst
        zip_buffer.seek(0)

        logger.info(f"Serving ZIP archive for workspace: {workspace_id} with files: {[f.name for f in valid_files]}")
        return StreamingResponse(
            zip_buffer,
            media_type="application/x-zip-compressed",
            headers={
                "Content-Disposition": f"attachment; filename=little_chemik_export_{workspace_id}.zip"
            }
        )
    else:
        # Uživatel vybral jen JEDEN soubor a cíleně NECHCE zip -> pošleme napřímo
        single_file = valid_files[0]
        logger.info(f"Serving single file {single_file.name} for workspace: {workspace_id}")

        # Určení správného media type
        media_type = "application/octet-stream"
        for key, val in FILE_MAPPING.items():
            if val["filename"] == single_file.name:
                media_type = val["media_type"]
                break

        return FileResponse(
            path=single_file,
            media_type=media_type,
            filename=f"{workspace_id}_{single_file.name}"
        )


@router.get("/{workspace_id}", summary="Starý endpoint pro stažení jednoho souboru (zpětná kompatibilita)")
async def download_workspace_single_file(
        workspace_id: str,
        format: str = Query("pdb", description="Formát (pdb, crd, top)")
):
    """
    Tento endpoint zachováváme, protože frontendové knihovny (jako Molstar)
    potřebují jednoduchou GET URL, aby si mohly zobrazit strukturu.
    """
    workspace_manager.require_workspace(workspace_id)

    if format not in FILE_MAPPING:
        raise BadRequestError(f"Unsupported format: {format}")

    file_config = FILE_MAPPING[format]
    file_path = workspace_manager.get_file_path(workspace_id, file_config["filename"])

    if not file_path.exists():
        raise NotFoundError(f"File {file_config['filename']} not found.")

    logger.info(f"Serving {format.upper()} file via GET for workspace: {workspace_id}")
    return FileResponse(
        path=file_path,
        media_type=file_config["media_type"],
        filename=f"{workspace_id}_{file_config['filename']}"
    )