import logging
import os
import shutil
from tempfile import NamedTemporaryFile
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse

from app.workspaces.manager import workspace_manager
from app.services.export_service import export_service

logger = logging.getLogger(__name__)

router = APIRouter()

SUPPORTED_FORMATS = {
    "pdb": {"filename": "structure.pdb", "media_type": "chemical/x-pdb"},
    "crd": {"filename": "structure.crd", "media_type": "text/plain"},
    "top": {"filename": "structure.prmtop", "media_type": "text/plain"},
}

@router.get("/{workspace_id}", summary="Stažení PDB, CRD, TOP, ZIP nebo konkrétního souboru")
async def download_workspace(
        workspace_id: str,
        # Default necháme pro zpětnou kompatibilitu, kdyby to nějaká stará část kódu používala
        format: str = Query("pdb", description="Formát (pdb, crd, top, zip)"),
        # PŘIDÁNO: Nový parametr filename, který teď posílá tvůj frontend
        filename: str = Query(None, description="Konkrétní jméno souboru ke stažení")
):
    try:
        if not workspace_manager.workspace_exists(workspace_id):
            logger.warning(f"Download attempt for non-existent workspace: {workspace_id}")
            raise HTTPException(status_code=404, detail="Workspace not found or expired.")

        # --- PŘIDÁNO: Zpracování parametru filename (Priorita) ---
        if filename:
            # DŮLEŽITÉ: Ochrana proti Path Traversal (aby někdo nestáhl např. ../../../etc/passwd)
            safe_filename = os.path.basename(filename)
            file_path = workspace_manager.get_file_path(workspace_id, safe_filename)

            if not file_path.exists():
                raise HTTPException(
                    status_code=404,
                    detail=f"File {safe_filename} not found in workspace."
                )

            logger.info(f"Serving explicit file {safe_filename} for workspace: {workspace_id}")
            return FileResponse(
                path=file_path,
                filename=f"{workspace_id}_{safe_filename}"
            )
        # ---------------------------------------------------------

        # --- Původní logika pro parametr format ---
        crd_path = workspace_manager.get_file_path(workspace_id, "structure.crd")
        if format in ["crd", "zip"] and not crd_path.exists():
            logger.info(f"CRD file missing for workspace {workspace_id}. Generating on demand...")
            success = export_service.generate_amber_crd_from_pdb(workspace_id)
            if not success and format == "crd":
                raise HTTPException(status_code=500, detail="Failed to generate CRD file from PDB.")

        if format == "zip":
            workspace_dir = workspace_manager.get_workspace_dir(workspace_id)
            tmp_zip = NamedTemporaryFile(delete=False, suffix=".zip")
            base_name = tmp_zip.name.replace('.zip', '')
            shutil.make_archive(base_name, 'zip', workspace_dir)

            logger.info(f"Serving ZIP archive for workspace: {workspace_id}")
            return FileResponse(
                path=tmp_zip.name,
                media_type="application/zip",
                filename=f"{workspace_id}_all_files.zip"
            )

        if format not in SUPPORTED_FORMATS:
            raise HTTPException(status_code=400, detail=f"Unsupported format: {format}")

        file_config = SUPPORTED_FORMATS[format]
        file_path = workspace_manager.get_file_path(workspace_id, file_config["filename"])

        if not file_path.exists():
            raise HTTPException(status_code=404,
                                detail=f"File {file_config['filename']} not found. Has topology been generated?")

        logger.info(f"Serving {format.upper()} file for workspace: {workspace_id}")
        return FileResponse(
            path=file_path,
            media_type=file_config["media_type"],
            filename=f"{workspace_id}_{file_config['filename']}"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Error serving file for workspace {workspace_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")