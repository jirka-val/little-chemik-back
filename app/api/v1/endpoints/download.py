import io
import re
import zipfile
import logging
from pathlib import Path
from typing import Optional
from fastapi import APIRouter, Query
from fastapi.responses import StreamingResponse, FileResponse, PlainTextResponse
from pydantic import BaseModel

from app.core.exceptions import BadRequestError, NotFoundError
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()

# Bezpečný "holý" název souboru pro explicitní ?filename= override níže -
# workspace_manager.get_file_path() dělá jen workspace_dir / filename beze
# sanitizace, takže tohle je jediná obrana proti path traversal (../..),
# jakmile filename poprvé přišel z query stringu (dřív šel jen přes pevný
# FILE_MAPPING, kde uživatelský vstup do cesty vůbec nevstupoval).
_SAFE_FILENAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*\.pdb$")

# Rezidua objemového rozpouštědla, která interaktivní 3D viewer (Molstar)
# nepotřebuje vidět atom po atomu - viz "light" parametr níže.
_BULK_SOLVENT_RESNAMES = {"WAT", "HOH", "SOL"}


def _strip_bulk_solvent(pdb_text: str) -> str:
    """
    Odstraní ATOM/HETATM záznamy objemové vody (WAT/HOH/SOL) z PDB textu.
    Používá se výhradně pro odlehčené načtení do interaktivního 3D vieweru u
    velkých solvatovaných struktur (u 1JJ2 to bylo ~780 000 z ~925 000 atomů,
    přes 84 % celého souboru) - ne pro export/stažení, kde uživatel čekává
    kompletní data. Ionty a krystalové ligandy zůstávají beze změny.
    """
    kept = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() in _BULK_SOLVENT_RESNAMES:
            continue
        kept.append(line)
    return "\n".join(kept)


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
        format: str = Query("pdb", description="Formát (pdb, crd, top)"),
        light: bool = Query(
            False,
            description="Pro pdb: vynechá objemovou vodu (WAT/HOH/SOL) - určeno pro "
                        "interaktivní 3D viewer u velkých solvatovaných struktur, ne pro export."
        ),
        filename: Optional[str] = Query(
            None,
            description="Explicitní přepis souboru uvnitř workspace (např. structure_preview.pdb "
                        "pro side-chain GUI náhled) místo defaultní cesty odvozené z 'format'."
        ),
):
    """
    Tento endpoint zachováváme, protože frontendové knihovny (jako Molstar)
    potřebují jednoduchou GET URL, aby si mohly zobrazit strukturu.
    """
    workspace_manager.require_workspace(workspace_id)

    if filename is not None:
        if not _SAFE_FILENAME_RE.match(filename):
            raise BadRequestError("Invalid filename.")
        file_path = workspace_manager.get_file_path(workspace_id, filename)
        media_type = "chemical/x-pdb"
    else:
        if format not in FILE_MAPPING:
            raise BadRequestError(f"Unsupported format: {format}")
        file_config = FILE_MAPPING[format]
        file_path = workspace_manager.get_file_path(workspace_id, file_config["filename"])
        media_type = file_config["media_type"]

    if not file_path.exists():
        raise NotFoundError(f"File {file_path.name} not found.")

    if file_path.suffix == ".pdb" and light:
        # Odlehčená verze se počítá za běhu (jeden průchod textem, řádově
        # desítky ms i na desítky MB) - nemá smysl ji cachovat na disk, "light"
        # je čistě prezentační ořez pro viewer, ne artefakt výpočtu.
        pdb_text = file_path.read_text(encoding="utf-8")
        stripped = _strip_bulk_solvent(pdb_text)
        logger.info(
            f"Serving lightweight (solvent-stripped) PDB via GET for workspace: {workspace_id} "
            f"({len(pdb_text.splitlines())} -> {len(stripped.splitlines())} lines)"
        )
        return PlainTextResponse(content=stripped, media_type=media_type)

    logger.info(f"Serving {file_path.name} via GET for workspace: {workspace_id}")
    return FileResponse(
        path=file_path,
        media_type=media_type,
        filename=f"{workspace_id}_{file_path.name}"
    )