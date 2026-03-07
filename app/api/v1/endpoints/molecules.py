import logging
from fastapi import APIRouter, HTTPException, File, UploadFile
from app.services.pdb_service import PDBService
from app.workspaces.manager import workspace_manager
from fastapi import Body
from app.services.structure.hydrogenation import HydrogenationService

# Obnovení vašeho loggeru
logger = logging.getLogger("api")

router = APIRouter()
pdb_service = PDBService()

# Vytvoření instance služby (pod router = APIRouter() a pdb_service = ...)
hydrogen_service = HydrogenationService()

@router.post("/upload")
async def upload_molecule(file: UploadFile = File(...)):
    """
    Přijme skutečný soubor z frontendu (multipart/form-data),
    vytvoří pro něj workspace a vrátí jeho ID.
    """
    if not file.filename.endswith('.pdb'):
        logger.warning(f"Uživatel se pokusil nahrát nepodporovaný formát: {file.filename}")
        raise HTTPException(status_code=400, detail="Zatím podporujeme pouze .pdb soubory")

    try:
        workspace_id = await workspace_manager.create_from_upload(file)
        logger.info(f"Úspěšně vytvořen workspace {workspace_id} ze souboru {file.filename}")

        return {
            "workspace_id": workspace_id,
            "filename": file.filename,
            "message": "Molekula úspěšně nahrána a Workspace vytvořen."
        }
    except Exception as e:
        logger.exception(f"Kritická chyba při ukládání souboru {file.filename}: {e}")
        raise HTTPException(status_code=500, detail="Nepodařilo se uložit nahraný soubor.")


@router.get("/fetch-pdb/{pdb_code}")
async def fetch_pdb_by_code(pdb_code: str):
    """
    Stáhne molekulu z Protein Data Bank,
    rovnou ji uloží na disk do workspace a vrátí její ID.
    """
    try:
        logger.info(f"Požadavek na stažení PDB kódu: {pdb_code}")
        pdb_content = await pdb_service.get_remote_pdb_content(pdb_code.lower())

        workspace_id = workspace_manager.create_from_string(pdb_content)
        logger.info(f"Úspěšně staženo a vytvořen workspace {workspace_id} pro {pdb_code}")

        return {
            "workspace_id": workspace_id,
            "filename": f"{pdb_code}.pdb",
            "message": f"Molekula {pdb_code} úspěšně stažena z PDB."
        }
    except FileNotFoundError as e:
        logger.error(f"PDB kód {pdb_code} nebyl nalezen.")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception(f"Chyba při stahování molekuly {pdb_code} z externí databáze.")
        raise HTTPException(status_code=500, detail=f"Chyba při stahování z PDB: {str(e)}")


@router.post("/add-hydrogens/{workspace_id}")
async def add_hydrogens(
        workspace_id: str,
        ph: float = Body(7.0, embed=True)
):
    """
    Přečte molekulu z disku podle workspace_id, přidá vodíky podle zadaného pH,
    a přepíše dočasný soubor novou verzí.
    """
    logger.info(f"Požadavek na přidání vodíků pro workspace: {workspace_id} s pH {ph}")

    # 0. Zkontrolujeme, jestli soubor stále existuje
    if not workspace_manager.workspace_exists(workspace_id):
        logger.error(f"Workspace {workspace_id} nebyl nalezen.")
        raise HTTPException(status_code=404, detail="Workspace nenalezen. Nahráli jste soubor?")

    try:
        file_path = workspace_manager.get_file_path(workspace_id)

        # 1. Přečtení starého souboru z disku do paměti (tvůj PDBFixer potřebuje text)
        with open(file_path, "r", encoding="utf-8") as f:
            pdb_text = f.read()

        # 2. Zavoláme tvou chemickou službu (ta zůstala úplně stejná)
        updated_pdb_text = hydrogen_service.add_hydrogen_atoms(pdb_text, ph=ph)

        # 3. PŘEPÍŠEME ten samý soubor na disku novou, opravenou verzí
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(updated_pdb_text)

        logger.info(f"Vodíky úspěšně přidány do {workspace_id}.")

        # Vracíme pouze potvrzení! Žádný obří textový string.
        return {
            "workspace_id": workspace_id,
            "message": "Vodíky byly úspěšně doplněny.",
            "ph": ph
        }

    except Exception as e:
        logger.exception(f"Chyba při doplňování vodíků pro {workspace_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Chyba při úpravě struktury: {str(e)}")