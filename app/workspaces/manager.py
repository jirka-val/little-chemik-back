import os
import uuid
import shutil
from fastapi import UploadFile

# Složka se vytvoří v hlavním adresáři projektu
WORKSPACE_DIR = os.path.join(os.getcwd(), "temp_workspaces")


class WorkspaceManager:
    def __init__(self):
        # Při startu aplikace (nebo importu) zkontroluje, zda složka existuje. Pokud ne, vytvoří ji.
        os.makedirs(WORKSPACE_DIR, exist_ok=True)

    async def create_from_upload(self, file: UploadFile) -> str:
        """Přijme soubor nahraný uživatelem, uloží ho a vrátí jeho UUID."""
        workspace_id = str(uuid.uuid4())
        file_path = self.get_file_path(workspace_id)

        # Uložení po částech (tzv. chunks) - nikdy nesežere paměť, ani u 1GB souboru
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)

        return workspace_id

    def create_from_string(self, content: str) -> str:
        """Přijme textový řetězec (např. stažený z Protein Data Bank), uloží ho a vrátí UUID."""
        workspace_id = str(uuid.uuid4())
        file_path = self.get_file_path(workspace_id)

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)

        return workspace_id

    def get_file_path(self, workspace_id: str) -> str:
        """Vrátí absolutní cestu k souboru podle jeho ID."""
        return os.path.join(WORKSPACE_DIR, f"{workspace_id}.pdb")

    def workspace_exists(self, workspace_id: str) -> bool:
        """Ověří, zda dočasný soubor stále existuje (zda už nebyl smazán)."""
        return os.path.exists(self.get_file_path(workspace_id))


# Vytvoříme jedinou instanci (tzv. Singleton), kterou budeme používat napříč aplikací
workspace_manager = WorkspaceManager()