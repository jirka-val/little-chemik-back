import os
import uuid
import shutil
from pathlib import Path
from fastapi import UploadFile

from app.core.exceptions import WorkspaceNotFoundError

# Složka pro dočasné pracovní prostory
WORKSPACE_DIR = os.path.join(os.getcwd(), "temp_workspaces")


class WorkspaceManager:
    def __init__(self):
        # Zajistí existenci hlavní složky pro workspaces
        os.makedirs(WORKSPACE_DIR, exist_ok=True)

    def get_workspace_dir(self, workspace_id: str) -> Path:
        """Vrátí Path k adresáři konkrétního workspace a zajistí jeho existenci."""
        workspace_path = Path(WORKSPACE_DIR) / workspace_id
        workspace_path.mkdir(parents=True, exist_ok=True)
        return workspace_path

    async def create_from_upload(self, file: UploadFile) -> str:
        workspace_id = str(uuid.uuid4())
        workspace_path = self.get_workspace_dir(workspace_id)

        # aby ostatní služby věděly, co mají číst.
        file_path = workspace_path / "structure.pdb"

        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        return workspace_id

    def create_from_string(self, content: str, filename: str = "structure.pdb") -> str:
        """Vytvoří složku, uloží do ní řetězec jako soubor a vrátí UUID složky."""
        workspace_id = str(uuid.uuid4())
        workspace_path = self.get_workspace_dir(workspace_id)

        file_path = workspace_path / filename

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)

        return workspace_id

    # Přidáme defaultní hodnotu None nebo prázdný string
    def get_file_path(self, workspace_id: str, filename: str = "structure.pdb") -> Path:
        """Vrátí cestu k souboru. Pokud filename chybí, předpokládá structure.pdb."""
        return self.get_workspace_dir(workspace_id) / filename

    def workspace_exists(self, workspace_id: str) -> bool:
        """Ověří, zda adresář workspace existuje."""
        return os.path.exists(os.path.join(WORKSPACE_DIR, workspace_id))

    def require_workspace(self, workspace_id: str) -> None:
        """
        Vyhodí WorkspaceNotFoundError (jednotná 404 obálka), pokud workspace
        neexistuje. Nahrazuje opakované
        `if not workspace_manager.workspace_exists(...): raise HTTPException(404, ...)`
        rozeseté po endpointech vlastními, mírně odlišnými texty chyby.
        """
        if not self.workspace_exists(workspace_id):
            raise WorkspaceNotFoundError()


# Singleton instance
workspace_manager = WorkspaceManager()