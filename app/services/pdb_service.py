import httpx
from pathlib import Path
import logging
from app.core.config import settings

logger = logging.getLogger(__name__)

class PDBService:
    def __init__(self):
        self.data_dir = settings.PDB_DATA_DIR
        self.data_dir.mkdir(parents=True, exist_ok=True)

    async def fetch_pdb_content(self, pdb_code: str) -> str:
        pdb_filename = self.data_dir / f"{pdb_code}.pdb"

        if not pdb_filename.exists():
            logger.info(f"Stahuji molekulu {pdb_code}...")
            async with httpx.AsyncClient() as client:
                response = await client.get(f"{settings.RCSB_PDB_URL}/{pdb_code}.pdb")
                if response.status_code != 200:
                    raise FileNotFoundError(f"Molekula {pdb_code} neexistuje.")

                pdb_filename.write_text(response.text)

        return pdb_filename.read_text()