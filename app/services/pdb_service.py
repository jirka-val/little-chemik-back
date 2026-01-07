import httpx
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class PDBService:
    DATA_DIR = Path("data/pdb_files")
    BASE_URL = "https://files.rcsb.org/download"

    def __init__(self):
        self.DATA_DIR.mkdir(parents=True, exist_ok=True)

    async def fetch_pdb_content(self, pdb_code: str) -> str:
        pdb_filename = self.DATA_DIR / f"{pdb_code}.pdb"

        if not pdb_filename.exists():
            logger.info(f"Stahuji molekulu {pdb_code} z RCSB...")
            async with httpx.AsyncClient() as client:
                response = await client.get(f"{self.BASE_URL}/{pdb_code}.pdb")
                if response.status_code != 200:
                    logger.error(f"Molekula {pdb_code} nebyla nalezena na RCSB serveru.")
                    raise FileNotFoundError(f"Molekula {pdb_code} neexistuje.")

                pdb_filename.write_text(response.text)

        return pdb_filename.read_text()