# app/services/pdb_service.py
import httpx
import logging
from app.core.config import settings

logger = logging.getLogger(__name__)


class PDBService:
    async def get_remote_pdb_content(self, pdb_code: str) -> str:
        logger.info(f"Stahuji molekulu {pdb_code} z RCSB pro frontend...")

        async with httpx.AsyncClient() as client:
            response = await client.get(f"{settings.RCSB_PDB_URL}/{pdb_code}.pdb")

            if response.status_code != 200:
                raise FileNotFoundError(f"Molekula {pdb_code} neexistuje v RCSB databázi.")

            return response.text