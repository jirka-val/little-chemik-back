"""
Server-side snapshot katalogu force fieldů stažených z IDA.

Dřív GET /api/forcefields/{workspace_id} volalo IDA synchronně při KAŽDÉM
požadavku (viz ForceFieldService.get_matching_forcefields) - to byl přímý
zdroj pomalosti FF panelu, na kterou upozornil Pavel. Tahle service místo
toho drží jeden perzistentní JSON snapshot na disku (data/ff_catalog.json);
request handler z něj vždy jen čte, a živé volání IDA dělá pouze
refresh_catalog() - buď ručně (tlačítko Refresh ve FF panelu), nebo z
nočního background jobu (app/workspaces/tasks/ff_catalog_refresher.py).
"""

import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.core.config import settings
from app.services.forcefield_service import ForceFieldService

logger = logging.getLogger(__name__)


class FFCatalogService:
    def __init__(self, path: Optional[Path] = None, ff_service: Optional[ForceFieldService] = None):
        self.path = path or settings.FF_CATALOG_SNAPSHOT_FILE
        self.ff_service = ff_service or ForceFieldService()

    def _load(self) -> Dict[str, Any]:
        if not self.path.exists():
            return {"fetched_at": None, "forcefields": []}
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.error(f"FF catalog snapshot at {self.path} is unreadable ({e}), treating as empty.")
            return {"fetched_at": None, "forcefields": []}

    def get_forcefields(self) -> List[Dict[str, Any]]:
        """Rychlá cesta pro request handlery - čte ze souboru, NIKDY nevolá IDA."""
        return self._load().get("forcefields", [])

    def fetched_at(self) -> Optional[str]:
        return self._load().get("fetched_at")

    def refresh_catalog(self) -> Dict[str, Any]:
        """
        Jediné místo, které smí zavolat IDA pro celý katalog. Zapíše nový
        snapshot atomicky (tmp soubor + replace) a spustí reconciliaci proti
        klasifikaci - nově objevené FF skončí v new_unclassified, ne rovnou
        v "supported" (viz ForceFieldClassificationService.reconcile).
        """
        from app.services.ff_classification_service import classification_service

        all_ffs = self.ff_service.fetch_all_forcefields()
        snapshot = {
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "forcefields": all_ffs,
        }

        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = self.path.with_suffix(".json.tmp")
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(snapshot, f, ensure_ascii=False)
        tmp_path.replace(self.path)

        added = classification_service.reconcile(all_ffs, self.ff_service.ff_name)
        logger.info(f"FF catalog refreshed: {len(all_ffs)} force field(s) from IDA, {added} newly unclassified.")
        return snapshot

    def ensure_catalog(self) -> Dict[str, Any]:
        """Bootstrap pro první spuštění po deploy - pokud na disku ještě není
        žádný snapshot, udělá jeden synchronní refresh, ať GET /api/forcefields
        nevrátí prázdno, než doběhne první noční job."""
        if self.path.exists():
            return self._load()
        return self.refresh_catalog()


catalog_service = FFCatalogService()
