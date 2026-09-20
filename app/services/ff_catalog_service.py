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

import base64
import json
import logging
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from app.core.config import settings
from app.services.forcefield_service import ForceFieldService

logger = logging.getLogger(__name__)

_ION_MOL_TYPES = ("I1", "I1+", "Im", "Im+")
_RESIDUE_SECTION_RE = re.compile(r"^\[\s*(\S+)\s*\]", re.MULTILINE)


class FFCatalogService:
    def __init__(self, path: Optional[Path] = None, ff_service: Optional[ForceFieldService] = None):
        self.path = path or settings.FF_CATALOG_SNAPSHOT_FILE
        self.ff_service = ff_service or ForceFieldService()
        # In-memory cache pro _load()/get_buildable_ion_resnames() - bez ní
        # se ff_catalog.json (~15 MB) četl a parsoval ze souboru při KAŽDÉM
        # volání get_forcefields() (naměřeno ~126 ms/volání), a
        # get_buildable_ion_resnames() nad tím navíc dekódovalo/regexovalo
        # base64 residue_lib obsah všech iontových FF znovu - u /prepare, kde
        # se volá 2x (cation+anion, viz validation.py _build_salt_specs), to
        # přidávalo ~260 ms na každý požadavek jen na opakované parsování
        # téhož souboru. `catalog_service` je proces-lokální singleton a
        # jediné místo, které soubor přepisuje, je refresh_catalog() (ruční
        # tlačítko i noční job běží nad stejnou instancí - viz
        # ff_catalog_refresher.py), takže cache je bezpečné invalidovat
        # výhradně tam.
        self._snapshot_cache: Optional[Dict[str, Any]] = None
        self._buildable_ions_cache: Optional[Dict[str, Set[str]]] = None

    def _load(self) -> Dict[str, Any]:
        if self._snapshot_cache is not None:
            return self._snapshot_cache
        if not self.path.exists():
            return {"fetched_at": None, "forcefields": []}
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                snapshot = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.error(f"FF catalog snapshot at {self.path} is unreadable ({e}), treating as empty.")
            return {"fetched_at": None, "forcefields": []}
        self._snapshot_cache = snapshot
        return snapshot

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

        self._snapshot_cache = snapshot
        self._buildable_ions_cache = None  # odvozené z obsahu, přepočítat líně při dalším dotazu

        added = classification_service.reconcile(all_ffs, self.ff_service.ff_name)
        logger.info(f"FF catalog refreshed: {len(all_ffs)} force field(s) from IDA, {added} newly unclassified.")
        return snapshot

    def get_buildable_ion_resnames(self) -> Dict[str, Set[str]]:
        """
        Pro každou iontovou mol_type skupinu (I1/I1+/Im/Im+) zjistí, které
        resnames mají v AKTUÁLNÍM katalogu reálně definované parametry
        (residue_lib_ff_file aspoň jednoho FF s tou skupinou), ne jen výskyt v
        converting_dictionary.json.

        Ten totiž zná chemickou identitu iontu (jaký mol_type by měl mít),
        ale ne, jestli pro něj v katalogu skutečně existuje FF - u Im+ je
        rozdíl reálný: ~8 iontů (Ca2+, Cd2+, Ce3+, Ce4+, Hg2+, U, U4+, V2+)
        converting_dictionary zná, ale žádný katalogový FF pro ně nemá
        parametry, takže by jejich výběr vždycky skončil KeyError hluboko v
        builderu. Navíc pro pár prvků katalog definuje jen starší
        dvoupísmenné jméno holého symbolu (CA/CD/CE/Ce/HG), ne
        nábojem-sufixovanou variantu (Ca2+/Cd2+/Ce3+.../Hg2+) - která z
        dvojice je reálně stavitelná, tak není možné odvodit jinak než
        přímým rozborem katalogu.

        Výsledek je cachovaný (base64 decode + regex nad ~60 FF souborů by se
        jinak opakovalo při každém volání - viz cache poznámka v __init__),
        invalidovaný jedině přes refresh_catalog(). Vrácený dict/set NENÍ
        kopie - volající ho čte, nikdy nemodifikuje (viz stávající použití v
        validation.py).
        """
        if self._buildable_ions_cache is not None:
            return self._buildable_ions_cache

        by_group: Dict[str, Set[str]] = {mt: set() for mt in _ION_MOL_TYPES}
        for ff in self.get_forcefields():
            mol_types = [mt for mt in (ff.get("molecule_type") or []) if mt in by_group]
            if not mol_types:
                continue
            raw = ff.get("residue_lib_ff_file")
            if not raw:
                continue
            try:
                content = base64.b64decode(raw).decode("utf-8", errors="replace")
            except Exception:
                logger.warning(f"FF '{ff.get('ff_name')}' has an undecodable residue_lib_ff_file, skipping.")
                continue
            defined = {m for m in _RESIDUE_SECTION_RE.findall(content) if m.lower() != "bondedtypes"}
            for mt in mol_types:
                by_group[mt] |= defined

        self._buildable_ions_cache = by_group
        return by_group

    def ensure_catalog(self) -> Dict[str, Any]:
        """Bootstrap pro první spuštění po deploy - pokud na disku ještě není
        žádný snapshot, udělá jeden synchronní refresh, ať GET /api/forcefields
        nevrátí prázdno, než doběhne první noční job."""
        if self.path.exists():
            return self._load()
        return self.refresh_catalog()


catalog_service = FFCatalogService()
