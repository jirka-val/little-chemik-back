"""
Wrapper nad data/force_fields.json - jediné místo v backendu, které tenhle
soubor čte i zapisuje.

force_fields.json rozděluje force fieldy do tierů (recommended / supported /
obsolete, plus nová new_unclassified) a definuje "profily" vody - svazek
konkrétního vodního modelu s kompatibilní sadou iontových FF pro I1/I1+/Im/
Im+ (viz solvent.profiles). Tenhle soubor je čistě deklarativní data; tahle
service nad ním poskytuje dotazy, které potřebuje /api/forcefields
(get_my_forcefields), a dvě zápisové operace:

  - reconcile(): po každém obnovení katalogu z IDA (FFCatalogService)
    zjistí, které FF force_fields.json vůbec nezná, a přidá je do
    new_unclassified - Pavel: nové FF se NEMAJÍ automaticky zařadit do
    "supported", mají počkat na ruční zařazení adminem.
  - set_solute_tier(): admin endpoint pro přeřazení D/R/P force fieldu do
    jiného tieru.
"""

import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from app.core.config import settings

logger = logging.getLogger(__name__)

SOLUTE_GROUPS = ("D", "R", "P")
ION_MOL_TYPES = ("I1", "I1+", "Im", "Im+")
WATER_MOL_TYPES = ("W3", "W4", "W5")
SOLUTE_TIERS = ("recommended", "supported", "obsolete")
TIER_RANK = {"recommended": 0, "supported": 1, "obsolete": 2, "new_unclassified": 3}


class ForceFieldClassificationService:
    def __init__(self, path: Optional[Path] = None):
        self.path = path or settings.FORCE_FIELDS_CLASSIFICATION_FILE
        self._data: Dict[str, Any] = {}
        self.reload()

    def reload(self) -> None:
        with open(self.path, "r", encoding="utf-8") as f:
            self._data = json.load(f)

    def _save(self) -> None:
        # Zápis přes dočasný soubor + replace, ať souběžný request nikdy
        # nepřečte napůl zapsaný JSON.
        tmp_path = self.path.with_suffix(".json.tmp")
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(self._data, f, indent=2, ensure_ascii=False)
            f.write("\n")
        tmp_path.replace(self.path)

    # ---------- solute (D/R/P) ----------

    def solute_group(self, group: str) -> Optional[Dict[str, Any]]:
        return self._data.get("solute_force_fields", {}).get(group)

    def classify_solute(self, group: str, ff_name: str) -> Dict[str, Any]:
        g = self.solute_group(group)
        if not g:
            return {"tier": "new_unclassified", "is_default": False}
        is_default = ff_name == g.get("default")
        for tier in SOLUTE_TIERS:
            if ff_name in (g.get(tier) or []):
                return {"tier": tier, "is_default": is_default}
        return {"tier": "new_unclassified", "is_default": False}

    # ---------- voda: profily ----------

    def water_types(self) -> Dict[str, Any]:
        return self._data.get("solvent", {}).get("water_types", {})

    def default_water_type(self) -> Optional[str]:
        return self._data.get("solvent", {}).get("default_water_type")

    def profiles(self) -> Dict[str, Any]:
        return self._data.get("solvent", {}).get("profiles", {})

    def default_profile_id(self) -> Optional[str]:
        wt = self.default_water_type()
        if not wt:
            return None
        return self.water_types().get(wt, {}).get("default_profile")

    def water_profiles_for(self, water_type: str) -> List[Dict[str, Any]]:
        """Profily patřící k danému vodnímu typu (W3/W4/W5), obohacené o
        tier/is_default - tohle je jednotka, kterou FF panel nabízí místo
        surového seznamu vodních FF (viz plán, sekce Profily vody a iontů).

        `is_default` je JEDINÝ globální default (solvent.default_water_type ->
        jeho default_profile), ne "default v rámci tohoto vodního typu" -
        guided mód smí nabídnout přesně jeden profil, ne jeden za každý
        vodní typ. `is_default_for_water_type` nese tu slabší, per-typovou
        informaci pro případ, že by ji frontend chtěl použít při přepnutí
        na jiný vodní typ ve standard/expert módu."""
        wt_data = self.water_types().get(water_type)
        if not wt_data:
            return []
        local_default_profile_id = wt_data.get("default_profile")
        global_default_profile_id = self.default_profile_id()
        result = []
        for profile_id in wt_data.get("profiles", []):
            profile = self.profiles().get(profile_id)
            if not profile:
                continue
            result.append({
                "profile_id": profile_id,
                "water_type": water_type,
                "tier": wt_data.get("tier", "supported"),
                "is_default": profile_id == global_default_profile_id,
                "is_default_for_water_type": profile_id == local_default_profile_id,
                "water": profile.get("water"),
                "ions": profile.get("ions", {}),
                "compatibility": profile.get("compatibility"),
            })
        return result

    def all_water_profiles(self) -> List[Dict[str, Any]]:
        result: List[Dict[str, Any]] = []
        for water_type in self.water_types():
            result.extend(self.water_profiles_for(water_type))
        return result

    # ---------- ionty: tier odvozený z profilů ----------

    def ion_tier_map(self, mol_type: str) -> Dict[str, Dict[str, Any]]:
        """
        Ionty nemají vlastní recommended/supported/obsolete seznam - jejich
        tier se odvozuje z toho, v jak dobrých profilech se používají (nejlepší
        tier napříč všemi profily, které daný iontový FF pro `mol_type`
        používají). is_default = True, pokud je součástí default profilu
        default vodního typu (solvent.default_water_type).
        """
        default_profile_id = self.default_profile_id()
        best_tier: Dict[str, str] = {}
        is_default: Dict[str, bool] = {}

        for water_type, wt_data in self.water_types().items():
            tier = wt_data.get("tier", "supported")
            for profile_id in wt_data.get("profiles", []):
                profile = self.profiles().get(profile_id)
                if not profile:
                    continue
                ion_name = (profile.get("ions") or {}).get(mol_type)
                if not ion_name:
                    continue
                current = best_tier.get(ion_name)
                if current is None or TIER_RANK[tier] < TIER_RANK[current]:
                    best_tier[ion_name] = tier
                if profile_id == default_profile_id:
                    is_default[ion_name] = True

        return {
            name: {"tier": tier, "is_default": is_default.get(name, False)}
            for name, tier in best_tier.items()
        }

    def classify_ion(self, mol_type: str, ff_name: str) -> Dict[str, Any]:
        entry = self.ion_tier_map(mol_type).get(ff_name)
        return entry or {"tier": "new_unclassified", "is_default": False}

    # ---------- reconciliace nových FF z IDA ----------

    def reconcile(self, catalog: List[Dict[str, Any]], ff_name_fn: Callable[[Dict[str, Any]], str]) -> int:
        """
        Projde stažený katalog a pro každou skupinu (D/R/P i vodní/iontové
        mol_type) zaznamená jména FF, která force_fields.json vůbec nezná, do
        new_unclassified - místo aby automaticky spadly do "supported".
        Zapíše na disk jen pokud se něco změnilo. Vrací počet nově přidaných
        záznamů (pro log FFCatalogService.refresh_catalog()).
        """
        added = 0
        solute = self._data.setdefault("solute_force_fields", {})
        solvent = self._data.setdefault("solvent", {})
        solvent_unclassified = solvent.setdefault("new_unclassified", {})
        ion_known_cache = {mt: set(self.ion_tier_map(mt).keys()) for mt in ION_MOL_TYPES}
        water_known_cache = {
            mt: {
                self.profiles().get(pid, {}).get("water")
                for pid in self.water_types().get(mt, {}).get("profiles", [])
            }
            for mt in WATER_MOL_TYPES
        }

        for ff in catalog:
            name = ff_name_fn(ff)
            for mol_type in ff.get("molecule_type") or []:
                if mol_type in SOLUTE_GROUPS:
                    group = solute.get(mol_type)
                    if not group:
                        continue
                    known = (
                        set(group.get("recommended", []))
                        | set(group.get("supported", []))
                        | set(group.get("obsolete", []))
                        | set(group.get("new_unclassified", []))
                    )
                    if name not in known:
                        group.setdefault("new_unclassified", []).append(name)
                        added += 1
                elif mol_type in ION_MOL_TYPES:
                    known = ion_known_cache[mol_type] | set(solvent_unclassified.get(mol_type, []))
                    if name not in known:
                        solvent_unclassified.setdefault(mol_type, []).append(name)
                        added += 1
                elif mol_type in WATER_MOL_TYPES:
                    known = water_known_cache[mol_type] | set(solvent_unclassified.get(mol_type, []))
                    if name not in known:
                        solvent_unclassified.setdefault(mol_type, []).append(name)
                        added += 1

        if added:
            self._save()
            logger.info(f"FF classification reconcile: {added} new unclassified force field(s) recorded.")
        return added

    # ---------- admin: přeřazení tieru ----------

    def set_solute_tier(self, group: str, ff_name: str, tier: str) -> None:
        """Přeřadí FF do jiného tieru ve skupině D/R/P (admin PATCH endpoint).

        Voda/ionty se v první verzi řadí ručně přímo v force_fields.json -
        přeřazení profilu je kurátorské rozhodnutí (jaká sada iontů k vodě
        patří), ne prosté přesunutí jména mezi seznamy, takže si nezaslouží
        stejné zjednodušené API."""
        if tier not in SOLUTE_TIERS:
            raise ValueError(f"Unknown tier '{tier}', expected one of {SOLUTE_TIERS}")
        g = self.solute_group(group)
        if g is None:
            raise ValueError(f"Unknown solute group '{group}', expected one of {SOLUTE_GROUPS}")

        for bucket in (*SOLUTE_TIERS, "new_unclassified"):
            if ff_name in (g.get(bucket) or []):
                g[bucket].remove(ff_name)

        g.setdefault(tier, []).append(ff_name)
        self._save()
        logger.info(f"FF classification: '{ff_name}' reclassified to '{tier}' in group '{group}'.")


classification_service = ForceFieldClassificationService()
