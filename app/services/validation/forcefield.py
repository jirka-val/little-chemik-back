import logging
import requests
from typing import List

logger = logging.getLogger(__name__)

class ForceFieldValidator:
    EXTERNAL_URL = "https://next.ida.4sims.eu/api/force_fields/"

    def __init__(self, dict_path: str = None):
        # Tvůj inicializační kód (načítání slovníku) může zůstat,
        # ale pro tuhle funkci ho teď nepotřebujeme
        pass

    # ODSTRANILI JSME @staticmethod, aby to šlo volat přes instanci ff_validator
    def get_matching_forcefields(self, molecule_types: List[str]):
        """
        Stáhne FF z externího API a vybere ty, které odpovídají typům v PDB.
        """
        # --- OPRAVA PRO VODU ---
        # Převedeme seznam na Set (množinu) a pokud je tam 'W',
        # rovnou k němu přihodíme i specifické varianty vody.
        search_types = set(molecule_types)
        if "W" in search_types:
            search_types.update(["W3", "W4", "W5"])
        # -----------------------

        try:
            # Přidáváme x-client-version header pro jistotu
            headers = {"x-client-version": "0.1.0"}
            response = requests.get(self.EXTERNAL_URL, headers=headers, timeout=10)
            response.raise_for_status()
            all_ffs = response.json()

            matched_ffs = []
            for ff in all_ffs:
                ff_types = ff.get("molecule_type") or []

                # Porovnáváme proti 'search_types', které teď obsahuje i W3, W4, W5
                if any(t in search_types for t in ff_types):
                    matched_ffs.append(ff)

            return matched_ffs
        except Exception as e:
            logger.error(f"Chyba při komunikaci s FF API: {e}")
            return []