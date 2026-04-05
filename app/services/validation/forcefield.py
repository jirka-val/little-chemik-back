import logging
import requests
from typing import List

logger = logging.getLogger(__name__)

class ForceFieldValidator:
    EXTERNAL_URL = "https://next.ida.4sims.eu/api/force_fields/"

    def __init__(self, dict_path: str = None):
        pass

    # ODSTRANILI JSME @staticmethod, aby to šlo volat přes instanci ff_validator
    def get_matching_forcefields(self, molecule_types: List[str]):
        """
        Stáhne FF z externího API a vybere ty, které odpovídají typům v PDB.
        """
        search_types = set(molecule_types)
        if "W" in search_types:
            search_types.update(["W3", "W4", "W5"])

        try:
            headers = {"x-client-version": "0.1.0"}
            response = requests.get(self.EXTERNAL_URL, headers=headers, timeout=10)
            response.raise_for_status()
            all_ffs = response.json()

            matched_ffs = []
            for ff in all_ffs:
                ff_types = ff.get("molecule_type") or []

                if any(t in search_types for t in ff_types):
                    matched_ffs.append(ff)

            return matched_ffs
        except Exception as e:
            logger.error(f"Chyba při komunikaci s FF API: {e}")
            return []