import logging
import requests
import os
import shutil
from typing import List, Dict, Any
from pathlib import Path

logger = logging.getLogger(__name__)


class ForceFieldService:
    # URL na externí databázi silových polí
    EXTERNAL_URL = "https://next.ida.4sims.eu/api/force_fields/"

    # Lokální cache pro uložení extrahovaných souborů (.rtp, .itp, .atp)
    CACHE_DIR = Path("data/ff_cache")

    def __init__(self):
        # Vytvoření cache adresáře při inicializaci
        try:
            self.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.error(f"Failed to create cache directory {self.CACHE_DIR}: {e}")

    def get_matching_forcefields(self, molecule_types: List[str]) -> List[Dict[str, Any]]:
        """
        Stáhne kompletní seznam FF z API a profiltruje ty, které odpovídají typům v molekule.
        Slouží pro nabídku uživateli v UI.
        """
        search_types = set(molecule_types)
        # Rozšíření pro různé typy vody, pokud je přítomna
        if "W" in search_types:
            search_types.update(["W3", "W4", "W5"])

        try:
            headers = {"x-client-version": "0.1.0"}
            logger.info(f"Fetching force fields from {self.EXTERNAL_URL}")

            response = requests.get(self.EXTERNAL_URL, headers=headers, timeout=10)
            response.raise_for_status()
            all_ffs = response.json()

            matched_ffs = []
            for ff in all_ffs:
                ff_types = ff.get("molecule_type") or []
                if any(t in search_types for t in ff_types):
                    matched_ffs.append(ff)

            logger.info(f"Successfully matched {len(matched_ffs)} force fields.")
            return matched_ffs

        except requests.exceptions.RequestException as e:
            logger.error(f"External API communication error: {e}")
            return []
        except Exception as e:
            logger.error(f"Unexpected error while matching force fields: {e}")
            return []

    def prepare_forcefield_files(self, ff_data: Dict[str, Any]) -> Path:
        """
        Extrahuje soubory z obřího JSONu a uloží je na disk pro FF_IDA.py.
        Vrací Path k adresáři s konkrétním silovým polem.
        """
        ff_name = ff_data.get("name", "unknown_ff").replace(" ", "_")
        target_dir = self.CACHE_DIR / ff_name

        # Kontrola, zda už FF v cache existuje (abychom neparsovali znovu)
        if target_dir.exists() and any(target_dir.iterdir()):
            logger.info(f"Force field '{ff_name}' found in local cache.")
            return target_dir

        try:
            target_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Extracting and saving files for force field: {ff_name}")

            # Mapování klíčů z vašeho "extrémního" JSONu na fyzické soubory
            file_map = {
                "rtp_file": f"{ff_name}.rtp",
                "atp_file": f"{ff_name}.atp",
                "bonded_itp": f"bonded_{ff_name}.itp",
                "nonbonded_itp": f"nonbonded_{ff_name}.itp"
            }

            files_saved = 0
            for json_key, filename in file_map.items():
                content = ff_data.get(json_key)
                if content:
                    file_path = target_dir / filename
                    with open(file_path, "w", encoding="utf-8") as f:
                        f.write(content)
                    files_saved += 1

            if files_saved == 0:
                logger.warning(f"No valid force field files found in JSON for '{ff_name}'")
            else:
                logger.info(f"Saved {files_saved} files to {target_dir}")

            return target_dir

        except Exception as e:
            logger.error(f"Failed to prepare force field files for '{ff_name}': {e}")
            raise

    def clear_cache(self):
        """Vymaže celou složku ff_cache."""
        if self.CACHE_DIR.exists():
            try:
                shutil.rmtree(self.CACHE_DIR)
                self.CACHE_DIR.mkdir()
                logger.info("Force field cache cleared successfully.")
            except Exception as e:
                logger.error(f"Failed to clear force field cache: {e}")