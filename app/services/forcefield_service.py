import logging
import requests
import shutil
import base64
from typing import List, Dict, Any
from pathlib import Path

logger = logging.getLogger(__name__)


class ForceFieldService:
    # URL na externí databázi silových polí
    EXTERNAL_URL = "https://next.ida.4sims.eu/api/force_fields/"

    # Lokální cache pro uložení extrahovaných souborů
    CACHE_DIR = Path("data/ff_cache")

    def __init__(self):
        try:
            self.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.error(f"Failed to create cache directory {self.CACHE_DIR}: {e}")

    def _decode_content(self, content: Any) -> str:
        """Bezpečně dekóduje obsah z API (Base64 nebo čistý text)."""
        if not content or not isinstance(content, str):
            return ""

        # Odstranění data URI prefixu, pokud existuje
        if "base64," in content:
            content = content.split("base64,")[1]

        try:
            # Pokus o dekódování base64
            decoded_bytes = base64.b64decode(content)
            return decoded_bytes.decode('utf-8')
        except Exception:
            # Pokud dekódování selže, předpokládáme, že je to už čistý text
            return content

    def prepare_forcefield_files(self, ff_data: Dict[str, Any]) -> Path:
        """
        Extrahuje, zploští a uloží soubory silového pole na disk.
        Spojuje hlavní RTP soubor s knihovnou reziduí a nahrazuje #include ostrými daty.
        """
        raw_name = ff_data.get('display_name') or ff_data.get('ff_name') or 'unknown_ff'
        ff_name = raw_name.replace(" ", "_")
        target_dir = self.CACHE_DIR / ff_name
        target_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Preparing force field: {ff_name}")

        # Načtení všech potřebných dat z API
        rtp_raw = ff_data.get('force_field_file')
        res_lib_raw = ff_data.get('residue_lib_ff_file')
        nb_raw = ff_data.get('nonbonded_ff_file')
        b_raw = ff_data.get('bonded_ff_file')
        atp_raw = ff_data.get('atom_type_ff_file')

        # Dekódování obsahu do textu
        rtp_content = self._decode_content(rtp_raw)
        res_lib_content = self._decode_content(res_lib_raw)
        nb_content = self._decode_content(nb_raw)
        b_content = self._decode_content(b_raw)
        atp_content = self._decode_content(atp_raw)

        if not res_lib_content:
            logger.warning(f"Residue library (residue_lib_ff_file) is missing for {ff_name}. This will likely cause KeyError.")

        # ZPLOŠTĚNÍ RTP (Kombinace a nahrazení #include ostrými daty z ITP)
        # Spojíme hlavní RTP wrapper s knihovnou reziduí (to je to "maso" silového pole).
        raw_combined_content = rtp_content + "\n" + res_lib_content

        flattened_rtp_lines = []
        for line in raw_combined_content.splitlines():
            clean_line = line.strip()
            if clean_line.lower().startswith('#include'):
                if "nonbonded" in clean_line.lower():
                    flattened_rtp_lines.append(nb_content)
                elif "bonded" in clean_line.lower():
                    flattened_rtp_lines.append(b_content)
            else:
                flattened_rtp_lines.append(line)

        final_rtp_content = "\n".join(flattened_rtp_lines)

        # Zápis finálních souborů na disk (názvy vyžadované knihovnou FF_IDA)
        files = {
            f"{ff_name}.rtp": final_rtp_content,
            f"nonbonded_{ff_name}.itp": nb_content,
            f"bonded_{ff_name}.itp": b_content,
            f"{ff_name}.atp": atp_content
        }

        for fname, data in files.items():
            with open(target_dir / fname, "w", encoding="utf-8") as f:
                f.write(data)

        logger.info(f"Force field {ff_name} successfully prepared and saved to disk.")
        return target_dir

    def get_matching_forcefields(self, molecule_types: List[str]) -> List[Dict[str, Any]]:
        """Stáhne seznam FF z API a vyfiltruje ty odpovídající molekule."""
        search_types = set(molecule_types)
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
            return matched_ffs
        except Exception as e:
            logger.error(f"External API communication error: {e}")
            return []

    def clear_cache(self):
        """Vymaže celou složku ff_cache."""
        if self.CACHE_DIR.exists():
            try:
                shutil.rmtree(self.CACHE_DIR)
                self.CACHE_DIR.mkdir()
                logger.info("Force field cache cleared successfully.")
            except Exception as e:
                logger.error(f"Failed to clear force field cache: {e}")