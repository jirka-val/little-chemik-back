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

    # Lokální cache pro uložení extrahovaných souborů (formát konzumovaný FF_IDA/TopologyService)
    CACHE_DIR = Path("data/ff_cache")

    # Druhý pohled na stejná data, v adresářové/souborové konvenci, kterou čeká
    # app/builder (SolvationVdwParameters.from_force_field_root): adresáře
    # pojmenované "{ff_name}_{mol_type}" a soubory obsahující "residue_lib"/
    # "forcefield" v názvu. Obsahově jde o stejné UMFFF soubory, jen zdvojené
    # pod jinými jmény, aby builder mohl žít vedle TopologyService beze změny.
    FORGE_CACHE_DIR = Path("data/ff_cache_forge")

    def __init__(self):
        try:
            self.CACHE_DIR.mkdir(parents=True, exist_ok=True)
            self.FORGE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
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

    def _ff_name(self, ff_data: Dict[str, Any]) -> str:
        raw_name = ff_data.get('display_name') or ff_data.get('ff_name') or 'unknown_ff'
        return raw_name.replace(" ", "_")

    def prepare_forcefield_files(self, ff_data: Dict[str, Any]) -> Path:
        """
        Extrahuje, zploští a uloží soubory silového pole na disk.
        Spojuje hlavní RTP soubor s knihovnou reziduí a nahrazuje #include ostrými daty.
        """
        ff_name = self._ff_name(ff_data)
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

    def prepare_forge_force_field_directory(self, ff_data: Dict[str, Any], mol_type: str) -> Path:
        """
        Zajistí, že pro dané silové pole existuje i druhá kopie souborů v layoutu,
        který čeká app/builder (SolvationVdwParameters.from_force_field_root):
        adresář "{ff_name}_{mol_type}" obsahující soubory s "residue_lib"/
        "nonbonded"/"bonded"/"forcefield" v názvu.

        Obsahově je to stejný UMFFF formát, jaký už produkuje prepare_forcefield_files()
        (potvrzeno na existujících cache souborech) - jde jen o duplicitní zápis pod
        jinými jmény, aby TopologyService/FF_IDA a builder mohly číst tu samou
        vyexportovanou FF data nezávisle na sobě.
        """
        base_dir = self.prepare_forcefield_files(ff_data)
        ff_name = self._ff_name(ff_data)

        target_dir = self.FORGE_CACHE_DIR / f"{ff_name}_{mol_type}"
        target_dir.mkdir(parents=True, exist_ok=True)

        rtp_path = base_dir / f"{ff_name}.rtp"
        nb_path = base_dir / f"nonbonded_{ff_name}.itp"
        b_path = base_dir / f"bonded_{ff_name}.itp"

        # Celý flattened .rtp obsahuje jak [ defaults ] (fudgeLJ/fudgeQQ), tak
        # sekce [ resname ] s [ atoms ]/[ bonds ] - vyhovuje tedy zároveň roli
        # "residue_lib" i "forcefield" souboru, jen pod dvěma různými jmény.
        shutil.copyfile(rtp_path, target_dir / f"residue_lib_{ff_name}.rtp")
        shutil.copyfile(rtp_path, target_dir / f"forcefield_{ff_name}.itp")
        shutil.copyfile(nb_path, target_dir / f"nonbonded_{ff_name}.itp")
        shutil.copyfile(b_path, target_dir / f"bonded_{ff_name}.itp")

        logger.info(f"Forge-compatible force field view ready: {target_dir}")
        return target_dir

    def fetch_all_forcefields(self) -> List[Dict[str, Any]]:
        """
        Stáhne KOMPLETNÍ, nefiltrovaný seznam FF z IDA. Jediné místo, které smí
        volat EXTERNAL_URL - používá ho jak get_matching_forcefields (filtruje
        výsledek pro jeden request), tak FFCatalogService.refresh_catalog()
        (ukládá celý výsledek na disk, aby GET /api/forcefields nemusel při
        každém požadavku čekat na IDA - viz ff_catalog_service.py).

        Vyhazuje výjimku při selhání (na rozdíl od get_matching_forcefields,
        která ji dřív polykala a vracela [] - to bylo v pořádku, dokud FF
        seznam nebyl nikde perzistovaný, ale FFCatalogService potřebuje umět
        selhání refreshe rozlišit od "IDA momentálně nemá žádné FF").
        """
        headers = {"x-client-version": "0.1.0"}
        logger.info(f"Fetching force fields from {self.EXTERNAL_URL}")
        response = requests.get(self.EXTERNAL_URL, headers=headers, timeout=30)
        response.raise_for_status()
        return response.json()

    @staticmethod
    def filter_forcefields(all_ffs: List[Dict[str, Any]], molecule_types: List[str]) -> List[Dict[str, Any]]:
        """
        Čistá filtrovací funkce bez síťového volání - použitelná jak nad
        čerstvě staženým seznamem, tak nad lokálním katalogovým snapshotem
        (viz FFCatalogService.get_forcefields v ff_catalog_service.py).
        """
        search_types = set(molecule_types)

        if "W" in search_types:
            search_types.update(["W3", "W4", "W5"])

        if "I" in search_types:
            search_types.update(["I1", "I1+", "Im", "Im+"])

        matched_ffs = []
        for ff in all_ffs:
            ff_types = ff.get("molecule_type") or []
            if any(t in search_types for t in ff_types):
                matched_ffs.append(ff)
        return matched_ffs

    def get_matching_forcefields(self, molecule_types: List[str]) -> List[Dict[str, Any]]:
        """Stáhne seznam FF živě z API a vyfiltruje ty odpovídající molekule.

        POZOR: toto dělá živé síťové volání na IDA při každém volání - GET
        /api/forcefields/{workspace_id} ho už NEPOUŽÍVÁ (viz forcefields.py),
        protože přesně tohle bylo zdrojem pomalosti FF panelu. Zůstává tu pro
        případ, kdy je potřeba čerstvá jednorázová odpověď mimo cache vrstvu.
        """
        try:
            all_ffs = self.fetch_all_forcefields()
            return self.filter_forcefields(all_ffs, molecule_types)
        except Exception as e:
            logger.error(f"External API communication error: {e}")
            return []

    def ff_name(self, ff_data: Dict[str, Any]) -> str:
        """Veřejná verze _ff_name - FFCatalogService/FFClassificationService
        potřebují stejné odvození jména, aby se dalo párovat s force_fields.json."""
        return self._ff_name(ff_data)

    def clear_cache(self):
        """Vymaže celou složku ff_cache i její builder-kompatibilní zrcadlo."""
        for cache_dir in (self.CACHE_DIR, self.FORGE_CACHE_DIR):
            if cache_dir.exists():
                try:
                    shutil.rmtree(cache_dir)
                    cache_dir.mkdir()
                    logger.info(f"Force field cache cleared successfully: {cache_dir}")
                except Exception as e:
                    logger.error(f"Failed to clear force field cache {cache_dir}: {e}")