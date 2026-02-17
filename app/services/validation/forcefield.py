import json
import logging
from typing import List, Set

logger = logging.getLogger(__name__)


class ForceFieldValidator:
    def __init__(self, dict_path: str = "converting_dictionary.json"):
        # Načtení slovníku zůstává stejné, abychom měli proti čemu porovnávat
        try:
            with open(dict_path, 'r', encoding='utf-8') as f:
                self.ff_dict = json.load(f)

            # Vytvoříme si plochý seznam všech klíčů ze všech kategorií (R, D, P, I1...)
            self.supported_names = set()
            for category_data in self.ff_dict.values():
                self.supported_names.update(category_data.keys())

        except Exception as e:
            logger.error(f"Nepodařilo se načíst konverzní slovník: {e}")
            self.supported_names = set()

    def is_supported(self, name: str) -> bool:
        """Jednoduchá otázka: Existuje toto jméno v mém JSONu?"""
        return name in self.supported_names

    def check_residue_compatibility(self, processed_tokens: List[dict]) -> List[str]:
        """
        Místo topologie kontrolujeme už hotové tokeny z analýzy.
        Tím využijeme fakt, že analysis_service už vyřešil RU5, GOL atd.
        """
        unsupported = set()
        for token in processed_tokens:
            if token.get("is_gap"):
                continue

            # Použijeme už vypočítané ff_resname (např. RU3, RA5, GOL)
            name_to_check = token.get("ff_resname") or token.get("pdb_resname")

            if not self.is_supported(name_to_check):
                unsupported.add(name_to_check)

        return sorted(list(unsupported))