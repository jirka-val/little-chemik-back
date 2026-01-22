import json
import logging
from typing import List, Set

logger = logging.getLogger(__name__)


class ForceFieldValidator:
    def __init__(self, dict_path: str = "converting_dictionary.json"):
        try:
            with open(dict_path, 'r') as f:
                self.ff_dict = json.load(f)
            self.supported_residues = set()
            for category in self.ff_dict.values():
                self.supported_residues.update(category.keys())
        except Exception as e:
            logger.error(f"Nepodařilo se načíst konverzní slovník: {e}")
            self.supported_residues = set()

    def check_residue_compatibility(self, topology) -> List[str]:
        """Kontroluje rezidua včetně terminálních variant (DA5, NALA atd.)."""
        unsupported = set()
        for chain in topology.chains():
            residues = list(chain.residues())
            for i, residue in enumerate(residues):
                res_name = residue.name
                possible_names = [res_name]

                # Inteligentní mapování terminálů
                if len(residues) > 1:
                    if i == 0:  # Začátek řetězce
                        possible_names.extend([f"{res_name}5", f"N{res_name}"])
                    elif i == len(residues) - 1:  # Konec řetězce
                        possible_names.extend([f"{res_name}3", f"C{res_name}"])
                else:
                    possible_names.append(f"{res_name}N")

                if not any(name in self.supported_residues for name in possible_names):
                    unsupported.add(res_name)
        return sorted(list(unsupported))