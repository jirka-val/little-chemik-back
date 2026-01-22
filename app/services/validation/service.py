from typing import Dict, Any
from .checker import StructureChecker
from .forcefield import ForceFieldValidator
from .conformations import ConformationManager


class ValidationService:
    def __init__(self):
        self.ff_validator = ForceFieldValidator()
        self.conf_manager = ConformationManager()

    def validate_pdb_content(self, pdb_content: str, label: str = "current_state") -> Dict[str, Any]:
        # 1. Kontrola AltLocs
        alt_locs = self.conf_manager.detect_alt_locs(pdb_content)

        # 2. Fyzická analýza (PDBFixer)
        checker = StructureChecker(pdb_content)
        results = checker.run_diagnostics()

        # 3. Kontrola Force Fieldu
        unsupported = self.ff_validator.check_residue_compatibility(checker.fixer.topology)

        # Kompletace chyb
        errors = results["errors"]
        if alt_locs:
            errors.append({
                "issue": "alt_locs_detected",
                "message": "Detekovány alternativní konformace.",
                "details": alt_locs
            })
        if unsupported:
            errors.append({
                "issue": "unsupported_residues",
                "details": unsupported,
                "action_required": "Rezidua nejsou v converting_dictionary.json"
            })

        return {
            "label": label,
            "is_ready_for_hpc": len(errors) == 0,
            "stats": results["stats"],
            "errors": errors,
            "warnings": []
        }

    def apply_alt_loc_selection(self, pdb_content: str, selections: Dict[str, str]) -> str:
        return self.conf_manager.filter_pdb_by_selection(pdb_content, selections)