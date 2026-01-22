import io
from .checker import StructureChecker
from .forcefield import ForceFieldValidator
from .conformations import ConformationManager


class ValidationService:
    def __init__(self):
        self.ff_validator = ForceFieldValidator()
        self.conf_manager = ConformationManager()

    def validate_pdb_content(self, pdb_content: str, label: str):
        # 1. Kontrola AltLocs
        alt_locs = self.conf_manager.detect_alt_locs(pdb_content)

        # 2. PDBFixer kontroly
        checker = StructureChecker(pdb_content)
        structure_results = checker.run_diagnostics()

        # 3. Force Field kontroly
        unsupported = self.ff_validator.check_residue_compatibility(checker.fixer.topology)

        # Sestavení finálního reportu
        errors = structure_results["errors"]
        if alt_locs:
            errors.append({"issue": "alt_locs_detected", "details": alt_locs})
        if unsupported:
            errors.append({"issue": "unsupported_residues", "details": unsupported})

        return {
            "label": label,
            "is_ready_for_hpc": len(errors) == 0,
            "errors": errors,
            "stats": structure_results["stats"]
        }