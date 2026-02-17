# app/services/validation/service.py

from typing import Dict, Any, List
from .checker import StructureChecker
from .forcefield import ForceFieldValidator
from .conformations import ConformationManager

class ValidationService:
    def __init__(self):
        self.ff_validator = ForceFieldValidator()
        self.conf_manager = ConformationManager()

    def validate_pdb_content(self, pdb_content: str, label: str = "current_state") -> Dict[str, Any]:
        # 1. Detekce AltLocs (Konformace)
        alt_locs = self.conf_manager.detect_alt_locs(pdb_content)

        # 2. Diagnostika struktury (PDBFixer)
        # Checker by měl uvnitř používat PDBFixer
        checker = StructureChecker(pdb_content)
        structure_results = checker.run_diagnostics() # Vrací chybějící atomy, rezidua, atd.

        # Místo: unsupported = self.ff_validator.check_residue_compatibility(checker.fixer.topology)
        # Použijeme výsledek z tvé analýzy, kterou už máš (např. z build_sequence_tokens)
        unsupported = self.ff_validator.check_residue_compatibility(structure_results["tokens"])

        # 4. Sestavení komplexního reportu
        errors = structure_results.get("errors", [])
        warnings = structure_results.get("warnings", [])

        # Zařazení AltLocs do reportu
        if alt_locs:
            errors.append({
                "type": "CONFORMATION",
                "issue": "alt_locs_detected",
                "message": f"Nalezeno {len(alt_locs)} reziduí s více konformacemi.",
                "details": alt_locs,
                "critical": True
            })

        # Zařazení chybějících parametrů do reportu
        if unsupported:
            errors.append({
                "type": "FORCEFIELD",
                "issue": "unsupported_residues",
                "message": "Některá rezidua nejsou definována v converting_dictionary.json.",
                "details": unsupported,
                "critical": True
            })

        # Kontrola vodíků (častý problém pro simulace)
        has_hydrogens = any(atom.element.symbol == 'H' for atom in checker.fixer.topology.atoms())
        if not has_hydrogens:
            warnings.append({
                "type": "STRUCTURE",
                "issue": "missing_hydrogens",
                "message": "V molekule chybí vodíky. Před simulací bude nutná protonace.",
                "critical": False
            })

        return {
            "label": label,
            "summary": {
                "is_ready_for_hpc": len([e for e in errors if e.get("critical")]) == 0,
                "has_warnings": len(warnings) > 0,
                "atom_count": structure_results["stats"]["atom_count"],
                "residue_count": structure_results["stats"]["residue_count"],
                "chain_ids": structure_results["stats"]["chain_ids"]
            },
            "analysis": {
                "errors": errors,
                "warnings": warnings,
                "metadata": {
                    "is_protonated": has_hydrogens,
                    "has_alt_locs": len(alt_locs) > 0
                }
            }
        }

    def apply_alt_loc_selection(self, pdb_content: str, selections: Dict[str, str]) -> str:
        return self.conf_manager.filter_pdb_by_selection(pdb_content, selections)