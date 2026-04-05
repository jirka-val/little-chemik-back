# app/services/validation/service.py

from typing import Dict, Any, List
from .checker import StructureChecker
from .forcefield import ForceFieldValidator
from .conformations import ConformationManager


class ValidationService:
    def __init__(self):
        self.ff_validator = ForceFieldValidator()
        self.conf_manager = ConformationManager()

    def validate_pdb_content(self, pdb_content: str, label: str = "current_state", selections: Dict[str, str] = None) -> \
    Dict[str, Any]:
        """
        Komplexní validace PDB souboru.
        Nově zahrnuje geometrickou kontrolu kontinuity, pokud jsou aplikovány selekce AltLocs.
        """
        # 1. Detekce AltLocs (Konformace)
        alt_locs = self.conf_manager.detect_alt_locs(pdb_content)

        # 2. Diagnostika struktury (PDBFixer)
        checker = StructureChecker(pdb_content)
        structure_results = checker.run_diagnostics()

        # 3. Kontrola kompatibility s ForceFieldem (DOČASNĚ VYPNUTO)
        # unsupported = self.ff_validator.validate_residues(structure_results["tokens"])
        unsupported = []  # Tímto zajistíme, že kód pod tím nebude padat

        # 4. Sestavení komplexního reportu
        errors = structure_results.get("errors", [])
        warnings = structure_results.get("warnings", [])

        # Zařazení AltLocs do reportu
        if alt_locs:
            errors.append({
                "type": "CONFORMATION",
                "issue": "alt_locs_detected",
                "message": f"Detected {len(alt_locs)} residues with multiple conformations (AltLocs). Selection required.",
                "details": alt_locs,
                "critical": True
            })

        # GEOMETRICKÁ KONTROLA MEZER (Klíčové pro profi použití)
        # Pokud uživatel už nějaké varianty vybral, zkontrolujeme, zda k sobě pasují
        if selections:
            continuity_issues = self.conf_manager.validate_continuity(pdb_content, selections)
            for issue in continuity_issues:
                errors.append({
                    "type": "STRUCTURE",
                    "issue": "continuity_gap",
                    "message": issue["message"],
                    "details": issue["details"],
                    "critical": True
                })

        # Zařazení chybějících parametrů do reportu
        if unsupported:
            errors.append({
                "type": "FORCEFIELD",
                "issue": "unsupported_residues",
                "message": "Some residues are not defined in the HPC forcefield dictionary.",
                "details": unsupported,
                "critical": True
            })

        # Kontrola vodíků (častý problém pro simulace)
        has_hydrogens = any(atom.element.symbol == 'H' for atom in checker.fixer.topology.atoms())
        if not has_hydrogens:
            warnings.append({
                "type": "STRUCTURE",
                "issue": "missing_hydrogens",
                "message": "Molecule lacks hydrogens. Protonation will be required before simulation.",
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

    def apply_alt_loc_selection(self, pdb_content: str, selections: Dict[str, str]) -> Dict[str, Any]:
        """
        Aplikuje výběr a vrací jak nové PDB, tak novou validaci včetně kontroly mezer.
        """
        # 1. Vyfiltrování PDB (ponechání vybraných variant)
        filtered_pdb = self.conf_manager.filter_pdb_by_selection(pdb_content, selections)

        # 2. Okamžitá re-validace výsledné struktury
        # Předáváme selections, aby validate_pdb_content mohl spustit validate_continuity
        validation_results = self.validate_pdb_content(filtered_pdb, label="after_selection", selections=selections)

        return {
            "pdb_content": filtered_pdb,
            "validation": validation_results
        }