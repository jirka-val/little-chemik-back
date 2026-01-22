import io
from pdbfixer import PDBFixer
from typing import Dict, Any, List


class StructureChecker:
    def __init__(self, pdb_content: str):
        self.pdb_content = pdb_content
        f = io.StringIO(pdb_content)
        self.fixer = PDBFixer(pdbfile=f)

    def run_diagnostics(self) -> Dict[str, Any]:
        """Spustí základní analýzu PDBFixer a vrátí nalezené chyby."""
        self.fixer.findMissingResidues()
        self.fixer.findNonstandardResidues()
        self.fixer.findMissingAtoms()

        errors = []

        # Kontrola těžkých atomů
        for residue, atoms in self.fixer.missingAtoms.items():
            atom_names = [atom.name for atom in atoms]
            errors.append({
                "resn": residue.name,
                "id": residue.index,
                "issue": f"Chybějící těžké atomy: {', '.join(atom_names)}"
            })

        # Kontrola diskontinuity (Missing Residues)
        for chain_res, res_names in self.fixer.missingResidues.items():
            chain_id = chain_res[0].id if hasattr(chain_res[0], 'id') else str(chain_res[0])
            errors.append({
                "chain": chain_id,
                "issue": f"V řetězci chybí úsek o délce {len(res_names)} reziduí."
            })

        # Kontrola vodíků
        has_hydrogens = any(atom.element.symbol == 'H' for atom in self.fixer.topology.atoms())
        if not has_hydrogens:
            errors.append({"issue": "Molekula neobsahuje žádné vodíky. Je nutná protonace."})

        return {
            "errors": errors,
            "stats": {
                "atom_count": sum(1 for _ in self.fixer.topology.atoms())
            }
        }