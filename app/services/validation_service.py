import io
import logging
from typing import Dict, Any
from pdbfixer import PDBFixer

logger = logging.getLogger(__name__)


class ValidationService:
    def validate_pdb_content(self, pdb_content: str, label: str = "current_state") -> Dict[str, Any]:

        try:
            f = io.StringIO(pdb_content)
            fixer = PDBFixer(pdbfile=f)

            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.findMissingAtoms()

            errors = []  # Blokační chyby
            warnings = []  # Informační varování

            # ---  KONTROLA TĚŽKÝCH ATOMŮ (Kritické) ---
            for residue, atoms in fixer.missingAtoms.items():
                atom_names = [atom.name for atom in atoms]
                chain_id = getattr(residue.chain, 'id', str(residue.chain))
                errors.append({
                    "resn": residue.name,
                    "id": residue.index,
                    "chain": chain_id,
                    "issue": f"Chybějící těžké atomy: {', '.join(atom_names)}"
                })

            # ---  KONTROLA PŘERUŠENÍ ŘETĚZCE (Kritické) ---
            for chain_res, res_names in fixer.missingResidues.items():
                chain_item = chain_res[0]

                # Zjištění ID řetězce bezpečně
                if hasattr(chain_item, 'id'):
                    chain_id = chain_item.id
                elif isinstance(chain_item, int):
                    try:
                        chains = list(fixer.topology.chains())
                        chain_id = chains[chain_item].id
                    except (IndexError, AttributeError):
                        chain_id = f"Index:{chain_item}"
                else:
                    chain_id = str(chain_item)

                errors.append({
                    "chain": chain_id,
                    "issue": f"V řetězci chybí úsek o délce {len(res_names)} reziduí (diskontinuita)."
                })

            # ---  KONTROLA VODÍKŮ (Kritické) ---
            has_hydrogens = any(atom.element.symbol == 'H' for atom in fixer.topology.atoms())
            if not has_hydrogens:
                errors.append({
                    "issue": "Molekula neobsahuje žádné vodíky. Je nutná protonace (Add Hydrogens)."
                })

            # ---  NESTANDARDNÍ REZIDUA (Kritické/Varování) ---
            for res in fixer.nonstandardResidues:
                errors.append({
                    "resn": res.name,
                    "id": res.index,
                    "issue": f"Nestandardní reziduum ({res.name}). Chybí parametry Force Field."
                })

            # ---  HETEROGENY  ---
            water_residues = [r for r in fixer.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3']]
            if water_residues:
                warnings.append({
                    "issue": f"Nalezeno {len(water_residues)} molekul krystalové vody (doporučeno odstranit)."
                })

            ions = [r.name for r in fixer.topology.residues() if r.name in ['NA', 'CL', 'MG', 'K', 'ZN']]
            if ions:
                warnings.append({
                    "issue": f"Detekovány ionty: {', '.join(set(ions))}."
                })

            is_ready = len(errors) == 0

            return {
                "label": label,
                "is_ready_for_hpc": is_ready,
                "stats": {
                    "total_errors": len(errors),
                    "total_warnings": len(warnings),
                    "atom_count": sum(1 for _ in fixer.topology.atoms())
                },
                "errors": errors,
                "warnings": warnings
            }

        except Exception as e:
            # Přidáno exc_info=True, abys v konzoli viděl přesné místo chyby
            logger.error(f"Chyba při komplexní diagnostice molekuly: {str(e)}", exc_info=True)
            return {"error": f"Diagnostika selhala: {str(e)}"}