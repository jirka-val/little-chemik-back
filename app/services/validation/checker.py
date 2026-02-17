# app/services/validation/checker.py
import io
from pdbfixer import PDBFixer
from typing import Dict, Any, List
# Importujeme tvou existující službu pro analýzu, abychom získali správně pojmenované tokeny
from app.services.analysis_service import build_sequence_tokens


class StructureChecker:
    def __init__(self, pdb_content: str):
        self.pdb_content = pdb_content
        # PDBFixer standardně maže alternativní polohy při načtení.
        # Aby frontend mohl nabídnout jejich výběr, musíme je detekovat v surovém textu.
        self.detected_alt_locs = self._scan_for_alt_locs(pdb_content)

        f = io.StringIO(pdb_content)
        self.fixer = PDBFixer(pdbfile=f)

    def _scan_for_alt_locs(self, pdb_text: str) -> List[Dict]:
        """
        Ruční scan PDB textu pro detekci alternativních poloh (AltLocs).
        """
        alt_data = {}
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                alt_id = line[16].strip()
                if alt_id:
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_id = line[22:26].strip()

                    key = (chain_id, res_id, res_name)
                    if key not in alt_data:
                        alt_data[key] = set()
                    alt_data[key].add(alt_id)

        return [
            {
                "chain": k[0],
                "res_id": k[1],
                "res_name": k[2],
                "variants": sorted(list(v))
            }
            for k, v in alt_data.items() if len(v) > 1
        ]

    def run_diagnostics(self) -> Dict[str, Any]:
        """
        Spustí kompletní diagnostiku struktury pro HPC přípravu.
        Vrací data kompatibilní s novým ValidationService.
        """
        # 1. Analýza sekvence (Získáme tokeny s ff_resname jako RU3, DA5 atd.)
        # Tato funkce vrací slovník s klíčem "chains"
        analysis_data = build_sequence_tokens(self.pdb_content)

        all_tokens = []
        if "chains" in analysis_data:
            for ch_id in analysis_data["chains"]:
                # Vytaháme tokeny ze všech řetězců do jednoho seznamu
                all_tokens.extend(analysis_data["chains"][ch_id]["tokens"])

        # 2. Spuštění interních kontrol PDBFixeru
        self.fixer.findMissingResidues()
        self.fixer.findNonstandardResidues()
        self.fixer.findMissingAtoms()

        errors = []
        warnings = []

        # 3. Kontrola chybějících těžkých atomů
        for residue, atoms in self.fixer.missingAtoms.items():
            atom_names = [atom.name for atom in atoms]
            errors.append({
                "type": "STRUCTURE",
                "resn": residue.name,
                "id": residue.id,
                "chain": residue.chain.id,
                "issue": "missing_atoms",
                "message": f"Chybějící těžké atomy: {', '.join(atom_names)}",
                "critical": True
            })

        # 4. Kontrola diskontinuity (Missing Residues)
        for chain_res, res_names in self.fixer.missingResidues.items():
            chain = chain_res[0]
            if isinstance(chain, int):
                chains_list = list(self.fixer.topology.chains())
                chain = chains_list[chain] if chain < len(chains_list) else chain

            chain_id = getattr(chain, 'id', str(chain))
            errors.append({
                "type": "STRUCTURE",
                "chain": chain_id,
                "issue": "missing_segment",
                "message": f"V řetězci {chain_id} chybí úsek o délce {len(res_names)} reziduí.",
                "critical": True
            })

        # 5. Kontrola vody
        has_water = any(r.name in ['HOH', 'WAT'] for r in self.fixer.topology.residues())
        if has_water:
            warnings.append({
                "type": "CLEANUP",
                "issue": "water_detected",
                "message": "Struktura obsahuje krystalovou vodu. Doporučujeme ji před simulací odstranit."
            })

        # Statistiky pro frontend
        chain_ids = sorted(list(set(c.id for c in self.fixer.topology.chains())))

        return {
            "errors": errors,
            "warnings": warnings,
            "tokens": all_tokens,  # Klíč, který chyběl a způsoboval KeyError
            "stats": {
                "atom_count": sum(1 for _ in self.fixer.topology.atoms()),
                "residue_count": sum(1 for _ in self.fixer.topology.residues()),
                "chain_ids": chain_ids
            }
        }