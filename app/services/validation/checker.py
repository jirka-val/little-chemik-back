# app/services/validation/checker.py
import io
from pdbfixer import PDBFixer
from typing import Dict, Any, List
# Importing your existing analysis service to obtain correctly named tokens
from app.services.analysis_service import build_sequence_tokens


class StructureChecker:
    def __init__(self, pdb_content: str):
        self.pdb_content = pdb_content
        # PDBFixer removes alternative locations by default upon loading.
        # To allow the frontend to offer a selection, we must detect them in the raw text.
        self.detected_alt_locs = self._scan_for_alt_locs(pdb_content)

        f = io.StringIO(pdb_content)
        self.fixer = PDBFixer(pdbfile=f)

    def _scan_for_alt_locs(self, pdb_text: str) -> List[Dict]:
        """
        Manual scan of PDB text to detect alternative locations (AltLocs).
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
        Runs complete structure diagnostics for HPC preparation.
        Returns data compatible with the new ValidationService.
        """
        # 1. Sequence analysis (Get tokens with ff_resnames like RU3, DA5 etc.)
        # This function returns a dictionary with the "chains" key
        analysis_data = build_sequence_tokens(self.pdb_content)

        all_tokens = []
        if "chains" in analysis_data:
            for ch_id in analysis_data["chains"]:
                # Extract tokens from all chains into a single list
                all_tokens.extend(analysis_data["chains"][ch_id]["tokens"])

        # 2. Run internal PDBFixer checks
        self.fixer.findMissingResidues()
        self.fixer.findNonstandardResidues()
        self.fixer.findMissingAtoms()

        errors = []
        warnings = []

        # 3. Check for missing heavy atoms
        for residue, atoms in self.fixer.missingAtoms.items():
            atom_names = [atom.name for atom in atoms]
            errors.append({
                "type": "STRUCTURE",
                "resn": residue.name,
                "id": residue.id,
                "chain": residue.chain.id,
                "issue": "missing_atoms",
                "message": f"Missing heavy atoms: {', '.join(atom_names)}",
                "critical": True
            })

        # 4. Check for discontinuity (Missing Residues)
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
                "message": f"Chain {chain_id} is missing a segment of {len(res_names)} residues.",
                "critical": True
            })

        # 5. Water check
        has_water = any(r.name in ['HOH', 'WAT'] for r in self.fixer.topology.residues())
        if has_water:
            warnings.append({
                "type": "CLEANUP",
                "issue": "water_detected",
                "message": "Structure contains crystal water. It is recommended to remove it before simulation."
            })

        # Statistics for frontend
        chain_ids = sorted(list(set(c.id for c in self.fixer.topology.chains())))

        return {
            "errors": errors,
            "warnings": warnings,
            "tokens": all_tokens,  # Key that was missing and caused KeyError
            "stats": {
                "atom_count": sum(1 for _ in self.fixer.topology.atoms()),
                "residue_count": sum(1 for _ in self.fixer.topology.residues()),
                "chain_ids": chain_ids
            }
        }