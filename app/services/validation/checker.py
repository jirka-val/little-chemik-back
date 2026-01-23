import io
from pdbfixer import PDBFixer
from typing import Dict, Any, List


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
        Ruční scan PDB textu pro detekci alternativních poloh (AltLocs),
        protože PDBFixer je automaticky odstraňuje.
        """
        alt_data = {}
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                alt_id = line[16].strip()  # Sloupec 17 v PDB (index 16)
                if alt_id:
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_id = line[22:26].strip()

                    key = (chain_id, res_id, res_name)
                    if key not in alt_data:
                        alt_data[key] = set()
                    alt_data[key].add(alt_id)

        # Vrátíme pouze rezidua, která mají více než jednu variantu
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
        Vrací data kompatibilní s frontendem Little Chemik.
        """
        self.fixer.findMissingResidues()
        self.fixer.findNonstandardResidues()
        self.fixer.findMissingAtoms()

        errors = []

        # 1. Detekované AltLocs (Předáme frontendu pro zobrazení panelu výběru)
        if self.detected_alt_locs:
            errors.append({
                "issue": "alt_locs_detected",
                "details": self.detected_alt_locs
            })

        # 2. Kontrola chybějících těžkých atomů
        for residue, atoms in self.fixer.missingAtoms.items():
            atom_names = [atom.name for atom in atoms]
            errors.append({
                "resn": residue.name,
                "id": residue.id,  # OPRAVENO: .id je skutečné PDB číslo, ne .index
                "chain": residue.chain.id,  # PŘIDÁNO: ID řetězce pro přesný zoom v Molstaru
                "issue": f"Chybějící těžké atomy: {', '.join(atom_names)}"
            })

        # 3. Kontrola diskontinuity (Missing Residues)
        for chain_res, res_names in self.fixer.missingResidues.items():
            chain = chain_res[0]

            # Pokud je chain pouze index (int), musíme získat skutečný objekt z topologie
            if isinstance(chain, int):
                chain = list(self.fixer.topology.chains())[chain]

            chain_id = chain.id
            errors.append({
                "chain": chain_id,
                "issue": f"V řetězci {chain_id} chybí úsek o délce {len(res_names)} reziduí."
            })

        # 4. Kontrola vodíků
        has_hydrogens = any(atom.element.symbol == 'H' for atom in self.fixer.topology.atoms())
        if not has_hydrogens:
            errors.append({
                "issue": "Molekula neobsahuje žádné vodíky. Je nutná protonace (Add Hydrogens).",
                "type": "protonation_needed"
            })

        # 5. Kontrola vody (Varování)
        has_water = any(r.name in ['HOH', 'WAT'] for r in self.fixer.topology.residues())
        warnings = []
        if has_water:
            warnings.append({
                "issue": "Struktura obsahuje krystalovou vodu. Pro čistou simulaci je vhodné ji odstranit."
            })

        # Celkový stav připravenosti
        # Považujeme za nepřipravené, pokud jsou tam chyby (kromě AltLocs, ty se vyřeší výběrem)
        real_errors = [e for e in errors if e['issue'] != 'alt_locs_detected']
        is_ready = len(real_errors) == 0 and not self.detected_alt_locs

        return {
            "errors": errors,
            "warnings": warnings,
            "is_ready_for_hpc": is_ready,
            "stats": {
                "atom_count": sum(1 for _ in self.fixer.topology.atoms()),
                "residue_count": sum(1 for _ in self.fixer.topology.residues())
            }
        }