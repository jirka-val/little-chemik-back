import logging
from typing import List, Dict, Any
from app.services.pdb_service import PDBService
from app.utils.alias import resn_alias, name_alias
from app.utils.NucleicAcidToolbox import get_nucleotides

logger = logging.getLogger(__name__)


# Pomocné třídy pro NucleicAcidToolbox
class Atom:
    def __init__(self, name, resn, resi, chain, coord):
        self.name = name
        self.resn = resn
        self.resi = resi
        self.chain = chain
        self.coord = coord


class PDBModel:
    def __init__(self, atoms: List[Atom]):
        self.atom = atoms


class ValidationService:
    def __init__(self):
        self.pdb_service = PDBService()

    def _parse_pdb_to_model(self, pdb_content: str) -> PDBModel:
        """Převede text PDB na objekt model, který očekává NucleicAcidToolbox."""
        atoms = []
        for line in pdb_content.splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Jednoduchý PDB parser podle pozic v řádku
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21:22].strip()
                resi = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # 1. Aplikace aliasů (standardizace názvů)
                standard_resn = resn_alias(resn)
                standard_name = name_alias(standard_resn, name)

                atoms.append(Atom(standard_name, standard_resn, resi, chain, [x, y, z]))
        return PDBModel(atoms)

    async def validate_molecule(self, pdb_code: str) -> Dict[str, Any]:
        """Hlavní metoda pro validaci molekuly."""
        try:
            # Získání obsahu souboru
            content = await self.pdb_service.fetch_pdb_content(pdb_code.lower())
            model = self._parse_pdb_to_model(content)

            # 2. Kontrola integrity pomocí NucleicAcidToolbox
            # get_nucleotides interně volá check_nucleotide() a check_bond()
            valid_nucleotides = get_nucleotides(model)

            # Analýza výsledků
            total_atoms = len(model.atom)
            num_valid_nucs = len(valid_nucleotides)

            is_ready = num_valid_nucs > 0

            return {
                "pdb_code": pdb_code,
                "is_ready_for_hpc": is_ready,
                "stats": {
                    "total_atoms_parsed": total_atoms,
                    "valid_nucleotides_found": num_valid_nucs
                },
                "details": [
                    {"id": n.id, "type": n.nuc, "chain": n.chain}
                    for n in valid_nucleotides
                ] if is_ready else "Žádné validní nukleové kyseliny nebyly nalezeny."
            }

        except FileNotFoundError:
            return {"error": f"Molekula {pdb_code} nebyla nalezena."}
        except Exception as e:
            logger.exception(f"Chyba při validaci {pdb_code}")
            return {"error": f"Interní chyba validace: {str(e)}"}