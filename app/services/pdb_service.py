# app/services/pdb_service.py
import httpx
import logging
from app.core.config import settings
from pathlib import Path

logger = logging.getLogger(__name__)


class PDBService:
    async def get_remote_pdb_content(self, pdb_code: str) -> str:
        logger.info(f"Stahuji molekulu {pdb_code} z RCSB pro frontend...")

        async with httpx.AsyncClient() as client:
            response = await client.get(f"{settings.RCSB_PDB_URL}/{pdb_code}.pdb")

            if response.status_code != 200:
                raise FileNotFoundError(f"Molekula {pdb_code} neexistuje v RCSB databázi.")

            return response.text

    def get_molecule_types(self, pdb_content: str) -> list[str]:
        """
        Rychlá detekce typů v PDB pro externí API.
        """
        found = set()
        # Mapování reziduí na kódy externího API
        mapping = {
            "D": {"DA", "DC", "DG", "DT"},
            "R": {"A", "C", "G", "U"},
            "P": {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"},
            "W": {"HOH", "WAT", "SOL"}
        }

        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                for api_code, residues in mapping.items():
                    if res_name in residues:
                        found.add(api_code)
        return list(found)

def parse_pdb_to_topology_dict(pdb_content: str, selected_force_fields: dict = None) -> dict:
    if selected_force_fields is None:
        selected_force_fields = {"R": "OL3"}

    box = [0.0, 0.0, 0.0]
    chains = {}
    seen_residues = set()

    for line in pdb_content.splitlines():
        if line.startswith("CRYST1"):
            try:
                box = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
            except ValueError:
                pass

        elif line.startswith("ATOM  ") or line.startswith("HETATM"):
            res_name = line[17:20].strip()
            chain_id = line[21]
            try:
                res_seq = int(line[22:26])
            except ValueError:
                continue

            if res_name in ['HOH', 'WAT', 'SOL']:
                continue

            res_key = (chain_id, res_seq)

            if res_key not in seen_residues:
                seen_residues.add(res_key)
                if chain_id not in chains:
                    chains[chain_id] = []
                chains[chain_id].append(res_name)

    final_residues = []

    for chain_id, res_list in chains.items():
        chain_length = len(res_list)

        for i, res_name in enumerate(res_list):
            mol_type = "R"

            is_nucleotide = len(res_name) == 1 and res_name in ['A', 'C', 'G', 'U'] or res_name.startswith('R')

            base_name = res_name
            if is_nucleotide and len(base_name) == 1:
                base_name = 'R' + base_name

            if is_nucleotide:
                if i == 0 and not base_name.endswith('5'):
                    final_resn = base_name + "5"
                elif i == chain_length - 1 and not base_name.endswith('3'):
                    final_resn = base_name + "3"
                else:
                    final_resn = base_name
            else:
                final_resn = base_name

            final_residues.append({
                "resn": final_resn,
                "mol_type": mol_type
            })

    return {
        "residues": final_residues,
        "force_fields": selected_force_fields,
        "box": box
    }


def remove_residue_from_pdb(pdb_path: Path, chain: str, resseq: int) -> bool:
    """
    Odstraní všechna data (ATOM i HETATM) pro dané reziduum z PDB souboru.
    Vrací True, pokud byla provedena změna.
    """
    if not pdb_path.exists():
        logger.error(f"PDB file not found at {pdb_path}")
        return False

    try:
        with open(pdb_path, "r") as f:
            lines = f.readlines()

        new_lines = []
        found = False

        for line in lines:
            # PDB záznamy atomů mají pevnou šířku, musí mít aspoň 26 znaků
            if len(line) >= 26 and line.startswith(("ATOM", "HETATM")):
                try:
                    line_chain = line[21].strip()
                    line_resseq = int(line[22:26].strip())

                    if line_chain == chain and line_resseq == resseq:
                        found = True
                        continue  # Mažeme tento atom
                except (ValueError, IndexError):
                    # Pokud narazíme na nečitelný ResSeq, řádek raději necháme
                    pass

            new_lines.append(line)

        if found:
            with open(pdb_path, "w") as f:
                f.writelines(new_lines)
            logger.info(f"Successfully removed residue {resseq} from chain {chain} in {pdb_path.name}")
            return True

        logger.warning(f"Residue {resseq} in chain {chain} not found in {pdb_path.name}")
        return False

    except Exception as e:
        logger.error(f"Error while editing PDB file: {e}")
        return False