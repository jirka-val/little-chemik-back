# app/services/structure/editor_service.py
import io
import logging
from pdbfixer import PDBFixer
from openmm.app import PDBFile

logger = logging.getLogger(__name__)


class StructureEditorService:
    def __init__(self):
        """
        Služba pro přímou editaci struktury na úrovni reziduí a atomů.
        """
        pass

    def rename_residue(self, pdb_content: str, chain_id: str, res_num: int, new_res_name: str) -> str:
        """
        Přepíše název rezidua (např. HIS -> HIE) beze změny souřadnic.
        Pracuje přímo se sloupci PDB formátu.
        """
        lines = pdb_content.split('\n')
        new_lines = []

        # PDB formát vyžaduje přesně 3 znaky pro název rezidua (sloupce 18-20)
        new_res_padded = f"{new_res_name:<3}"[:3]

        changes_made = 0

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
                try:
                    # Indexy v Pythonu jsou o 1 menší než čísla sloupců v PDB standardu
                    line_chain = line[21].strip()
                    line_res_num = int(line[22:26].strip())

                    if line_chain == chain_id and line_res_num == res_num:
                        # Rozřízneme řádek, vložíme nový název a zase slepíme
                        line = line[:17] + new_res_padded + line[20:]
                        changes_made += 1
                except ValueError:
                    pass  # Ignorujeme řádky, kde chybí čísla (špatný formát)

            new_lines.append(line)

        logger.info(
            f"Přejmenováno reziduum {res_num} v řetězci {chain_id} na {new_res_padded}. Změněno atomů: {changes_made}")
        return '\n'.join(new_lines)

    def rename_atom(self, pdb_content: str, chain_id: str, res_num: int, old_atom_name: str, new_atom_name: str) -> str:
        """
        Přepíše jméno konkrétního atomu v daném reziduu (např. CD1 -> CG2).
        """
        lines = pdb_content.split('\n')
        new_lines = []

        # PDB formát vyžaduje 4 znaky pro atom (sloupce 13-16).
        # Standardně: pokud má atom < 4 znaky, začíná mezerou (např. " CA ")
        if len(new_atom_name) < 4:
            new_atom_padded = f" {new_atom_name:<3}"
        else:
            new_atom_padded = f"{new_atom_name:<4}"[:4]

        changes_made = 0

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
                try:
                    line_chain = line[21].strip()
                    line_res_num = int(line[22:26].strip())
                    line_atom = line[12:16].strip()

                    if line_chain == chain_id and line_res_num == res_num and line_atom == old_atom_name:
                        line = line[:12] + new_atom_padded + line[16:]
                        changes_made += 1
                except ValueError:
                    pass

            new_lines.append(line)

        logger.info(f"Přejmenován atom '{old_atom_name}' na '{new_atom_padded}' v reziduu {res_num} ({chain_id}).")
        return '\n'.join(new_lines)

    def remove_atom(self, pdb_content: str, chain_id: str, res_num: int, atom_name: str) -> str:
        """
        Odstraní konkrétní atom z daného rezidua.
        Tato funkce jednoduše vynechá příslušný řádek ATOM/HETATM v PDB souboru.
        """
        lines = pdb_content.split('\n')
        new_lines = []

        # Očistíme název od mezer pro spolehlivé porovnání (např. " OXT" -> "OXT")
        target_atom = atom_name.strip()
        removed_count = 0

        for line in lines:
            # Zajímá nás jen fyzická struktura, hlavičky ignorujeme
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
                try:
                    # Vytažení hodnot ze správných PDB sloupců
                    line_chain = line[21].strip()
                    line_res_num = int(line[22:26].strip())
                    line_atom = line[12:16].strip()

                    # Pokud je to náš hledaný atom v daném reziduu, přeskočíme ho
                    if line_chain == chain_id and line_res_num == res_num and line_atom == target_atom:
                        removed_count += 1
                        continue  # "continue" způsobí, že se řádek nepřidá do new_lines

                except ValueError:
                    pass  # Špatně zformátovaný řádek ignorujeme

            # Všechny ostatní řádky normálně opíšeme
            new_lines.append(line)

        logger.info(f"Odstraněno {removed_count} atomů s názvem '{target_atom}' v reziduu {res_num} ({chain_id}).")

        if removed_count == 0:
            logger.warning(f"Atom '{target_atom}' nebyl v reziduu {res_num} nalezen. Nic nebylo smazáno.")

        return '\n'.join(new_lines)

    def mutate_residue(self, pdb_content: str, chain_id: str, res_num: int, mutate_to: str) -> str:
        """
        Provede bodovou mutaci rezidua pomocí PDBFixeru.
        Odstraní starý postranní řetězec a dogeneruje 3D souřadnice nového.
        """
        logger.info(f"Zahajuji mutaci: Řetězec {chain_id}, Reziduum {res_num} -> {mutate_to}")

        input_stream = io.StringIO(pdb_content)
        fixer = PDBFixer(pdbfile=input_stream)

        # PDBFixer syntaxe pro mutaci: "NovýNázev-ČísloRezidua-IDŘetězce" (např. "TYR-42-A")
        mutation_query = f"{mutate_to}-{res_num}-{chain_id}"

        try:
            # 1. Aplikace mutace (PDBFixer v paměti "odřízne" staré atomy)
            fixer.applyMutations([mutation_query])

            # 2. Dopočítání chybějících (nových) atomů
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            # 3. Zápis do stringu
            output_stream = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, output_stream, keepIds=True)

            logger.info("Mutace úspěšně dokončena.")
            return output_stream.getvalue()

        except Exception as e:
            logger.error(f"Chyba při mutaci struktury: {str(e)}")
            raise ValueError(f"Nepodařilo se mutovat reziduum: {str(e)}")