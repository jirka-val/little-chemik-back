import io
from pdbfixer import PDBFixer
from openmm.app import PDBFile


class HydrogenationService:
    def __init__(self):
        pass

    def add_hydrogen_atoms(self, pdb_content: str, ph: float = 7.0) -> str:
        """
        Načte PDB obsah, doplní chybějící vodíky pomocí PDBFixeru
        při zadaném pH a vrátí nový PDB řetězec.
        """
        # 1. Načtení PDB obsahu do PDBFixeru
        input_stream = io.StringIO(pdb_content)
        fixer = PDBFixer(pdbfile=input_stream)

        # 2. Identifikace chybějících atomů a reziduí (nutné před přidáním vodíků)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()

        # 3. Přidání vodíků
        # PDBFixer automaticky určí protonační stavy na základě pH
        fixer.addMissingHydrogens(ph)

        # 4. Export výsledku zpět do PDB formátu (string)
        output_stream = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, output_stream, keepIds=True)

        return output_stream.getvalue()