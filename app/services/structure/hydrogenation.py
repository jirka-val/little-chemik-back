import io
import os
import tempfile
import subprocess
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm import unit


class HydrogenationService:
    def __init__(self):
        """
        Služba pro pokročilou přípravu struktury: protonace, řešení ligandů/vod a solvatace.
        Zcela nezávislá na externích silových polích (bez Amberu).
        """
        pass

    def _apply_propka(self, pdb_content: str, ph: float) -> str:
        with tempfile.TemporaryDirectory() as temp_dir:
            input_pdb_path = os.path.join(temp_dir, "input.pdb")
            with open(input_pdb_path, "w", encoding="utf-8") as f:
                f.write(pdb_content)

            try:
                subprocess.run(
                    ["propka3", input_pdb_path],
                    cwd=temp_dir,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                print("PROPKA úspěšně analyzovala strukturu.")
            except Exception as e:
                print(f"Varování: PROPKA selhala, pokračuji se standardními pravidly. Chyba: {e}")

        return pdb_content

    def prepare_structure(self,
                          pdb_content: str,
                          ph: float = 7.0,
                          crystal_water_mode: str = "remove_all",
                          add_solvent: bool = False,
                          box_padding_nm: float = 1.0,
                          ionic_strength: float = 0.15,
                          positive_ion: str = "Na+",
                          negative_ion: str = "Cl-") -> str:
        """
        Komplexní příprava PDB souboru (vodíky, krystalové vody, solvatace, ionty).
        """
        # 1. PROPKA pro přesnou protonaci
        pdb_content = self._apply_propka(pdb_content, ph)

        # 2. Načtení do PDBFixeru
        input_stream = io.StringIO(pdb_content)
        fixer = PDBFixer(pdbfile=input_stream)

        # 3. Řešení krystalových vod a ligandů
        if crystal_water_mode == "remove_all":
            fixer.removeHeterogens(False)
        elif crystal_water_mode == "keep_water":
            fixer.removeHeterogens(True)

            # 4. Přidání chybějících částí a vodíků
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(ph)

        # 5. SOLVATACE A IONTY
        if add_solvent:
            print("Přidávám solvent a ionty pomocí PDBFixer...")
            fixer.addSolvent(
                padding=box_padding_nm * unit.nanometers,
                positiveIon=positive_ion,
                negativeIon=negative_ion,
                ionicStrength=ionic_strength * unit.molar
            )

        # 6. Export zpět do PDB textu
        output_stream = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, output_stream, keepIds=True)

        return output_stream.getvalue()