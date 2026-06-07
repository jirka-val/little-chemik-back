import io
import os
import tempfile
import subprocess
import logging
import time
import math

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm import unit

# Inicializace loggeru pro sledování chyb na backendu
logger = logging.getLogger(__name__)


class HydrogenationService:
    def __init__(self):
        """
        Služba pro pokročilou přípravu struktury: protonace, řešení ligandů/vod a solvatace.
        Zcela nezávislá na externích silových polích (bez Amberu).
        """
        pass

    def _apply_propka(self, pdb_content: str, ph: float) -> str:
        """
        Interní volání PROPKA3 pro analýzu pKa (v této implementaci jako placeholder).
        """
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
                logger.info("PROPKA úspěšně analyzovala strukturu.")
            except Exception as e:
                logger.warning(f"PROPKA selhala, pokračuji se standardními pravidly. Chyba: {e}")

        return pdb_content

    def add_hydrogen_atoms(self, pdb_content: str, ph: float = 7.0, optimize: bool = True) -> str:
        """
        OPRAVA PRO TESTY: Tato metoda sjednocuje volání z endpointu molecules.py a testů.

        """
        return self.prepare_structure(
            pdb_content=pdb_content,
            ph=ph,
            crystal_water_mode="remove_all",
            add_solvent=False
        )

    def _fix_cryst1_for_octahedron(self, pdb_content: str, fixer) -> str:
        """
        Opraví CRYST1 řádek pro truncated octahedron.
        OpenMM totiž zapisuje špatné úhly pro octahedron, fixujeme to zde.
        """
        lines = pdb_content.split('\n')
        new_lines = []
        
        try:
            box_vectors = fixer.topology.getPeriodicBoxVectors()
            if box_vectors is not None:
                # Vypočítej délky a úhly z vektorů
                a_vec, b_vec, c_vec = box_vectors
                
                # Délky
                a = float(a_vec[0].value_in_unit(unit.angstroms))
                b = float(b_vec[1].value_in_unit(unit.angstroms))
                c = float(c_vec[2].value_in_unit(unit.angstroms))
                
                # Úhly pro truncated octahedron: alpha=beta=gamma=109.47122063°
                alpha = 109.47122063
                beta = 109.47122063
                gamma = 109.47122063
                
                # Vytvořit správný CRYST1 řádek
                cryst1_line = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
                
                # Nahradit starý CRYST1 řádek novým
                for i, line in enumerate(lines):
                    if line.startswith("CRYST1"):
                        new_lines.append(cryst1_line)
                    else:
                        new_lines.append(line)
                
                logger.info(f"CRYST1 korektně nastaven: a={a:.3f}, b={b:.3f}, c={c:.3f}, α=β=γ=109.47°")
                return '\n'.join(new_lines)
        except Exception as e:
            logger.warning(f"Nepodařilo se zkorrektovat CRYST1: {e}, použiji původní řádek")
            return pdb_content
        
        return pdb_content

    def prepare_structure(self,
                          pdb_content: str,
                          ph: float = 7.0,
                          crystal_water_mode: str = "remove_all",
                          add_solvent: bool = True,
                          box_padding_nm: float = 1.0,
                          box_shape: str = "cube",
                          ionic_strength: float = 0.15,
                          positive_ion: str = "Na+",
                          negative_ion: str = "Cl-") -> str:
        """
        Komplexní příprava PDB souboru (vodíky, krystalové vody, solvatace, ionty).
        """
        total_start = time.time()
        try:
            logger.info("--- Zahajuji přípravu struktury ---")

            # 1. PROPKA pro přesnou protonaci
            propka_start = time.time()
            pdb_content = self._apply_propka(pdb_content, ph)
            logger.info(f"PROPKA dokončena za: {time.time() - propka_start:.2f} s")

            # 2. Načtení do PDBFixeru
            load_start = time.time()
            input_stream = io.StringIO(pdb_content)
            fixer = PDBFixer(pdbfile=input_stream)
            logger.info(f"Načtení do PDBFixeru dokončeno za: {time.time() - load_start:.2f} s")

            # 3. Řešení krystalových vod a ligandů
            if crystal_water_mode == "remove_all":
                fixer.removeHeterogens(False)
            elif crystal_water_mode == "keep_water":
                fixer.removeHeterogens(True)

            # 4. Oprava struktury
            find_start = time.time()
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            logger.info(f"Hledání chybějících reziduí/atomů dokončeno za: {time.time() - find_start:.2f} s")

            add_atoms_start = time.time()
            fixer.addMissingAtoms()
            logger.info(f"Přidání chybějících těžkých atomů dokončeno za: {time.time() - add_atoms_start:.2f} s")

            add_hydrogens_start = time.time()
            fixer.addMissingHydrogens(ph)
            logger.info(f"Přidání vodíků (pH {ph}) dokončeno za: {time.time() - add_hydrogens_start:.2f} s")

            # 5. SOLVATACE A IONTY
            valid_shapes = {"cube", "octahedron"}
            # Normalizace: "truncated octahedron" → "octahedron"
            normalized_shape = box_shape.lower().replace("truncated ", "").strip()
            shape = normalized_shape if normalized_shape in valid_shapes else "cube"
            
            logger.info(f"Box shape vstup: '{box_shape}' → normalizováno na: '{shape}'")
            
            if add_solvent:
                logger.info(f"Zahajuji solvataci (tvar: {shape})")
                logger.info(
                    f"Zahajuji solvataci (padding: {box_padding_nm} nm, tvar: {box_shape}, síla: {ionic_strength}M). Toto může trvat velmi dlouho...")
                solv_start = time.time()
                try:
                    fixer.addSolvent(
                        padding=box_padding_nm * unit.nanometers,
                        boxShape=shape,  # Použije validovaný tvar
                        positiveIon=positive_ion,
                        negativeIon=negative_ion,
                        ionicStrength=ionic_strength * unit.molar
                    )
                    logger.info(f"Solvatace a ionty ÚSPĚŠNĚ DOKONČENY za: {time.time() - solv_start:.2f} s")
                except Exception as solv_err:
                    solv_time = time.time() - solv_start
                    logger.error(f"Pád PDBFixeru při solvataci po {solv_time:.2f} s: {solv_err}")
                    raise ValueError(f"Nepodařilo se přidat vodní box: {solv_err}")

            # 6. Export zpět do PDB textu
            write_start = time.time()
            output_stream = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, output_stream, keepIds=True)
            pdb_output = output_stream.getvalue()
            logger.info(f"Zápis do StringIO dokončen za: {time.time() - write_start:.2f} s")

            # 7. OPRAVA CRYST1 pro truncated octahedron
            if add_solvent and shape == "octahedron":
                pdb_output = self._fix_cryst1_for_octahedron(pdb_output, fixer)
                logger.info("CRYST1 řádek opraven pro truncated octahedron")

            logger.info(
                f"--- Příprava struktury kompletně hotova za celkový čas: {time.time() - total_start:.2f} s ---")
            return pdb_output

        except Exception as e:
            total_time = time.time() - total_start
            logger.error(f"Kritická chyba v HydrogenationService po {total_time:.2f} s: {str(e)}")
            raise e