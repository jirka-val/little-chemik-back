import json
import logging
from pathlib import Path
from typing import Dict, Any
import traceback

from app.services.pdb_service import PDBService, parse_pdb_to_topology_dict
from app.services.forcefield_service import ForceFieldService
from app.services.topology_patches import apply_topology_patches
from app.utils.adams4sims_processing_library.utils import AMBER_topology
from app.utils.adams4sims_processing_library import FF_IDA
from app.workspaces.manager import WorkspaceManager
from app.utils.adams4sims_processing_library.utils.alias import resn_alias

logger = logging.getLogger(__name__)


class TopologyService:
    def __init__(self):
        self.ff_service = ForceFieldService()
        self.workspace_manager = WorkspaceManager()
        self.pdb_service = PDBService()

    def generate_topology(self, workspace_id: str, pdb_filename: str, ff_selections: Dict[str, Any]) -> Dict[str, str]:
        """
        Main pipeline for generating AMBER topology (.prmtop) and updated PDB from a PDB file.
        """
        try:
            # --- Ladění (Debug): Uložení příchozích dat ---
            workspace_dir = self.workspace_manager.get_workspace_dir(workspace_id)
            debug_json_path = workspace_dir / "received_ff_selections.json"

            with open(debug_json_path, "w", encoding="utf-8") as json_file:
                json.dump(ff_selections, json_file, indent=4, ensure_ascii=False)

            logger.info(f"Saved incoming forcefield selections to {debug_json_path}")
            # -----------------------------------------------

            # 1. Načtení původního PDB obsahu
            pdb_path = self.workspace_manager.get_file_path(workspace_id, pdb_filename)
            with open(pdb_path, "r", encoding="utf-8") as f:
                pdb_content = f.read()

            logger.info(f"Starting topology generation for workspace {workspace_id}")

            # 2. Definice jmen pro parser
            rna_names = ["RU5", "RU3", "RU", "RA5", "RA3", "RA", "RC5", "RC3", "RC",
                         "C5", "C3", "C", "RG5", "RG3", "RG", "G5", "G3", "G"]
            dna_names = ["DA", "DA5", "DA3", "DC", "DC5", "DC3", "DG", "DG5", "DG3", "DT", "DT5", "DT3"]
            water_names = ["HOH", "WAT", "SOL", "W", "W3", "W4", "W5"]
            ion_names = ["NA", "Na+", "CL", "Cl-", "K", "K+", "MG", "Mg2+", "CA", "Ca2+"]

            ff_mapping = {}
            for mol_type, data in ff_selections.items():
                raw_name = data.get('display_name') or data.get('ff_name') or 'unknown_ff'
                ff_name = raw_name.replace(" ", "_")
                ff_mapping[mol_type] = ff_name

                if mol_type == 'R':
                    for name in rna_names: ff_mapping[name] = ff_name
                elif mol_type == 'D':
                    for name in dna_names: ff_mapping[name] = ff_name
                elif mol_type == 'W':
                    for name in water_names: ff_mapping[name] = ff_name
                elif mol_type == 'I':
                    for name in ion_names: ff_mapping[name] = ff_name

            # 3. PŘEDNAČTENÍ SILOVÝCH POLÍ (Potřebujeme je pro výpočet Extra Points)
            preloaded_ff_data = {}
            water_ff_instance = None

            for mol_type, ff_data in ff_selections.items():
                raw_name = ff_data.get('display_name') or ff_data.get('ff_name') or 'unknown_ff'
                ff_name = raw_name.replace(" ", "_")

                ff_path = self.ff_service.prepare_forcefield_files(ff_data)

                ff_instance = FF_IDA.ff(
                    str(ff_path / f"{ff_name}.rtp"),
                    str(ff_path / f"nonbonded_{ff_name}.itp"),
                    str(ff_path / f"bonded_{ff_name}.itp"),
                    str(ff_path / f"{ff_name}.atp")
                )

                # Aplikace záplat (Ionty, Voda)
                apply_topology_patches(ff_instance, mol_type)

                # Uložení instancí pro pozdější použití
                preloaded_ff_data[mol_type] = ff_instance
                preloaded_ff_data[ff_name] = ff_instance

                # OPRAVA 1: Propojení DNA ('D') na vnitřní tag Amberu ('R')
                if mol_type == 'D':
                    preloaded_ff_data['R'] = ff_instance
                elif mol_type == 'R':
                    preloaded_ff_data['D'] = ff_instance

                # Pokud je to voda, uložíme si ji bokem pro parser
                if mol_type == 'W':
                    water_ff_instance = ff_instance

            # 4. OPRAVA PDB (Srovnání iontů a vstříknutí EP)
            logger.info("Reordering PDB and calculating Extra Points...")
            fixed_pdb_content = self.pdb_service.reorder_and_inject_eps(pdb_content, water_ff_instance)

            # 5. ROZPARSOVÁNÍ UŽ OPRAVENÉHO PDB DO TOPOLOGIE
            mol = parse_pdb_to_topology_dict(fixed_pdb_content, ff_mapping)
            mol['force_field_data'] = preloaded_ff_data

            # Aplikace aliasů (HOH -> WAT atd.)
            for res in mol['residues']:
                res['resn'] = resn_alias(res['resn'])

            # 6. Výpočet AMBER topologie z opravených dat
            logger.info("Running AMBER topology calculation...")
            topology_data = AMBER_topology.create_AMBER_topology(mol)

            # 7. Uložení souborů
            # A) Uložení .prmtop (Topologie)
            prmtop_filename = pdb_filename.replace(".pdb", ".prmtop")
            prmtop_path = workspace_dir / prmtop_filename
            AMBER_topology.write_AMBER_topology(str(prmtop_path), topology_data)

            # B) Uložení opraveného PDB (Souřadnice)
            # Přepíšeme pouze ten původní soubor ve workspace,
            # aby frontend (Molstar) hned viděl ty správné atomy a Extra Pointy.
            original_pdb_path = workspace_dir / pdb_filename
            with open(original_pdb_path, "w", encoding="utf-8") as f:
                f.write(fixed_pdb_content)

            logger.info(f"Topology successfully saved to {prmtop_path}")
            logger.info(f"PDB coordinates updated at {original_pdb_path}")

            # Vracíme slovník - coordinates_file teď odkazuje na ten samý soubor (structure.pdb)
            return {
                "topology_file": prmtop_filename,
                "coordinates_file": pdb_filename
            }

        except Exception as e:
            logger.error(f"Topology generation failed: {e}")
            logger.error(f"=== TRACEBACK ===\n{traceback.format_exc()}")
            raise e