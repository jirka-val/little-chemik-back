import json
import logging
from pathlib import Path
from typing import Dict, Any
import traceback

from app.services.pdb_service import parse_pdb_to_topology_dict
from app.services.forcefield_service import ForceFieldService
from app.utils.adams4sims_processing_library.utils import AMBER_topology
from app.utils.adams4sims_processing_library import FF_IDA
from app.workspaces.manager import WorkspaceManager
from app.utils.adams4sims_processing_library.utils.alias import resn_alias

logger = logging.getLogger(__name__)


class TopologyService:
    def __init__(self):
        self.ff_service = ForceFieldService()
        self.workspace_manager = WorkspaceManager()

    def generate_topology(self, workspace_id: str, pdb_filename: str, ff_selections: Dict[str, Any]) -> str:
        """
        Main pipeline for generating AMBER topology (.prmtop) from a PDB file.
        """
        try:
            # --- PŘIDÁNO: Uložení příchozího JSONu pro ladění (Debug) ---
            workspace_dir = self.workspace_manager.get_workspace_dir(workspace_id)
            debug_json_path = workspace_dir / "received_ff_selections.json"

            with open(debug_json_path, "w", encoding="utf-8") as json_file:
                json.dump(ff_selections, json_file, indent=4, ensure_ascii=False)

            logger.info(f"Saved incoming forcefield selections to {debug_json_path}")
            # ------------------------------------------------------------

            # 1. Načtení PDB obsahu
            pdb_path = self.workspace_manager.get_file_path(workspace_id, pdb_filename)
            with open(pdb_path, "r") as f:
                pdb_content = f.read()

            logger.info(f"Starting topology generation for workspace {workspace_id}")

            # 2. Mapování pro parser
            # Ponecháváme tvoje moderní názvosloví (RU5, RA...), aby parser věděl,
            # že jde o molekuly typu RNA ('R').
            rna_names = ["RU5", "RU3", "RU", "RA5", "RA3", "RA", "RC5", "RC3", "RC",
                         "C5", "C3", "C", "RG5", "RG3", "RG", "G5", "G3", "G"]

            ff_mapping = {}
            for mol_type, data in ff_selections.items():
                raw_name = data.get('display_name') or data.get('ff_name') or 'unknown_ff'
                ff_name = raw_name.replace(" ", "_")
                ff_mapping[mol_type] = ff_name
                if mol_type == 'R':
                    for name in rna_names:
                        ff_mapping[name] = ff_name

            # Rozparsování PDB do vnitřní struktury 'mol'
            mol = parse_pdb_to_topology_dict(pdb_content, ff_mapping)

            # 3. Načtení Force Fieldů
            mol['force_field_data'] = {}
            for mol_type, ff_data in ff_selections.items():
                raw_name = ff_data.get('display_name') or ff_data.get('ff_name') or 'unknown_ff'
                ff_name = raw_name.replace(" ", "_")

                # ForceFieldService teď zploští RTP a přidá residue_lib (RU5, RG...)
                ff_path = self.ff_service.prepare_forcefield_files(ff_data)

                # Vytvoření instance silového pole pomocí šéfovy knihovny
                ff_instance = FF_IDA.ff(
                    str(ff_path / f"{ff_name}.rtp"),
                    str(ff_path / f"nonbonded_{ff_name}.itp"),
                    str(ff_path / f"bonded_{ff_name}.itp"),
                    str(ff_path / f"{ff_name}.atp")
                )

                # --- ZÁPLATA PRO RIGIDNÍ VODU (TIP3P, TIP5P, SPCE) ---
                if mol_type == 'W':
                    if hasattr(ff_instance, 'b') and isinstance(ff_instance.b, dict):
                        if 'angletypes' not in ff_instance.b:
                            ff_instance.b['angletypes'] = {}

                        # 1. Varianta pro modely rodiny TIP (TIP3P, TIP4P, TIP5P)
                        ff_instance.b['angletypes'][('WH', 'WH', 'WO')] = [0.0, 37.74, 0.0]
                        ff_instance.b['angletypes'][('WO', 'WH', 'WH')] = [0.0, 37.74, 0.0]
                        ff_instance.b['angletypes'][('WH', 'WO', 'WH')] = [0.0, 104.52, 0.0]

                        # Virtuální body pro TIP4P/TIP5P
                        ff_instance.b['angletypes'][('WEP', 'WO', 'WEP')] = [0.0, 109.47, 0.0]
                        ff_instance.b['angletypes'][('WEP', 'WO', 'WH')] = [0.0, 109.47, 0.0]
                        ff_instance.b['angletypes'][('WH', 'WO', 'WEP')] = [0.0, 109.47, 0.0]

                        # 2. Varianta pro modely rodiny SPC (SPC, SPCE)
                        ff_instance.b['angletypes'][('HW', 'HW', 'OW')] = [0.0, 37.74, 0.0]
                        ff_instance.b['angletypes'][('OW', 'HW', 'HW')] = [0.0, 37.74, 0.0]
                        ff_instance.b['angletypes'][('HW', 'OW', 'HW')] = [0.0, 109.47, 0.0]
                # -----------------------------------------------------


                # Registrace instance k typu molekuly (např. 'R' nebo 'W')
                mol['force_field_data'][mol_type] = ff_instance
                mol['force_field_data'][ff_name] = ff_instance

            # Aplikace aliasů před generováním (HOH -> WAT atd.)
            for res in mol['residues']:
                res['resn'] = resn_alias(res['resn'])

            # 4. Samotný výpočet AMBER topologie
            logger.info("Running AMBER topology calculation...")
            topology_data = AMBER_topology.create_AMBER_topology(mol)

            # 5. Uložení výsledného .prmtop souboru
            output_filename = pdb_filename.replace(".pdb", ".prmtop")
            output_path = self.workspace_manager.get_workspace_dir(workspace_id) / output_filename

            AMBER_topology.write_AMBER_topology(str(output_path), topology_data)

            logger.info(f"Topology successfully saved to {output_path}")
            return output_filename

        except Exception as e:
            # Detailní logování chyby
            logger.error(f"Topology generation failed: {e}")
            logger.error(f"=== TRACEBACK ===\n{traceback.format_exc()}")
            raise