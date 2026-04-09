import logging
from pathlib import Path
from typing import Dict, List, Any

from app.services.pdb_service import parse_pdb_to_topology_dict
from app.services.forcefield_service import ForceFieldService
from app.utils.adams4sims_processing_library.utils import AMBER_topology
from app.utils.adams4sims_processing_library import FF_IDA
from app.workspaces.manager import WorkspaceManager

logger = logging.getLogger(__name__)


class TopologyService:
    def __init__(self):
        self.ff_service = ForceFieldService()
        self.workspace_manager = WorkspaceManager()

    def generate_topology(self, workspace_id: str, pdb_filename: str, ff_selections: Dict[str, Any]):
        """
        Hlavní pipeline pro generování topologie.
        ff_selections: např. {"R": ff_data_z_api_pro_OL3}
        """
        try:
            # Načtení PDB obsahu z workspace
            pdb_path = self.workspace_manager.get_file_path(workspace_id, pdb_filename)
            with open(pdb_path, "r") as f:
                pdb_content = f.read()

            logger.info(f"Starting topology generation for workspace {workspace_id}")

            # 2. PDB -> JSON (vytvoření 'mol' struktury)
            # Mapování typů (R, D...) na názvy FF
            ff_mapping = {mol_type: data['name'] for mol_type, data in ff_selections.items()}
            mol = parse_pdb_to_topology_dict(pdb_content, ff_mapping)

            # 3. Příprava Force Field souborů
            mol['force_field_data'] = {}
            for mol_type, ff_data in ff_selections.items():
                ff_name = ff_data['name']
                # prepare_forcefield_files uloží soubory na disk a vrátí cestu
                ff_path = self.ff_service.prepare_forcefield_files(ff_data)

                # Načtení do FF_IDA instance
                mol['force_field_data'][ff_name] = FF_IDA.ff(
                    str(ff_path / f"{ff_name}.rtp"),
                    str(ff_path / f"nonbonded_{ff_name}.itp"),
                    str(ff_path / f"bonded_{ff_name}.itp"),
                    str(ff_path / f"{ff_name}.atp")
                )

            # Samotný výpočet AMBER topologie
            logger.info("Running AMBER topology calculation...")
            topology_data = AMBER_topology.create_AMBER_topology(mol)

            # Uložení výsledku (.prmtop) do workspace
            output_filename = pdb_filename.replace(".pdb", ".prmtop")
            output_path = self.workspace_manager.get_workspace_path(workspace_id) / output_filename

            AMBER_topology.write_AMBER_topology(topology_data, str(output_path))

            logger.info(f"Topology successfully saved to {output_path}")
            return output_filename

        except Exception as e:
            logger.error(f"Topology generation failed: {e}")
            raise