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

# --- NOVÝ IMPORT NAŠÍ KALKULAČKY PRO VODU ---
from app.utils.water_models import water_extra_points

logger = logging.getLogger(__name__)


class TopologyService:
    def __init__(self):
        self.ff_service = ForceFieldService()
        self.workspace_manager = WorkspaceManager()
        self.pdb_service = PDBService()  # Přidáno pro práci se souřadnicemi vod

    def _inject_eps_into_pdb(self, pdb_content: str, ff_water_instance: Any) -> str:
        """
        Pomocná metoda: Najde vody, spočítá EP a přidá je na konec PDB stringu.
        """
        waters = self.pdb_service.extract_water_coordinates(pdb_content)
        if not waters:
            return pdb_content

        try:
            # Otestujeme na první vodě, jestli FF vůbec EP vyžaduje (např. 4-site/5-site)
            test_eps = water_extra_points(ff_water_instance, waters[0]["crd"])
            if not test_eps:
                logger.info("Water model does not require Extra Points (e.g., 3-site). Skipping injection.")
                return pdb_content
        except Exception as e:
            logger.error(f"Error checking water model EPs: {e}")
            return pdb_content

        logger.info(f"Adding Extra Points for {len(waters)} water molecules...")
        new_atom_lines = []
        atom_offset = 90000  # Bezpečný počáteční index pro EP, aby nekolidoval s existujícími atomy

        for water in waters:
            try:
                eps = water_extra_points(ff_water_instance, water["crd"])
            except ValueError as e:
                logger.warning(f"Skipping water {water['chain']}:{water['resseq']} - {e}")
                continue

            for i, ep_crd in enumerate(eps):
                atom_name = f"EP{i + 1}"
                # Standardní PDB formátování pro HETATM s pevnou šířkou
                line = (
                    f"HETATM{atom_offset:5d} {atom_name:<4} {water['res_name']:>3} {water['chain']}{water['resseq']:4d}    "
                    f"{ep_crd[0]:8.3f}{ep_crd[1]:8.3f}{ep_crd[2]:8.3f}  1.00  0.00          {atom_name[0]:>2}")
                new_atom_lines.append(line)
                atom_offset += 1

        # Vložíme nové atomy těsně před END
        lines = pdb_content.splitlines()
        final_lines = [line for line in lines if line.strip() != "END"]
        final_lines.extend(new_atom_lines)
        final_lines.append("END")

        return "\n".join(final_lines)

    def generate_topology(self, workspace_id: str, pdb_filename: str, ff_selections: Dict[str, Any]) -> Dict[str, str]:
        """
        Main pipeline for generating AMBER topology (.prmtop) and updated PDB from a PDB file.
        Vrací slovník s názvy vygenerovaných souborů.
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
            with open(pdb_path, "r") as f:
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

            # Rozparsování PDB do vnitřní struktury 'mol'
            mol = parse_pdb_to_topology_dict(pdb_content, ff_mapping)

            # 3. Načtení Force Fieldů a příprava instancí
            mol['force_field_data'] = {}
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

                # Registrace k typu molekuly
                mol['force_field_data'][mol_type] = ff_instance
                mol['force_field_data'][ff_name] = ff_instance

                # OPRAVA 1: Propojení DNA ('D') na vnitřní tag Amberu ('R')
                if mol_type == 'D':
                    mol['force_field_data']['R'] = ff_instance
                elif mol_type == 'R':
                    mol['force_field_data']['D'] = ff_instance

                # --- NOVÝ KROK: PŘIDÁNÍ EXTRA POINTS ---
                if mol_type == 'W':
                    pdb_content = self._inject_eps_into_pdb(pdb_content, ff_instance)

            # --- OPRAVA 2: Terminální rezidua (5' a 3') pro DNA/RNA ---
            chains = {}
            for res in mol['residues']:
                c_id = res.get('chain', 'A')
                if c_id not in chains: chains[c_id] = []
                chains[c_id].append(res)

            nucleic_bases = ["DA", "DC", "DG", "DT", "RA", "RC", "RG", "RU", "A", "G", "C", "T", "U"]

            for c_id, r_list in chains.items():
                if not r_list: continue
                first_res = r_list[0]
                last_res = r_list[-1]

                # 5' konec (start)
                if any(first_res['resn'] == base for base in nucleic_bases):
                    orig = first_res['resn']
                    first_res['resn'] = f"{orig}5"
                    logger.info(f"Chain {c_id}: Fixed 5' terminal {orig} -> {first_res['resn']}")

                # 3' konec (end)
                if len(r_list) > 1 and any(last_res['resn'] == base for base in nucleic_bases):
                    orig = last_res['resn']
                    last_res['resn'] = f"{orig}3"
                    logger.info(f"Chain {c_id}: Fixed 3' terminal {orig} -> {last_res['resn']}")
            # ------------------------------------------------------------

            # Aplikace aliasů (HOH -> WAT atd.)
            for res in mol['residues']:
                res['resn'] = resn_alias(res['resn'])

            # 4. Výpočet AMBER topologie
            logger.info("Running AMBER topology calculation...")
            topology_data = AMBER_topology.create_AMBER_topology(mol)

            # 5. Uložení souborů
            # A) Uložení .prmtop (Topologie)
            prmtop_filename = pdb_filename.replace(".pdb", ".prmtop")
            prmtop_path = workspace_dir / prmtop_filename
            AMBER_topology.write_AMBER_topology(str(prmtop_path), topology_data)

            # B) Uložení rozšířeného PDB (Souřadnice)
            # Nejdříve uložíme verzi _ready.pdb (pro explicitní stažení)
            extended_pdb_filename = pdb_filename.replace(".pdb", "_ready.pdb")
            extended_pdb_path = workspace_dir / extended_pdb_filename
            with open(extended_pdb_path, "w", encoding="utf-8") as f:
                f.write(pdb_content)

            # !!! KLÍČOVÝ KROK: Přepíšeme původní soubor (např. structure.pdb) !!!
            # Díky tomu Molstar po updateSilently okamžitě uvidí Extra Pointy
            original_pdb_path = workspace_dir / pdb_filename
            with open(original_pdb_path, "w", encoding="utf-8") as f:
                f.write(pdb_content)

            logger.info(f"Topology successfully saved to {prmtop_path}")
            logger.info(f"Original and Ready PDB updated at {original_pdb_path}")

            # Vracíme slovník (frontend může stále použít _ready.pdb pro stažení)
            return {
                "topology_file": prmtop_filename,
                "coordinates_file": extended_pdb_filename
            }

        except Exception as e:
            logger.error(f"Topology generation failed: {e}")
            logger.error(f"=== TRACEBACK ===\n{traceback.format_exc()}")
            raise e