# app/services/pdb_service.py
import httpx
import logging
from app.core.config import settings
from pathlib import Path

# --- NOVÝ IMPORT ---
from app.utils.water_models import water_extra_points

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
            "W": {"HOH", "WAT", "SOL"},
            "I": {"NA", "Na+", "CL", "Cl-", "K", "K+", "MG", "Mg2+", "CA", "Ca2+", "LI", "Li+", "RB", "Rb+", "CS",
                  "Cs+", "ZN", "Zn2+", "F", "F-", "BR", "Br-", "I-"}
        }

        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                for api_code, residues in mapping.items():
                    if res_name in residues:
                        found.add(api_code)
        return list(found)

    def extract_water_coordinates(self, pdb_content: str) -> list[dict]:
        """
        Projdede PDB soubor a pro každou molekulu vody získá souřadnice [O, H1, H2].
        Vrací seznam slovníků obsahující identifikaci vody a její souřadnice.
        """
        waters = {}

        for line in pdb_content.splitlines():
            # Zajímá nás jen ATOM nebo HETATM záznam
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()

                # Zjištění, zda jde o vodu
                if res_name in ["HOH", "WAT", "SOL"]:
                    chain_id = line[21]
                    try:
                        res_seq = int(line[22:26].strip())
                    except ValueError:
                        continue

                    atom_name = line[12:16].strip()

                    # PDB formát má fixní sloupce pro souřadnice
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue

                    key = (chain_id, res_seq, res_name)
                    if key not in waters:
                        waters[key] = {}

                    # Detekce kyslíku a vodíků (názvy v PDB mohou být různé: O, OW, H1, HW1...)
                    if atom_name.startswith("O"):
                        waters[key]["O"] = [x, y, z]
                    elif atom_name.startswith("H"):
                        # Pokud ještě nemáme první vodík, uložíme jako H1, jinak jako H2
                        if "H1" not in waters[key]:
                            waters[key]["H1"] = [x, y, z]
                        else:
                            waters[key]["H2"] = [x, y, z]

        valid_waters = []
        # Vyfiltrujeme jen ty vody, které mají všechny 3 atomy (O, H, H)
        for (chain, res_seq, res_name), atoms in waters.items():
            if "O" in atoms and "H1" in atoms and "H2" in atoms:
                valid_waters.append({
                    "chain": chain,
                    "resseq": res_seq,
                    "res_name": res_name,
                    "crd": [atoms["O"], atoms["H1"], atoms["H2"]]
                })

        return valid_waters

    # --- NOVÁ FUNKCE PRO OPRAVU A VSTŘIKOVÁNÍ EP ---
    def reorder_and_inject_eps(self, pdb_content: str, ff_water_instance) -> str:
        """
        1. Rozdělí PDB na solut, ionty a vodu.
        2. Spočítá a rovnou vloží Extra Points k vodám pomocí water_models.
        3. Získá čisté PDB seřazené jako: Solut -> Ionty -> Voda (O, H1, H2, EP...).
        4. Přečísluje atomy.
        """
        solute_lines = []
        ion_lines = []
        waters = {}

        # Rozšířený seznam iontů podle mappingu
        ion_resnames = {"NA", "Na+", "CL", "Cl-", "K", "K+", "MG", "Mg2+", "CA", "Ca2+", "LI", "Li+", "RB", "Rb+", "CS",
                        "Cs+", "ZN", "Zn2+", "F", "F-", "BR", "Br-", "I", "I-"}

        # 1. ČTENÍ A ROZTŘÍDĚNÍ
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21]
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue

                res_key = (chain, resid, resname)

                if resname in ["HOH", "WAT", "SOL"]:
                    if res_key not in waters:
                        waters[res_key] = []
                    waters[res_key].append(line)
                elif resname in ion_resnames:
                    ion_lines.append(line)
                else:
                    solute_lines.append(line)
            elif line.startswith("CRYST1"):
                # Udržíme si informace o boxu na začátku
                solute_lines.insert(0, line)

        # Pomocné funkce pro parsování a formátování uvnitř metody
        def extract_coords(l: str) -> list[float]:
            return [float(l[30:38]), float(l[38:46]), float(l[46:54])]

        def format_ep_line(a_id: int, ep_n: str, r_n: str, c: str, r_id: int, crd: list[float]) -> str:
            x, y, z = crd
            return f"HETATM{a_id:>5} {ep_n:<4} {r_n:>3} {c}{r_id:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {ep_n[0]:>2}"

        # 2. ZÁPIS A GENEROVÁNÍ EPs S NOVÝM ČÍSLOVÁNÍM
        new_lines = []
        atom_id = 1

        # Zápis solutu
        for line in solute_lines:
            if line.startswith("CRYST1"):
                new_lines.append(line)
            else:
                new_lines.append(f"{line[:6]}{atom_id:>5}{line[11:]}")
                atom_id += 1

        # Zápis iontů (ionty jsou nyní před vodou)
        for line in ion_lines:
            new_lines.append(f"{line[:6]}{atom_id:>5}{line[11:]}")
            atom_id += 1

        # Zápis vod a výpočet EPs
        for (chain, resid, resname), atom_lines in waters.items():
            water_crd = []
            o_crd, h1_crd, h2_crd = None, None, None

            for line in atom_lines:
                atom_name = line[12:16].strip()
                crd = extract_coords(line)

                if atom_name.startswith("O"):
                    o_crd = crd
                elif atom_name.startswith("H"):
                    if not h1_crd:
                        h1_crd = crd
                    else:
                        h2_crd = crd

                new_lines.append(f"{line[:6]}{atom_id:>5}{line[11:]}")
                atom_id += 1

            # Jakmile máme atomy vody zapsané, spočítáme EP, pokud máme kompletní molekulu
            if o_crd and h1_crd and h2_crd and ff_water_instance:
                try:
                    eps = water_extra_points(ff_water_instance, [o_crd, h1_crd, h2_crd])
                    for i, ep_crd in enumerate(eps):
                        ep_name = f"EP{i + 1}"
                        ep_line = format_ep_line(atom_id, ep_name, resname, chain, resid, ep_crd)
                        new_lines.append(ep_line)
                        atom_id += 1
                except Exception as e:
                    logger.warning(f"Nepodařilo se přidat EP pro vodu {chain}{resid}: {e}")

        new_lines.append("END")
        return "\n".join(new_lines) + "\n"

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

            res_key = (chain_id, res_seq)

            if res_key not in seen_residues:
                seen_residues.add(res_key)
                if chain_id not in chains:
                    chains[chain_id] = []
                chains[chain_id].append(res_name)

    final_residues = []

    # Seznam běžných iontů pro lepší klasifikaci
    ion_resnames = {"NA", "Na+", "CL", "Cl-", "K", "K+", "MG", "Mg2+", "CA", "Ca2+", "LI", "Li+", "RB", "Rb+", "CS",
                    "Cs+", "ZN", "Zn2+", "F", "F-", "BR", "Br-", "I", "I-"}

    for chain_id, res_list in chains.items():
        # 1. Zjistíme, kde přesně začíná a končí Nukleová kyselina
        nuc_indices = []
        for idx, name in enumerate(res_list):
            # Detekce RNA/DNA
            if (len(name) == 1 and name in ['A', 'C', 'G', 'U', 'T']) or name.startswith('R') or name.startswith('D'):
                nuc_indices.append(idx)

        first_nuc = nuc_indices[0] if nuc_indices else -1
        last_nuc = nuc_indices[-1] if nuc_indices else -1

        for i, res_name in enumerate(res_list):
            mol_type = "R"  # Výchozí

            if res_name in ['HOH', 'WAT', 'SOL']:
                mol_type = "W"
            elif res_name in ion_resnames:
                mol_type = "I"

            is_nucleotide = i in nuc_indices
            base_name = res_name

            if is_nucleotide and len(base_name) == 1:
                base_name = 'R' + base_name

            if is_nucleotide:
                # Přiřadíme 5' a 3' konce bezpečně POUZE na první a poslední nukleotid!
                if i == first_nuc and not base_name.endswith('5'):
                    final_resn = base_name + "5"
                elif i == last_nuc and not base_name.endswith('3'):
                    final_resn = base_name + "3"
                else:
                    final_resn = base_name
            else:
                final_resn = base_name

            final_residues.append({
                "resn": final_resn,
                "mol_type": mol_type,
                "chain": chain_id  # Předáme dál informaci o řetězci
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