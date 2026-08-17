from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Set, Any

from app.utils.adams4sims_processing_library.utils.alias import resn_alias, name_alias


@lru_cache(maxsize=1)
def load_converting_dictionary() -> dict:
    """
    NAČTE KONVERZNÍ SLOVNÍK ZE SOUBORU JSON V KOŘENOVÉM ADRESÁŘI A ZAJIŠŤUJE JEHO CACHOVÁNÍ PRO RYCHLÝ PŘÍSTUP.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(os.path.dirname(current_dir))
    dict_path = os.path.join(root_dir, "data", "converting_dictionary.json")

    try:
        with open(dict_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _parse_residues_from_pdb(pdb_text: str, chain: Optional[str]) -> List[Tuple[str, int, str, str, List[str]]]:
    """
    EXTRAHUJE DATA O REZIDUÍCH A JEJICH ATOMECH Z PDB FORMÁTU, PŘIČEMŽ PROVÁDÍ ALIASING NÁZVŮ ATOMŮ PODLE SLOVNÍKU.
    """
    residue_data = {}
    ordered_keys = []

    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        resname = line[17:20].strip()
        ch = (line[21] or "").strip() or "?"
        resseq_raw = line[22:26].strip()
        icode = (line[26] or " ").strip()
        atom_name_raw = line[12:16].strip()

        atom_name = name_alias(resname, atom_name_raw)

        if not resseq_raw:
            continue

        try:
            resseq = int(resseq_raw)
        except ValueError:
            continue

        if chain and ch != chain:
            continue

        key = (ch, resseq, icode)
        if key not in residue_data:
            residue_data[key] = {"resname": resname, "atoms": []}
            ordered_keys.append(key)

        if atom_name not in residue_data[key]["atoms"]:
            residue_data[key]["atoms"].append(atom_name)

    ordered_keys.sort(key=lambda x: (x[0], x[1], x[2]))

    return [
        (k[0], k[1], k[2], residue_data[k]["resname"], residue_data[k]["atoms"])
        for k in ordered_keys
    ]


def _infer_group(resname: str, conv: Dict) -> Optional[str]:
    """
    IDENTIFIKUJE CHEMICKOU KATEGORII REZIDUA PROHLEDÁVÁNÍM KLÍČŮ V KONVERZNÍM SLOVNÍKU NEBO POMOCÍ PREFIXŮ.
    """
    aliased = resn_alias(resname)
    for category in conv.keys():
        if isinstance(conv[category], dict):
            if resname in conv[category] or aliased in conv[category]:
                return category

    if resname.startswith("D"):
        return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"):
        return "R"

    return None


def _get_res_def(group: Optional[str], ff_name: str, conv: Dict) -> Optional[Dict]:
    """
    VYHLEDÁ DEFINICI REZIDUA (ATOMY A KONEKTIVITU) V KONKRÉTNÍ KATEGORII NEBO PROHLEDÁNÍM CELÉHO SLOVNÍKU.
    """
    if group and group in conv and ff_name in conv[group]:
        return conv[group][ff_name]

    for category in conv.values():
        if isinstance(category, dict) and ff_name in category:
            return category[ff_name]
    return None


def _pick_variant(group: Optional[str], pdb_resname: str, atoms: List[str], conv: Dict, terminal: str) -> Tuple[
    Optional[str], bool, Optional[str]]:
    """
    URČUJE VHODNOU FF VARIANTU NA ZÁKLADĚ TERMINÁLNÍ POZICE A PŘÍTOMNÝCH ATOMŮ (NAPŘ. PRO HISTIDIN).
    """
    search_group = resn_alias(pdb_resname)

    # 1. INTELIGENTNÍ DETEKCE HISTIDINU
    # Pokud máme v PDB 'HIS' (nebo aliasovaný HID, HIE, HIP), zkusíme variantu potvrdit podle reálných vodíků
    if pdb_resname == "HIS" or search_group in ["HID", "HIE", "HIP"]:
        # HIP má oba vodíky (HD1 na delta-dusíku a HE2 na epsilon-dusíku)
        if "HD1" in atoms and "HE2" in atoms:
            search_group = "HIP"
        # HIE má vodík jen na epsilon dusíku
        elif "HE2" in atoms:
            search_group = "HIE"
        # HID má vodík na delta dusíku (nebo je to výchozí stav, pokud vodíky chybí)
        else:
            search_group = "HID"

    candidates: List[str] = []

    # 2. SESTAVENÍ KANDIDÁTŮ PRO KONCE ŘETĚZCŮ
    if group == "P":
        # Proteiny (skupina P) mají N-konec a C-konec (např. NALA, CALA, NHID, CHID)
        if terminal == "5":
            candidates.append(f"N{search_group}")
        elif terminal == "3":
            candidates.append(f"C{search_group}")
    else:
        # Nukleové kyseliny mají koncovky 5 a 3 (např. RU5, RU3)
        if terminal == "5":
            candidates.append(f"{search_group}5")
        elif terminal == "3":
            candidates.append(f"{search_group}3")

    # Vždy přidáme jako fallback základní variantu uprostřed řetězce
    candidates.append(search_group)

    # 3. OVĚŘENÍ PROTI SLOVNÍKU
    for k in candidates:
        res_def = _get_res_def(group, k, conv)
        if res_def:
            return k, True, search_group

    return search_group, False, search_group


def _reterminate_as_gap_end(token: Dict[str, Any], conv: Dict, warnings: List[str], next_chain: str, next_resseq: int) -> None:
    """
    Přepíše už zapsaný token (poslední residuum PŘED sekvenční dírou v hlavním řetězci)
    na jeho umělou terminální variantu (C-konec pro protein / 3'-konec pro RNA-DNA).

    Bez tohoto kroku by ff_resname zůstal ve "středové" variantě, která v konverzním
    slovníku očekává navazující sousední residuum - a downstream builder (app/builder)
    by pak mohl přes chybějící úsek vytvořit nesmyslnou vazbu. Builder sám gap
    nedostavuje (viz INTEGRATION_CONTRACT.md), takže tohle rozhodnutí musí padnout tady.
    """
    group = token["group"]
    new_ff_resname, new_known, _ = _pick_variant(group, token["pdb_resname"], token["atoms"], conv, "3")
    token["ff_resname"] = new_ff_resname
    token["known"] = new_known
    token["missing_atoms"] = _check_missing_atoms(group, new_ff_resname, token["atoms"], conv)
    conn_info = _check_connectivity_integrity(group, new_ff_resname, token["atoms"], conv)
    token["is_broken"] = conn_info["is_broken"]
    token["connectivity_parts"] = conn_info["components"]
    token["terminus_reason"] = "gap"

    label = "C-terminus" if group == "P" else "3'-terminus"
    warnings.append(
        f"{token['chain']}:{token['resseq']} ({token['pdb_resname']}) treated as artificial {label} "
        f"— sequence gap before {next_chain}:{next_resseq}."
    )


def build_sequence_tokens(pdb_text: str, chain: Optional[str] = None, fill_gaps: bool = True):
    """
    SESTAVUJE KOMPLETNÍ SEZNAM TOKENS PRO DANÝ ŘETĚZEC, PROVÁDÍ ANALÝZU VARIANT A DETEKCI CHYBĚJÍCÍCH ČÁSTÍ STRUKTURY.
    """
    conv = load_converting_dictionary()
    all_residues = _parse_residues_from_pdb(pdb_text, chain)

    if not all_residues:
        return {"chain": chain, "tokens": [], "warnings": ["No residues found in PDB."]}

    residues_to_process = [r for r in all_residues if r[0] == chain] if chain else all_residues

    chains_dict = {}
    for r in residues_to_process:
        ch_id = r[0]
        if ch_id not in chains_dict:
            chains_dict[ch_id] = []
        chains_dict[ch_id].append(r)

    tokens = []
    warnings = []
    global_pos = 0

    for ch_id, residues in chains_dict.items():
        main_chain = []
        ligands = []

        for r in residues:
            resname = r[3]
            group = _infer_group(resname, conv)
            if group in ["R", "D", "P"]:
                main_chain.append(r)
            else:
                ligands.append(r)

        main_chain.sort(key=lambda x: (x[1], x[2]))
        first_main_seq = main_chain[0][1] if main_chain else None
        last_main_seq = main_chain[-1][1] if main_chain else None
        processed_ordered = main_chain + ligands
        prev_resseq = None
        prev_main_token = None

        for ch, resseq, icode, resname, atoms in processed_ordered:
            is_main = any(r[1] == resseq and r[3] == resname for r in main_chain)
            after_gap = False

            if fill_gaps and is_main and prev_resseq is not None and resseq > prev_resseq + 1:
                # OPRAVA: Zkontroluj, zda GAP je OPRAVDU prázdný nebo tam jen je neznámé reziduum
                gap_found = False
                for missing_seq in range(prev_resseq + 1, resseq):
                    # Hledej, zda existuje JAKÉKOLI reziduum se sekvencí missing_seq v PDB
                    residue_exists_in_pdb = any(r[1] == missing_seq and r[0] == ch for r in residues)

                    # Jen pokud OPRAVDU chybí v PDB -> vytvoř GAP token
                    if not residue_exists_in_pdb:
                        gap_found = True
                        global_pos += 1
                        tokens.append({
                            "position": global_pos, "chain": ch, "resseq": None, "icode": None,
                            "pdb_resname": "0", "is_gap": True, "group": None, "ff_resname": None,
                            "known": False, "atoms": [], "missing_atoms": []
                        })

                if gap_found:
                    # Residuum před dírou i residuum za dírou se stávají uměle
                    # terminálními, ať se přes chybějící úsek nepočítá žádná vazba.
                    after_gap = True
                    if prev_main_token is not None:
                        _reterminate_as_gap_end(prev_main_token, conv, warnings, ch, resseq)

            global_pos += 1
            group = _infer_group(resname, conv)

            terminal = ""
            terminus_reason = None
            if is_main:
                if after_gap:
                    terminal = "5"
                    terminus_reason = "gap"
                elif resseq == first_main_seq:
                    terminal = "5"
                    terminus_reason = "chain_end"
                elif resseq == last_main_seq:
                    terminal = "3"
                    terminus_reason = "chain_end"

            ff_resname, known, search_group = _pick_variant(group, resname, atoms, conv, terminal)
            missing_atoms = _check_missing_atoms(group, ff_resname, atoms, conv)

            if not known:
                warnings.append(f"Unknown residue '{resname}' at {ch}:{resseq}{icode or ''}")
            elif missing_atoms:
                warnings.append(f"Incomplete residue '{resname}' at {ch}:{resseq}: Missing {missing_atoms}")

            if after_gap:
                label = "N-terminus" if group == "P" else "5'-terminus"
                warnings.append(
                    f"{ch}:{resseq} ({resname}) treated as artificial {label} "
                    f"— sequence gap before this residue."
                )

            conn_info = _check_connectivity_integrity(group, ff_resname, atoms, conv)

            token = {
                "position": global_pos,
                "chain": ch,
                "resseq": resseq,
                "icode": icode or "",
                "pdb_resname": resname,
                "is_gap": False,
                "group": group,
                "ff_resname": ff_resname,
                "known": known,
                "atoms": atoms,
                "missing_atoms": missing_atoms,
                "is_broken": conn_info["is_broken"],
                "connectivity_parts": conn_info["components"],
                "terminus_reason": terminus_reason
            }
            tokens.append(token)

            if is_main:
                prev_resseq = resseq
                prev_main_token = token

    return {
        "chains": {
            ch_id: {
                "chain": ch_id,
                "tokens": [t for t in tokens if t["chain"] == ch_id],
                "warnings": [w for w in warnings if f"chain {ch_id}" in w or f"{ch_id}:" in w]
            }
            for ch_id in chains_dict.keys()
        }
    }


def _check_missing_atoms(group: Optional[str], ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    """
    POROVNÁVÁ SEZNAM ATOMŮ Z PDB SE ŠABLONOU V KONVERZNÍM SLOVNÍKU A VRACÍ SEZNAM VŠECH CHYBĚJÍCÍCH ELEMENTŮ.
    """
    res_def = _get_res_def(group, ff_name, conv)
    if not res_def:
        return []

    required_atoms = set(res_def.get("atom", {}).keys())
    actual_atoms = set(atoms)

    return sorted([a for a in required_atoms if a not in actual_atoms])


def _check_connectivity_integrity(group: Optional[str], ff_name: str, atoms: List[str], conv: Dict) -> Dict[str, Any]:
    """
    ANALYZUJE GEOMETRICKOU INTEGRITU REZIDUA POMOCÍ GRAFU KONEKTIVITY A IDENTIFIKUJE IZOLOVANÉ SKUPINY ATOMŮ.
    """
    res_def = _get_res_def(group, ff_name, conv)
    if not res_def:
        return {"is_broken": False, "components": [atoms] if atoms else []}

    conn_map = res_def.get("connectivity", {})
    if not conn_map:
        return {"is_broken": False, "components": [atoms] if atoms else []}

    present_atoms = set(atoms)
    graph = {atom: set() for atom in present_atoms}
    for u, neighbors in conn_map.items():
        if u in present_atoms:
            for v in neighbors:
                if v in present_atoms:
                    graph[u].add(v)
                    graph[v].add(u)

    visited = set()
    components = []
    for start_node in sorted(list(present_atoms)):
        if start_node not in visited:
            component = []
            queue = [start_node]
            visited.add(start_node)
            while queue:
                u = queue.pop(0)
                component.append(u)
                for v in graph.get(u, []):
                    if v not in visited:
                        visited.add(v)
                        queue.append(v)
            components.append(sorted(component))

    return {
        "is_broken": len(components) > 1,
        "components": components
    }


from typing import Dict, Any


def analyze_pdb_altlocs(pdb_text: str) -> Dict[str, Any]:
    """
    PROJDE PDB SOUBOR A IDENTIFIKUJE VŠECHNY ALTERNATIVNÍ POZICE (ALTLOCS),
    JEJICH OBSAZENOST A B-FAKTOR. VRACÍ STRUKTUROVANÝ DICT (JSON) PRO FRONTEND.
    NAVÍC ANALYZUJE KONEKTIVITU (BLOKY NA SEBE NAVAZUJÍCÍCH AMINOKYSELIN)
    A DOPORUČUJE NEJLEPŠÍ TRASU PRO ZACHOVÁNÍ PEPTIDOVÉ VAZBY.

    NOVĚ: DETEKUJE PŘÍTOMNOST VÍCE MODELŮ A SYMETRIE (REMARK 350).
    """
    altloc_data = {}
    models = []
    has_symmetry = False

    # Projdeme soubor řádek po řádku
    for line in pdb_text.splitlines():

        # --- NOVÉ: Detekce více modelů ---
        if line.startswith("MODEL "):
            try:
                # Ořízneme slovo "MODEL" a zkusíme získat číslo
                model_num = int(line[6:].strip())
                if model_num not in models:
                    models.append(model_num)
            except ValueError:
                pass
            continue

        # --- NOVÉ: Detekce Biological Assembly (Symetrie) ---
        # Hledáme řádek REMARK 350, který obsahuje transformační matici BIOMT
        if line.startswith("REMARK 350") and "BIOMT" in line:
            has_symmetry = True
            continue

        # --- PŮVODNÍ: Detekce AltLocs ---
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Index 16 je sloupec 17 v PDB (AltLoc)
            alt_loc = line[16]

            # Pokud to není mezera, našli jsme alternativní pozici
            if alt_loc != ' ':
                chain = line[21].strip() or "?"
                resseq_raw = line[22:26].strip()
                resname = line[17:20].strip()

                try:
                    resseq = int(resseq_raw)
                except ValueError:
                    continue

                # Sloupce 55-60 (index 54:60) pro Occupancy
                try:
                    occupancy = float(line[54:60].strip())
                except ValueError:
                    occupancy = 1.0

                # Sloupce 61-66 (index 60:66) pro B-faktor
                try:
                    b_factor = float(line[60:66].strip())
                except ValueError:
                    b_factor = 0.0

                # Unikátní klíč pro konkrétní aminokyselinu
                key = (chain, resseq, resname)

                if key not in altloc_data:
                    altloc_data[key] = {}

                if alt_loc not in altloc_data[key]:
                    altloc_data[key][alt_loc] = {
                        "occupancy": round(occupancy * 100, 1),
                        "bFactor": b_factor
                    }

    # Nyní to přetavíme do pole pro Frontend
    result_residues = []
    for (chain, resseq, resname), alt_locs in altloc_data.items():
        if len(alt_locs) > 0:
            result_residues.append({
                "chain": chain,
                "resseq": resseq,
                "resname": resname,
                "altLocs": alt_locs
            })

    # Seřadíme podle řetězce a čísla zbytku
    result_residues.sort(key=lambda x: (x["chain"], x["resseq"]))

    if result_residues:
        blocks = []
        current_block = [result_residues[0]]

        # 1. Seskládáme rezidua do bloků (pokud po sobě následují v číslování)
        for i in range(1, len(result_residues)):
            prev_res = current_block[-1]
            curr_res = result_residues[i]

            if curr_res["chain"] == prev_res["chain"] and curr_res["resseq"] == prev_res["resseq"] + 1:
                current_block.append(curr_res)
            else:
                blocks.append(current_block)
                current_block = [curr_res]

        blocks.append(current_block)

        # 2. Pro každý blok určíme vítěznou trasu
        for block in blocks:
            paths_stats = {}

            # Nasbíráme součty z celého bloku
            for res in block:
                for alt, info in res["altLocs"].items():
                    if alt not in paths_stats:
                        paths_stats[alt] = {"occ": 0.0, "bfact": 0.0, "count": 0}
                    paths_stats[alt]["occ"] += info["occupancy"]
                    paths_stats[alt]["bfact"] += info["bFactor"]
                    paths_stats[alt]["count"] += 1

            best_alt = None
            best_occ = -1.0
            best_bfact = float('inf')

            # Spočítáme průměry a vybereme vítěze pro daný řetězec
            for alt, stats in paths_stats.items():
                avg_occ = stats["occ"] / stats["count"]
                avg_bfact = stats["bfact"] / stats["count"]

                # Vyhrává vyšší occupancy. Pokud je 50 na 50 (remíza), vyhrává nižší B-faktor!
                if avg_occ > best_occ:
                    best_occ = avg_occ
                    best_bfact = avg_bfact
                    best_alt = alt
                elif avg_occ == best_occ:
                    if avg_bfact < best_bfact:
                        best_bfact = avg_bfact
                        best_alt = alt

            # 3. Zapíšeme vítěznou volbu do dat pro frontend
            for res in block:
                # Ošetření: zkontrolujeme, zda tuto trasu reziduum reálně obsahuje
                if best_alt in res["altLocs"]:
                    res["recommended_alt"] = best_alt
                else:
                    # Fallback na lokální maximum, pokud by vítězná trasa u tohoto rezidua nebyla
                    local_best = max(res["altLocs"].keys(),
                                     key=lambda k: (res["altLocs"][k]["occupancy"], -res["altLocs"][k]["bFactor"]))
                    res["recommended_alt"] = local_best
    # -------------------------------------------------------------------

    # Vracíme obohacený JSON
    return {
        "models": models,  # Přidáno: Pole s čísly modelů (např. [1, 2, 3])
        "hasSymmetry": has_symmetry,  # Přidáno: True/False, pokud existuje BIOMT matice
        "hasAltLocs": len(result_residues) > 0,
        "residues": result_residues
    }


def clean_pdb_altlocs(pdb_text: str, user_selection: dict) -> str:
    """
    PROJDE PDB SOUBOR A SMAŽE VŠECHNY ALTERNATIVNÍ ATOMY I JEJICH ANISOU ZÁZNAMY,
    KTERÉ UŽIVATEL NEVYBRAL. VYBRANÝM ATOMŮM ODSTRANÍ ALTLOC INDIKÁTOR
    A NASTAVÍ OBSAZENOST (OCCUPANCY) ZPĚT NA 1.00.
    """
    cleaned_lines = []

    for line in pdb_text.splitlines():
        # OPRAVA 1: Přidali jsme kontrolu i pro ANISOU řádky
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
            alt_loc = line[16]

            # Pokud má atom/anisou alternativní pozici
            if alt_loc != ' ':
                chain = line[21].strip() or "?"
                resseq_raw = line[22:26].strip()
                resname = line[17:20].strip()

                # Vygenerujeme stejný unikátní klíč, jaký nám posílá frontend (např. "A_45_TYR")
                unique_key = f"{chain}_{resseq_raw}_{resname}"

                # Pokud uživatel pro tento zbytek provedl volbu
                if unique_key in user_selection:
                    chosen_altloc = user_selection[unique_key]

                    # Pokud se AltLoc atomu neshoduje s volbou uživatele -> Smazat řádek (nepřidá se do cleaned_lines)
                    if alt_loc != chosen_altloc:
                        continue
                    else:
                        # Pokud je to ten správný atom/anisou, nahradíme jeho AltLoc znak ('A'/'B') za mezeru
                        line = line[:16] + ' ' + line[17:]

                        # OPRAVA 2: Pokud jde o ATOM/HETATM, musíme vrátit obsazenost (Occupancy) na 1.00
                        # (ANISOU řádky pole Occupancy nemají, proto ta podmínka)
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if len(line) >= 60:
                                # Obsazenost je v PDB formátu přesně na pozicích 55-60 (index 54:60)
                                line = line[:54] + "  1.00" + line[60:]

        # Přidáme upravený (nebo nedotčený) řádek do nového čistého souboru
        cleaned_lines.append(line)

    return '\n'.join(cleaned_lines)


def process_structure(pdb_text: str, target_model: int, apply_symmetry: bool, selection: dict) -> str:
    """
    Kombinovaná funkce pro kompletní fyzickou přípravu PDB souboru:
    1. Ponechá pouze vybraný MODEL (např. NMR ensemble).
    2. Vymaže nevybrané alternativní pozice a upraví obsazenost na 1.0.
    3. Pokud apply_symmetry=True, vybuduje plnou biologickou jednotku (Biological Assembly)
       pomocí BIOMT matic a matematicky dopočítá atomy.
    """
    lines = pdb_text.splitlines()

    # --- KROK 1: Přečtení BIOMT matic z REMARK 350 ---
    matrices = {}
    for line in lines:
        if line.startswith("REMARK 350   BIOMT"):
            try:
                row_idx = int(line[18:19])  # Řádek matice (1, 2 nebo 3)
                mat_num = int(line[20:24].strip())
                parts = line[24:].split()
                if len(parts) >= 4:
                    if mat_num not in matrices:
                        matrices[mat_num] = [[0.0] * 4 for _ in range(3)]
                    matrices[mat_num][row_idx - 1] = [float(x) for x in parts[:4]]
            except Exception:
                continue

    biomt_list = list(matrices.values())
    # Pokud v PDB není BIOMT matice, dodáme základní (Identitu - x*1 = x)
    if not biomt_list:
        biomt_list = [[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]]

    # --- KROK 2: Extrakce konkrétního modelu a AltLocs filtrace ---
    in_model = False
    has_models = any(l.startswith("MODEL ") for l in lines)

    # Pokud soubor nemá vůbec tagy MODEL (typická krystalografie), čteme atomy rovnou
    if not has_models:
        in_model = True

    base_atoms = []
    for line in lines:
        if line.startswith("MODEL "):
            try:
                current_model = int(line[6:].strip())
                in_model = (current_model == target_model)
            except ValueError:
                pass
            continue
        elif line.startswith("ENDMDL"):
            in_model = False
            continue

        if in_model and (line.startswith("ATOM  ") or line.startswith("HETATM")):
            alt_loc = line[16]
            chain = line[21]
            resseq = line[22:26].strip()
            resname = line[17:20].strip()

            # Vytvoření unikátního klíče z frontendu
            key = f"{chain.strip() or '?'}_{resseq}_{resname}"

            if alt_loc != ' ':
                chosen_alt = selection.get(key)
                # Pokud na této pozici sedí AltLoc a uživatel vybral něco jiného, přeskočíme atom
                if chosen_alt and alt_loc != chosen_alt:
                    continue
                # Úprava řádku: Vymazání AltLoc písmene (mezera) a fixace obsazenosti na 1.00
                line = line[:16] + ' ' + line[17:54] + "  1.00" + line[60:]

            base_atoms.append(line)

    # --- KROK 3: Fyzické budování biologické jednotky a přečíslování ---
    final_lines = []

    if apply_symmetry and len(biomt_list) > 1:
        # Pokud násobíme řetězce (A), musíme zkopírované kusy přejmenovat na nová písmena (B, C atd.)
        unique_chains = []
        for line in base_atoms:
            ch = line[21]
            if ch not in unique_chains:
                unique_chains.append(ch)

        chain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
        used_chains = set(unique_chains)
        available_chains = [c for c in chain_alphabet if c not in used_chains]

        # Mapa pro dynamické přidělování názvů řetězců
        matrix_chain_maps = {0: {ch: ch for ch in unique_chains}}  # 1. matice nemění název
        for i in range(1, len(biomt_list)):
            matrix_chain_maps[i] = {}
            for ch in unique_chains:
                if available_chains:
                    new_c = available_chains.pop(0)
                    matrix_chain_maps[i][ch] = new_c
                else:
                    matrix_chain_maps[i][ch] = ch  # Došly znaky, použijeme starý

        # Násobení: Vezmeme každou transformační matici a převalíme přes ni všechny vyčištěné atomy
        atom_serial = 1
        for i, matrix in enumerate(biomt_list):
            for line in base_atoms:
                try:
                    # Rozparsování původních souřadnic
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # Maticové násobení + translace
                    nx = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3]
                    ny = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3]
                    nz = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]

                    # Nové ID řetězce
                    orig_ch = line[21]
                    new_ch = matrix_chain_maps[i].get(orig_ch, orig_ch)

                    # Složení nového řádku zpět podle striktního PDB formátu
                    new_line = (
                            line[:6] +
                            f"{atom_serial:5d}" +
                            line[11:21] +
                            new_ch +
                            line[22:30] +
                            f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" +
                            line[54:]
                    )
                    final_lines.append(new_line)
                    atom_serial += 1
                except Exception:
                    final_lines.append(line)
    else:
        # Přečíslování atomů
        atom_serial = 1
        for line in base_atoms:
            try:
                new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
                final_lines.append(new_line)
                atom_serial += 1
            except Exception:
                final_lines.append(line)

    final_lines.append("END")
    return "\n".join(final_lines)