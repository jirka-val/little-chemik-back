from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Set, Any

from app.utils.alias import resn_alias, name_alias


@lru_cache(maxsize=1)
def load_converting_dictionary() -> dict:
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dict_path = os.path.join(base_dir, "data", "converting_dictionary.json")

    try:
        with open(dict_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def _parse_residues_from_pdb(pdb_text: str, chain: Optional[str]) -> List[Tuple[str, int, str, str, List[str]]]:
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
    aliased = resn_alias(resname)
    for category in conv.keys():
        if resname in conv[category] or aliased in conv[category]:
            return category

    if resname.startswith("D"):
        return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"):
        return "R"

    return None


def _pick_variant(group: Optional[str], pdb_resname: str, conv: Dict, terminal: str) -> Tuple[
    Optional[str], bool, Optional[str]]:
    """
    Vrací: (ff_resname, known, search_group)
    search_group je klíč v JSONu (např. 'RU'), pod kterým se mají hledat atomy.
    """
    # Normalizace (U -> RU)
    aliased = resn_alias(pdb_resname)

    # Určení základního klíče pro hledání v JSONu
    # Pokud je v JSONu klíč 'RU', použijeme ho jako search_group
    search_group = aliased if aliased in conv else pdb_resname

    candidates: List[str] = []
    # Definice konců (RU5, RU3)
    if terminal == "5":
        candidates += [f"{search_group}5"]
    elif terminal == "3":
        candidates += [f"{search_group}3"]

    candidates += [search_group]

    # Hledáme shodu v JSONu
    for k in candidates:
        if k in conv:
            return k, True, k
        # Fallback pokud je to vnořené v kategorii (např. R -> RU)
        if group and group in conv and k in conv[group]:
            return k, True, k

    # Pokud nic nenajde, vrátí aspoň název, ale known=False
    return search_group, search_group in conv, search_group


def build_sequence_tokens(pdb_text: str, chain: Optional[str] = None, fill_gaps: bool = True):
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

        for ch, resseq, icode, resname, atoms in processed_ordered:
            is_main = any(r[1] == resseq and r[3] == resname for r in main_chain)

            if fill_gaps and is_main and prev_resseq is not None and resseq > prev_resseq + 1:
                for missing_seq in range(prev_resseq + 1, resseq):
                    global_pos += 1
                    tokens.append({
                        "position": global_pos, "chain": ch, "resseq": None, "icode": None,
                        "pdb_resname": "0", "is_gap": True, "group": None, "ff_resname": None,
                        "known": False, "atoms": [], "missing_atoms": []
                    })

            global_pos += 1
            group = _infer_group(resname, conv)

            terminal = ""
            if is_main:
                if resseq == first_main_seq:
                    terminal = "5"
                elif resseq == last_main_seq:
                    terminal = "3"

            # Najdi v build_sequence_tokens tento blok a nahraď ho:

            # 1. Získáme variantu a hlavně search_group
            ff_resname, known, search_group = _pick_variant(group, resname, conv, terminal)

            # 2. ZMĚNA: Kontrolujeme missing_atoms vždy, když máme search_group,
            # nehledě na to, jestli 'known' dopadlo jako True/False
            missing_atoms = _check_missing_atoms(search_group, ff_resname, atoms, conv) if search_group else []

            if not known:
                warnings.append(f"Unknown residue '{resname}' at {ch}:{resseq}{icode or ''}")
            elif missing_atoms:
                warnings.append(f"Incomplete residue '{resname}' at {ch}:{resseq}: Missing {missing_atoms}")

            # 3. Kontrola konektivity taky se search_group
            conn_info = _check_connectivity_integrity(search_group, ff_resname, atoms, conv)

            tokens.append({
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
                "connectivity_parts": conn_info["components"]
            })

            if is_main:
                prev_resseq = resseq

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


def _check_missing_atoms(search_group: str, ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    if not search_group or search_group not in conv:
        return []

    # Najde definici v plochém JSONu (RU5) nebo vnořeném (RU -> RU5)
    res_def = conv.get(ff_name) or conv[search_group].get(ff_name, conv[search_group])

    if not isinstance(res_def, dict): return []

    # Získání očekávaných atomů
    required_atoms_dict = res_def.get("atom", {})
    if not required_atoms_dict and "connectivity" in res_def:
        required_atoms_dict = res_def.get("atom", {})  # Fallback

    required_atoms = set(required_atoms_dict.keys())
    actual_atoms_set = set(atoms)

    return sorted([a for a in required_atoms if a not in actual_atoms_set])


def _check_connectivity_integrity(search_group: str, ff_name: str, atoms: List[str], conv: Dict) -> Dict[str, Any]:
    if not search_group or search_group not in conv:
        return {"is_broken": False, "components": []}

    res_def = conv.get(ff_name) or conv[search_group].get(ff_name, conv[search_group])
    if not isinstance(res_def, dict):
        return {"is_broken": False, "components": []}

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