# app/services/analysis_service.py

from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Set, Any

from app.utils.alias import resn_alias, name_alias


@lru_cache(maxsize=1)
def load_converting_dictionary() -> Dict:
    candidates = [
        os.path.join(os.getcwd(), "converting_dictionary.json"),
        os.path.join(os.getcwd(), "app", "resources", "converting_dictionary.json"),
        os.path.join(os.getcwd(), "app", "data", "converting_dictionary.json"),
        os.path.join(os.getcwd(),
                     "jirka-val/little-chemik-back/little-chemik-back-a58f320e7b9715efb9053322f92c5fed67c5cd19/converting_dictionary.json"),
    ]
    for path in candidates:
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    raise FileNotFoundError("converting_dictionary.json not found.")


def _parse_residues_from_pdb(pdb_text: str, chain: Optional[str]) -> List[Tuple[str, int, str, str, List[str]]]:
    """
    Vrací unikátní residua se seznamem jejich atomů: (chain, resseq, icode, resname, atoms)
    """
    residue_data = {}  # (ch, resseq, icode) -> {"resname": str, "atoms": []}
    ordered_keys = []

    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        resname = line[17:20].strip()
        ch = (line[21] or "").strip() or "?"
        resseq_raw = line[22:26].strip()
        icode = (line[26] or " ").strip()
        atom_name_raw = line[12:16].strip()

        # Normalizace jména atomu podle aliasů (např. O5* -> O5')
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
    Zjišťuje kategorii rezidua (R, D, P, W3, I1, atd.) přímo ze slovníku.
    Prochází všechny dostupné kategorie a využívá aliasy pro přesnější detekci.
    """
    # Získání aliasu (např. 'A' -> 'RA') pro porovnání se slovníkem
    aliased = resn_alias(resname)

    # Procházíme všechny hlavní klíče ve slovníku (automatická podpora celé šíře)
    for category in conv.keys():
        # Kontrolujeme původní název i alias
        if resname in conv[category] or aliased in conv[category]:
            return category

    # Bezpečný fallback pro nukleové kyseliny, pokud nejsou přímo v kategorii
    if resname.startswith("D"):
        return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"):
        return "R"

    return None

def _validate_connectivity(group: str, ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    """
    Kontroluje, zda jsou přítomny všechny atomy vyžadované mapou konektivity.
    """
    if not group or not ff_name or group not in conv or ff_name not in conv[group]:
        return []

    res_def = conv[group][ff_name]
    conn_map = res_def.get("connectivity", {})

    # Množina všech atomů, které konektivita očekává (klíče i sousedé)
    required_atoms = set(conn_map.keys())
    for neighbors in conn_map.values():
        required_atoms.update(neighbors)

    # Vrátíme ty, které v PDB chybí
    return [a for a in required_atoms if a not in atoms]


def _pick_variant(group: Optional[str], pdb_resname: str, conv: Dict, terminal: str) -> Tuple[Optional[str], bool]:
    if not group or group not in conv:
        return None, False

    aliased = resn_alias(pdb_resname)
    candidates: List[str] = []

    # Mapování konců přesně podle tvého JSONu (např. RU + 3 = RU3)
    if terminal == "5":
        candidates += [f"{aliased}5", f"{pdb_resname}5"]
    elif terminal == "3":
        candidates += [f"{aliased}3", f"{pdb_resname}3"]

    # Základní názvy jako fallback
    candidates += [aliased, pdb_resname]

    for k in candidates:
        if k in conv[group]:
            return k, True

    # Pokud nic nenajde, vrátí první pokus, ale označí jako neznámé
    return candidates[0] if candidates else None, False


def build_sequence_tokens(pdb_text: str, chain: Optional[str] = None, fill_gaps: bool = True):
    conv = load_converting_dictionary()
    all_residues = _parse_residues_from_pdb(pdb_text, chain)

    if not all_residues:
        return {"chain": chain, "tokens": [], "warnings": ["No residues found in PDB."]}

    # Filtrace řetězců
    if chain is not None:
        residues_to_process = [r for r in all_residues if r[0] == chain]
    else:
        residues_to_process = all_residues

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
        # 1. ROZDĚLENÍ: Najdeme hlavní polymerní řetězec a ligandy zvlášť
        main_chain = []
        ligands = []

        for r in residues:
            resname = r[3]
            group = _infer_group(resname, conv)
            # Skupiny R (RNA), D (DNA) a P (Proteiny) tvoří hlavní řetězec
            if group in ["R", "D", "P"]:
                main_chain.append(r)
            else:
                ligands.append(r)

        # Seřadíme hlavní řetězec podle číslování PDB
        main_chain.sort(key=lambda x: (x[1], x[2]))

        # 2. KONCE: Určíme reálný začátek a konec pouze z hlavního řetězce
        first_main_seq = main_chain[0][1] if main_chain else None
        last_main_seq = main_chain[-1][1] if main_chain else None

        # 3. POŘADÍ: Zpracujeme nejdříve řetězec (s mezerami) a pak přidáme ligandy na konec
        processed_ordered = main_chain + ligands

        prev_resseq = None

        for ch, resseq, icode, resname, atoms in processed_ordered:
            # Mezery vyplňujeme pouze v rámci souvislého hlavního řetězce
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

            # 4. TERMINÁLY: Přiřadíme 5/3 pouze pokud je to součást hlavního řetězce a je na kraji
            terminal = ""
            if is_main:
                if resseq == first_main_seq:
                    terminal = "5"
                elif resseq == last_main_seq:
                    terminal = "3"

            # Vyhledání varianty (např. RU3) v JSONu - bez diakritiky
            ff_resname, known = _pick_variant(group, resname, conv, terminal)

            missing_atoms = _check_missing_atoms(group, ff_resname, atoms, conv) if known else []

            if not known:
                warnings.append(f"Unknown residue '{resname}' at {ch}:{resseq}{icode or ''}")
            elif missing_atoms:
                warnings.append(f"Incomplete residue '{resname}' at {ch}:{resseq}: Missing {missing_atoms}")

            conn_info = _check_connectivity_integrity(group, ff_resname, atoms, conv)
            if conn_info["is_broken"]:
                warnings.append(f"Residue {resname} ({resseq}) at chain {ch} is BROKEN!")

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


def _check_missing_atoms(group: str, ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    """
    Kontroluje, zda jsou přítomny všechny atomy definované v sekci 'atom' pro dané reziduum.
    """
    if not group or not ff_name or group not in conv or ff_name not in conv[group]:
        return []

    # Získání seznamu povinných atomů ze slovníku
    res_def = conv[group][ff_name]
    required_atoms = set(res_def.get("atom", {}).keys())

    # Porovnání se seznamem atomů nalezených v PDB
    actual_atoms_set = set(atoms)
    missing = [a for a in required_atoms if a not in actual_atoms_set]

    return sorted(missing)


def _check_connectivity_integrity(group: str, ff_name: str, atoms: List[str], conv: Dict) -> Dict[str, Any]:
    """
    Analyzuje integritu rezidua na základě mapy konektivity ze slovníku.
    Prověřuje, zda jsou všechny přítomné atomy propojeny v jeden celek.
    Pokud je reziduum rozděleno na více kusů, označí ho jako 'broken'.
    """
    # Pokud reziduum neznáme nebo nemá definovanou konektivitu, považujeme ho za celistvé
    if not group or not ff_name or group not in conv or ff_name not in conv[group]:
        return {"is_broken": False, "components": []}

    res_def = conv[group][ff_name]
    conn_map = res_def.get("connectivity", {})

    # Pro ionty nebo jednoduché molekuly bez definovaných vazeb v JSONu nehlásíme chybu
    if not conn_map:
        return {"is_broken": False, "components": [atoms] if atoms else []}

    present_atoms = set(atoms)
    if not present_atoms:
        return {"is_broken": False, "components": []}

    # Sestavení neorientovaného grafu (adjacency list)
    # Connectivity v JSONu je často definována jednosměrně (rodič -> potomci),
    # my potřebujeme obousměrné vazby pro hledání souvislých komponent.
    graph = {atom: set() for atom in present_atoms}
    for u, neighbors in conn_map.items():
        if u in present_atoms:
            for v in neighbors:
                if v in present_atoms:
                    graph[u].add(v)
                    graph[v].add(u)

    # Algoritmus pro nalezení komponent souvislosti (BFS - Prohledávání do šířky)
    visited = set()
    components = []

    # Seřadíme atomy, aby byl výsledek vždy stejný (deterministický)
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

    # Reziduum je rozbité ("broken"), pokud se rozpadlo na více nepropojených částí
    is_broken = len(components) > 1

    return {
        "is_broken": is_broken,
        "components": components,
        "component_count": len(components)
    }