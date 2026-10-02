"""Alternative locations (AltLocs), models and symmetry: what a structure
contains and which conformation to recommend."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple


from .chain_copies import find_identical_chain_copies
from .sequence import _BOND_DISTANCE_LIMIT_ANGSTROM

# ---------------------------------------------------------------------------
# AltLoc doporučení podle geometrie
# ---------------------------------------------------------------------------
# Písmena altlocu ve dvou sousedních reziduích na sebe nemusí navazovat (a
# stejné písmeno to taky nezaručuje) - rozhoduje skutečná vzdálenost
# vazebných atomů přes rozhraní reziduí (O3'-P u nukleových kyselin, C-N u
# proteinů, stejné meze jako _is_chemically_impossible_bond). Rezidua, jejichž
# altloc atomy leží přímo na takové vazbě, tvoří jeden "kus", o kterém se
# rozhoduje najednou: nejdřív co nejméně zlomů (i vůči sousedům bez altlocu),
# teprve potom occupancy a B-faktor sečtené přes celý kus.

_LINK_ATOMS_BY_GROUP = {"R": ("O3'", "P"), "P": ("C", "N")}

AltlocResidueKey = Tuple[str, int, str, str]  # (chain, resseq, icode, resname)
AltlocAtoms = Dict[AltlocResidueKey, Dict[str, Dict[str, Tuple[float, float, float]]]]


def _altloc_selection_key(key: AltlocResidueKey) -> str:
    # Stejný klíč, jaký posílá frontend a čte clean_pdb_altlocs.
    chain, resseq, _icode, resname = key
    return f"{chain}_{resseq}_{resname}"


def _parse_altloc_geometry(pdb_text: str) -> Tuple[Dict[str, List[AltlocResidueKey]], AltlocAtoms]:
    """
    Souřadnice atomů prvního modelu po reziduích: {key: {atom: {altloc:
    coord}}}, altloc " " = atom bez alternativ. Vrací i pořadí reziduí v
    každém řetězci tak, jak jdou v souboru.
    """
    order: Dict[str, List[AltlocResidueKey]] = {}
    atoms: AltlocAtoms = {}
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            resseq = int(line[22:26])
            coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        chain = line[21].strip() or "?"
        key = (chain, resseq, line[26].strip(), line[17:20].strip())
        if key not in atoms:
            atoms[key] = {}
            order.setdefault(chain, []).append(key)
        atom_name = line[12:16].strip().replace("*", "'")
        atoms[key].setdefault(atom_name, {}).setdefault(line[16], coord)
    return order, atoms


def _atom_coord(variants: Dict[str, Tuple[float, float, float]], letter: Optional[str]):
    if letter is not None and letter in variants:
        return variants[letter]
    if " " in variants:
        return variants[" "]
    return variants[sorted(variants)[0]]


def _residue_links(order, atoms: AltlocAtoms) -> List[Tuple[AltlocResidueKey, AltlocResidueKey, str, str, float]]:
    """(prev, curr, atom na prev, atom na curr, mez) pro sousedy v řetězci, které tvoří polymer."""
    links = []
    for keys in order.values():
        for prev, curr in zip(keys, keys[1:]):
            for group, (prev_atom, curr_atom) in _LINK_ATOMS_BY_GROUP.items():
                if prev_atom in atoms[prev] and curr_atom in atoms[curr]:
                    links.append((prev, curr, prev_atom, curr_atom, _BOND_DISTANCE_LIMIT_ANGSTROM[group]))
                    break
    return links


def _link_is_broken(atoms: AltlocAtoms, link, prev_letter: Optional[str], curr_letter: Optional[str]) -> bool:
    prev, curr, prev_atom, curr_atom, limit = link
    p = _atom_coord(atoms[prev][prev_atom], prev_letter)
    c = _atom_coord(atoms[curr][curr_atom], curr_letter)
    return sum((a - b) ** 2 for a, b in zip(p, c)) ** 0.5 > limit


def _link_letters(atoms: AltlocAtoms, key: AltlocResidueKey, atom_name: str) -> List[str]:
    return sorted(letter for letter in atoms[key][atom_name] if letter != " ")


def _recommend_altlocs(
    residue_stats: Dict[AltlocResidueKey, Dict[str, Dict[str, float]]],
    order,
    atoms: AltlocAtoms,
) -> Dict[AltlocResidueKey, str]:
    """Doporučené písmeno pro každé altloc reziduum (viz komentář nad sekcí)."""
    links = _residue_links(order, atoms)
    links_by_residue: Dict[AltlocResidueKey, List] = {}
    for link in links:
        links_by_residue.setdefault(link[0], []).append(link)
        links_by_residue.setdefault(link[1], []).append(link)

    # Kusy: sousední altloc rezidua spojená vazbou, na které má altloc aspoň
    # jeden z obou vazebných atomů.
    parent = {key: key for key in residue_stats}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    for prev, curr, prev_atom, curr_atom, _limit in links:
        if prev in parent and curr in parent and (
            _link_letters(atoms, prev, prev_atom) or _link_letters(atoms, curr, curr_atom)
        ):
            parent[find(prev)] = find(curr)

    chain_position = {key: index for keys in order.values() for index, key in enumerate(keys)}
    pieces: Dict[AltlocResidueKey, List[AltlocResidueKey]] = {}
    for key in residue_stats:
        pieces.setdefault(find(key), []).append(key)

    recommended: Dict[AltlocResidueKey, str] = {}
    for members in pieces.values():
        members.sort(key=lambda k: (k[0], chain_position.get(k, 0)))
        member_set = set(members)

        def node_cost(key, letter):
            stats = residue_stats[key][letter]
            breaks = 0
            for link in links_by_residue.get(key, ()):
                other = link[1] if link[0] == key else link[0]
                if other in member_set:
                    continue
                letters = (letter, None) if link[0] == key else (None, letter)
                breaks += _link_is_broken(atoms, link, *letters)
            return (breaks, -stats["occupancy"], stats["bFactor"])

        def edge_breaks(prev, prev_letter, curr, curr_letter):
            return sum(
                _link_is_broken(atoms, link, prev_letter, curr_letter)
                for link in links_by_residue.get(prev, ())
                if link[0] == prev and link[1] == curr
            )

        # Viterbi po řetězci: lexikografický součet (zlomy, -occupancy, B).
        best: Dict[str, Tuple[Tuple[float, float, float], List[str]]] = {}
        for index, key in enumerate(members):
            new_best = {}
            for letter in sorted(residue_stats[key]):
                own = node_cost(key, letter)
                if index == 0:
                    new_best[letter] = (own, [letter])
                    continue
                prev_key = members[index - 1]
                candidates = []
                for prev_letter, (cost, path) in best.items():
                    extra = edge_breaks(prev_key, prev_letter, key, letter)
                    total = (cost[0] + own[0] + extra, cost[1] + own[1], cost[2] + own[2])
                    candidates.append((total, path + [letter]))
                new_best[letter] = min(candidates, key=lambda item: (item[0], item[1]))
            best = new_best
        _cost, path = min(best.values(), key=lambda item: (item[0], item[1]))
        recommended.update(zip(members, path))
    return recommended


def find_altloc_selection_breaks(pdb_text: str, selection: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Zlomy řetězce, které vytvoří daný výběr altloců a kterým by šlo jinou
    volbou u dotčených reziduí předejít. Zlomy, které vzniknou při jakékoli
    volbě (skutečná mezera ve struktuře), se nehlásí.
    """
    order, atoms = _parse_altloc_geometry(pdb_text)

    def chosen(key):
        return selection.get(_altloc_selection_key(key))

    def options(key, atom_name):
        return _link_letters(atoms, key, atom_name) or [None]

    breaks = []
    for link in _residue_links(order, atoms):
        prev, curr, prev_atom, curr_atom, _limit = link
        if not (_link_letters(atoms, prev, prev_atom) or _link_letters(atoms, curr, curr_atom)):
            continue
        if not _link_is_broken(atoms, link, chosen(prev), chosen(curr)):
            continue
        avoidable = any(
            not _link_is_broken(atoms, link, p, c)
            for p in options(prev, prev_atom)
            for c in options(curr, curr_atom)
        )
        if avoidable:
            breaks.append({
                "chain": prev[0],
                "prev": {"resseq": prev[1], "icode": prev[2], "resname": prev[3], "altloc": chosen(prev)},
                "next": {"resseq": curr[1], "icode": curr[2], "resname": curr[3], "altloc": chosen(curr)},
            })
    return breaks


def analyze_pdb_altlocs(pdb_text: str) -> Dict[str, Any]:
    """
    PROJDE PDB SOUBOR A IDENTIFIKUJE VŠECHNY ALTERNATIVNÍ POZICE (ALTLOCS),
    JEJICH OBSAZENOST A B-FAKTOR (PRŮMĚR PŘES ATOMY DANÉ VARIANTY). VRACÍ
    STRUKTUROVANÝ DICT (JSON) PRO FRONTEND. DOPORUČENÍ (recommended_alt) SE
    ŘÍDÍ NEJDŘÍV GEOMETRIÍ - VIZ _recommend_altlocs.

    NOVĚ: DETEKUJE PŘÍTOMNOST VÍCE MODELŮ A SYMETRIE (REMARK 350).
    """
    models = []
    has_symmetry = False
    # key -> altloc -> [součet occupancy, součet B, počet atomů]
    sums: Dict[AltlocResidueKey, Dict[str, List[float]]] = {}

    for line in pdb_text.splitlines():
        if line.startswith("MODEL "):
            try:
                model_num = int(line[6:].strip())
                if model_num not in models:
                    models.append(model_num)
            except ValueError:
                pass
            continue

        # Biological Assembly: REMARK 350 s transformační maticí BIOMT
        if line.startswith("REMARK 350") and "BIOMT" in line:
            has_symmetry = True
            continue

        if line.startswith(("ATOM", "HETATM")) and line[16] != " ":
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            try:
                occupancy = float(line[54:60].strip())
            except ValueError:
                occupancy = 1.0
            try:
                b_factor = float(line[60:66].strip())
            except ValueError:
                b_factor = 0.0
            key = (line[21].strip() or "?", resseq, line[26].strip(), line[17:20].strip())
            acc = sums.setdefault(key, {}).setdefault(line[16], [0.0, 0.0, 0])
            acc[0] += occupancy
            acc[1] += b_factor
            acc[2] += 1

    residue_stats = {
        key: {
            letter: {"occupancy": occ / count, "bFactor": b / count}
            for letter, (occ, b, count) in letters.items()
        }
        for key, letters in sums.items()
    }

    order, atoms = _parse_altloc_geometry(pdb_text)
    # Rezidua jen z dalších modelů geometrii prvního modelu nemají - rozhodne occupancy/B.
    for key in residue_stats:
        if key not in atoms:
            atoms[key] = {}
            order.setdefault(key[0], []).append(key)
    recommended = _recommend_altlocs(residue_stats, order, atoms)

    result_residues = []
    for key in sorted(residue_stats, key=lambda k: (k[0], k[1], k[2])):
        chain, resseq, _icode, resname = key
        result_residues.append({
            "chain": chain,
            "resseq": resseq,
            "resname": resname,
            "altLocs": {
                letter: {"occupancy": round(stats["occupancy"] * 100, 1), "bFactor": round(stats["bFactor"], 2)}
                for letter, stats in residue_stats[key].items()
            },
            "recommended_alt": recommended[key],
        })

    return {
        "models": models,  # Pole s čísly modelů (např. [1, 2, 3])
        "hasSymmetry": has_symmetry,  # True/False, pokud existuje BIOMT matice
        "hasAltLocs": len(result_residues) > 0,
        "residues": result_residues,
        "copyGroups": find_identical_chain_copies(pdb_text),
    }
