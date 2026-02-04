# app/services/analysis_service.py

from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Set

from app.utils.alias import resn_alias


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
    Vrací unikátní residua s jejich seznamem atomů: (chain, resseq, icode, resname, atoms)
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
        atom_name = line[12:16].strip()

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
    """Zjišťuje kategorii rezidua (R, D, P, W3 atd.) přímo ze slovníku."""
    for category in conv.keys():
        if resname in conv[category]:
            return category

    # Fallback pro nukleové kyseliny
    if resname.startswith("D"): return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"): return "R"
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
    # Mapování konců (HPC označení)
    if terminal == "5":
        candidates += [f"{aliased}5", f"{pdb_resname}5", f"N{aliased}", f"N{pdb_resname}"]
    elif terminal == "3":
        candidates += [f"{aliased}3", f"{pdb_resname}3", f"C{aliased}", f"C{pdb_resname}"]

    candidates += [aliased, pdb_resname]

    for k in candidates:
        if k in conv[group]:
            return k, True
    return candidates[0] if candidates else None, False


def build_sequence_tokens(pdb_text: str, chain: Optional[str] = None, fill_gaps: bool = True):
    conv = load_converting_dictionary()
    residues = _parse_residues_from_pdb(pdb_text, chain)

    if not residues:
        return {"chain": chain, "tokens": [], "warnings": ["No residues found in PDB."]}

    if chain is None:
        chain = residues[0][0]
        residues = [r for r in residues if r[0] == chain]

    residues.sort(key=lambda x: (x[1], x[2]))
    first = residues[0][1]
    last = residues[-1][1]

    tokens = []
    warnings = []
    pos = 0
    prev_resseq = None

    for ch, resseq, icode, resname, atoms in residues:
        # Logika pro vyplnění mezer
        if fill_gaps and prev_resseq is not None and resseq > prev_resseq + 1:
            for missing in range(prev_resseq + 1, resseq):
                pos += 1
                tokens.append({
                    "position": pos,
                    "chain": chain,
                    "resseq": None,
                    "icode": None,
                    "pdb_resname": "0",
                    "is_gap": True,
                    "group": None,
                    "ff_resname": None,
                    "known": False,
                    "atoms": [],
                    "missing_atoms": []
                })

        pos += 1
        group = _infer_group(resname, conv)
        # Určení konce (5' nebo 3') pro DNA/RNA i Proteiny
        terminal = "5" if resseq == first else ("3" if resseq == last else "")
        ff_resname, known = _pick_variant(group, resname, conv, terminal)

        # NOVINKA: Kontrola chybějících atomů podle konektivity
        missing_atoms = []
        if known:
            missing_atoms = _validate_connectivity(group, ff_resname, atoms, conv)

        if not known:
            warnings.append(f"Unknown residue '{resname}' at {chain}:{resseq}{icode or ''}")
        elif missing_atoms:
            warnings.append(f"Incomplete residue '{resname}' at {chain}:{resseq}: Missing {len(missing_atoms)} atoms.")

        tokens.append({
            "position": pos,
            "chain": chain,
            "resseq": resseq,
            "icode": icode or "",
            "pdb_resname": resname,
            "is_gap": False,
            "group": group,
            "ff_resname": ff_resname,
            "known": known,
            "atoms": atoms,
            "missing_atoms": missing_atoms  # Seznam atomů, které chybí pro úplnou konektivitu
        })

        prev_resseq = resseq

    return {"chain": chain, "tokens": tokens, "warnings": warnings}