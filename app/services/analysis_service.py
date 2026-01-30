# app/services/analysis_service.py

from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from app.utils.alias import resn_alias

@lru_cache(maxsize=1)
def load_converting_dictionary() -> Dict:
    candidates = [
        os.path.join(os.getcwd(), "converting_dictionary.json"),
        os.path.join(os.getcwd(), "app", "resources", "converting_dictionary.json"),
        os.path.join(os.getcwd(), "app", "data", "converting_dictionary.json"),
    ]
    for path in candidates:
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    raise FileNotFoundError("converting_dictionary.json not found.")

def _parse_residues_from_pdb(pdb_text: str, chain: Optional[str]) -> List[Tuple[str, int, str, str]]:
    """
    Vrací unikátní residua v pořadí podle resseq: (chain, resseq, icode, resname)
    """
    seen = set()
    residues = []

    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        resname = line[17:20].strip()
        ch = (line[21] or "").strip() or "?"
        resseq_raw = line[22:26].strip()
        icode = (line[26] or " ").strip()  # insertion code
        if not resseq_raw:
            continue

        try:
            resseq = int(resseq_raw)
        except ValueError:
            continue

        if chain and ch != chain:
            continue

        key = (ch, resseq, icode)
        if key in seen:
            continue
        seen.add(key)
        residues.append((ch, resseq, icode, resname))

    residues.sort(key=lambda x: (x[0], x[1], x[2]))
    return residues

def _infer_group(resname: str, conv: Dict) -> Optional[str]:
    # converting_dictionary má top-level "D" a "R" :contentReference[oaicite:0]{index=0}
    for g in ("D", "R"):
        if g in conv and resname in conv[g]:
            return g
    if resname.startswith("D"):
        return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"):
        return "R"
    return None

def _pick_variant(group: Optional[str], pdb_resname: str, conv: Dict, terminal: str) -> Tuple[Optional[str], bool]:
    """
    terminal: "5" | "3" | ""
    """
    if not group or group not in conv:
        return None, False

    aliased = resn_alias(pdb_resname)  # např. A -> RA, A3 -> RA3 :contentReference[oaicite:1]{index=1}

    candidates: List[str] = []
    if terminal == "5":
        candidates += [f"{aliased}5", f"{pdb_resname}5"]
    elif terminal == "3":
        candidates += [f"{aliased}3", f"{pdb_resname}3"]

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

    # pokud chain nebyl specifikován, vezmeme první
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

    for ch, resseq, icode, resname in residues:
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
                })

        pos += 1
        group = _infer_group(resname, conv)
        terminal = "5" if resseq == first else ("3" if resseq == last else "")
        ff_resname, known = _pick_variant(group, resname, conv, terminal)

        if not known:
            warnings.append(f"Unknown residue '{resname}' at {chain}:{resseq}{icode or ''}")

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
        })

        prev_resseq = resseq

    return {"chain": chain, "tokens": tokens, "warnings": warnings}
