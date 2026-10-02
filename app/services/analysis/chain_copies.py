"""Identical copies of one molecule in a single PDB entry (e.g. 3SKR)."""

from __future__ import annotations

from typing import Any, Dict, List, Set, Tuple


# ---------------------------------------------------------------------------
# Identické kopie molekuly v jednom PDB (např. 3SKR: řetězce A a B jsou dvě
# kopie téhož riboswitche v asymetrické jednotce, REMARK 350 je uvádí jako
# dvě samostatné monomerní biologické jednotky). Inverzní operace k
# rozbalení symetrie - uživatel si může ponechat jen jednu kopii.
# ---------------------------------------------------------------------------

_NON_POLYMER_RESNAMES = {"HOH", "WAT", "SOL", "DOD"}


def _chain_sequences(pdb_text: str) -> Tuple[Dict[str, Tuple[str, ...]], bool]:
    """Sekvence řetězců: SEQRES, pokud je v souboru, jinak pozorovaná rezidua."""
    seqres: Dict[str, List[str]] = {}
    for line in pdb_text.splitlines():
        if line.startswith("SEQRES"):
            seqres.setdefault(line[11].strip() or "?", []).extend(line[19:].split())
    if seqres:
        return {chain: tuple(names) for chain, names in seqres.items()}, True

    observed: Dict[str, List[str]] = {}
    seen = set()
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        chain = line[21].strip() or "?"
        key = (chain, line[22:27])
        if key in seen:
            continue
        seen.add(key)
        observed.setdefault(chain, []).append(line[17:20].strip())
    return {chain: tuple(names) for chain, names in observed.items()}, False


def _biomolecule_chains(pdb_text: str) -> List[Set[str]]:
    """Řetězce jednotlivých BIOMOLECULE z REMARK 350 (APPLY THE FOLLOWING TO CHAINS)."""
    units: List[Set[str]] = []
    for line in pdb_text.splitlines():
        if not line.startswith("REMARK 350"):
            continue
        text = line[10:].strip()
        if text.startswith("BIOMOLECULE:"):
            units.append(set())
        elif units and ("APPLY THE FOLLOWING TO CHAINS:" in text or text.startswith("AND CHAINS:")):
            listed = text.split(":", 1)[1]
            units[-1].update(c.strip() for c in listed.split(",") if c.strip())
    return units


def find_identical_chain_copies(pdb_text: str) -> List[Dict[str, Any]]:
    """
    Skupiny řetězců se stejnou sekvencí. Ke každému řetězci statistiky
    polymeru (pozorovaná rezidua, atomy, průměrný B-faktor a occupancy) a
    počet navázaných heteroskupin (ionty/ligandy/vody se stejným chain ID).
    Doporučená kopie: nejúplnější, pak nejnižší průměrný B-faktor, pak
    nejvyšší occupancy.
    """
    sequences, from_seqres = _chain_sequences(pdb_text)
    groups: Dict[Tuple[str, ...], List[str]] = {}
    for chain, sequence in sequences.items():
        if len(sequence) >= 2:
            groups.setdefault(sequence, []).append(chain)
    groups = {seq: chains for seq, chains in groups.items() if len(chains) > 1}
    if not groups:
        return []

    stats: Dict[str, Dict[str, Any]] = {}
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip() or "?"
        if chain not in sequences:
            continue
        resname = line[17:20].strip()
        entry = stats.setdefault(chain, {"residues": set(), "atoms": 0, "b": 0.0, "occ": 0.0, "hetero": set()})
        is_polymer = line.startswith("ATOM") or (resname in sequences[chain] and resname not in _NON_POLYMER_RESNAMES)
        if not is_polymer:
            entry["hetero"].add(line[17:27])
            continue
        entry["residues"].add(line[22:27])
        entry["atoms"] += 1
        try:
            entry["b"] += float(line[60:66])
        except ValueError:
            pass
        try:
            entry["occ"] += float(line[54:60])
        except ValueError:
            entry["occ"] += 1.0

    units = _biomolecule_chains(pdb_text)
    result = []
    for sequence, chains in groups.items():
        rows = []
        for chain in chains:
            entry = stats.get(chain)
            if not entry or not entry["atoms"]:
                continue
            rows.append({
                "chain": chain,
                "observedResidues": len(entry["residues"]),
                "atoms": entry["atoms"],
                "bFactor": round(entry["b"] / entry["atoms"], 1),
                "occupancy": round(entry["occ"] / entry["atoms"] * 100, 1),
                "heteroGroups": len(entry["hetero"]),
            })
        if len(rows) < 2:
            continue
        best = min(rows, key=lambda r: (-r["observedResidues"], -r["atoms"], r["bFactor"], -r["occupancy"], r["chain"]))
        in_units = [next((i for i, unit in enumerate(units) if r["chain"] in unit), None) for r in rows]
        result.append({
            "chains": rows,
            "recommended": best["chain"],
            "sequenceLength": len(sequence),
            # REMARK 350 řadí každou kopii do jiné biologické jednotky.
            "separateBiologicalUnits": None not in in_units and len(set(in_units)) == len(rows),
            "fromSeqres": from_seqres,
        })
    return result
