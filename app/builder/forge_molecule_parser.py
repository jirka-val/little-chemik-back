#!/usr/bin/env python3
"""
First-pass parser for FORGE4SIMS-style structure JSON files.

Purpose
-------
Create a lightweight internal Molecule/Chain/Residue/Atom representation from:
  1) a structure JSON with keys `pdb_text` and `missing_atoms`,
  2) converting_dictionary.json.

Design contract
---------------
- PDB text is used as the coordinate source only.
- The `missing_atoms` JSON section plus converting_dictionary are the authority for
  supported residue identity, expected atom list, and atom order.
- Builder-supported residues are currently polymer mol_types R, D, P by default.
- Waters, ions, unsupported ligands and other non-builder residues are preserved as
  HETATM passthrough records, but are not standardized and not built.
- No atoms are built here; missing supported atoms have coord=None.

This module intentionally does not depend on MDAnalysis/RDKit/Biopython.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

Coord = Tuple[float, float, float]
ResidueKey = Tuple[str, int, str]
AtomKey = Tuple[str, int, str, str]


# -----------------------------------------------------------------------------
# Data model
# -----------------------------------------------------------------------------

@dataclass
class Atom:
    name: str
    element: Optional[str]
    coord: Optional[Coord]
    serial: Optional[int] = None
    built: bool = False
    build_source: Optional[str] = None  # "input", "template_rule", "user_torsion", ...
    build_rule_index: Optional[int] = None
    original_name: Optional[str] = None
    altloc: str = ""
    occupancy: Optional[float] = None
    bfactor: Optional[float] = None


@dataclass
class Residue:
    chain_id: str
    resseq: int
    icode: str
    ff_resname: str
    atoms: Dict[str, Atom]
    index_in_chain: int
    original_resname: Optional[str] = None
    protonation_source: Optional[str] = None
    state_assignment_source: Optional[str] = None
    group: Optional[str] = None
    is_gap: bool = False
    is_broken: bool = False
    # "gap" / "ter" / "geometry" / "chain_end" - proč byl tenhle konec rezidua
    # (5'/3', N-/C-terminus) vůbec vytvořen (viz analysis_service.py
    # terminus_reason). "chain_end" = reálný konec řetězce; ostatní hodnoty =
    # umělá terminalita na okraji mezery, kde je potřeba mezi tohle reziduum
    # a další v témže chain_id vypsat TER (viz forge_service.molecule_to_pdb).
    terminus_reason: Optional[str] = None
    connectivity_parts: List[List[str]] = field(default_factory=list)
    # Coordinate records present in the input but not expected by the residue's
    # current converting-dictionary state.  State-assignment layers may need to
    # restore them after changing residue identity (for example CYX -> CYS/HG).
    observed_extra_atoms: Dict[str, Atom] = field(default_factory=dict)


@dataclass
class Chain:
    chain_id: str
    residues: List[Residue] = field(default_factory=list)


@dataclass
class PDBAtomRecord:
    record_name: str
    serial: Optional[int]
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    resseq: int
    icode: str
    coord: Coord
    occupancy: Optional[float] = None
    bfactor: Optional[float] = None
    element: Optional[str] = None
    charge: str = ""
    # Filled from the corresponding FORGE token for passthrough waters/ions.
    # Raw PDB records do not carry this force-field identity themselves.
    group: Optional[str] = None
    ff_resname: Optional[str] = None


@dataclass
class PeriodicBox:
    """Periodic unit-cell metadata in Angstrom and Cartesian coordinates.

    Vectors are stored as three row vectors in GROMACS-compatible reduced
    orientation (a along x, b in xy).  `origin` identifies the origin of the
    primary primitive cell; solvated FORGE systems use the zero vector.
    """

    vectors: Tuple[Coord, Coord, Coord]
    origin: Coord = (0.0, 0.0, 0.0)
    shape: str = "orthorhombic"


@dataclass
class Molecule:
    chains: Dict[str, Chain]
    # Explicit non-polymer covalent links.  AtomKey is
    # (chain_id, resseq, icode, atom_name); endpoints are stored canonically.
    covalent_bonds: List[Tuple[AtomKey, AtomKey]] = field(default_factory=list)
    passthrough_atoms: List[PDBAtomRecord] = field(default_factory=list)
    unassigned_records: List[PDBAtomRecord] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    periodic_box: Optional[PeriodicBox] = None


# -----------------------------------------------------------------------------
# PDB parsing and formatting
# -----------------------------------------------------------------------------

def _safe_int(s: str) -> Optional[int]:
    s = s.strip()
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        return None


def _safe_float(s: str) -> Optional[float]:
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def infer_element(atom_name: str, pdb_element: Optional[str] = None) -> Optional[str]:
    """Infer element from a PDB atom name with a conservative fallback."""
    if pdb_element and pdb_element.strip():
        return pdb_element.strip().upper()
    name = atom_name.strip()
    if not name:
        return None
    # Remove leading digits, common for H names; keep first alphabetic character.
    letters = "".join(ch for ch in name if ch.isalpha())
    if not letters:
        return None
    # PDB biomolecular atom names in our supported residues are mostly one-letter elements.
    # Two-letter elements can be handled later for ions/ligands; passthrough keeps PDB element.
    return letters[0].upper()


def parse_pdb_atoms(pdb_text: str) -> List[PDBAtomRecord]:
    records: List[PDBAtomRecord] = []
    for line_number, line in enumerate(pdb_text.splitlines(), 1):
        rec = line[0:6].strip()
        if rec not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            raise ValueError(
                f"Malformed {rec} record at PDB line {line_number}: "
                f"expected at least 54 columns, found {len(line)}"
            )
        try:
            serial = _safe_int(line[6:11])
            atom_name = line[12:16].strip()
            altloc = line[16:17].strip()
            resname = line[17:20].strip()
            chain_id = line[21:22].strip()
            resseq = _safe_int(line[22:26])
            icode = line[26:27].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occupancy = _safe_float(line[54:60])
            bfactor = _safe_float(line[60:66])
            if line[54:60].strip() and occupancy is None:
                raise ValueError("invalid occupancy")
            if line[60:66].strip() and bfactor is None:
                raise ValueError("invalid B-factor")
            if not all(math.isfinite(value) for value in (x, y, z)):
                raise ValueError("non-finite coordinate")
            element = line[76:78].strip() if len(line) >= 78 else ""
            charge = line[78:80].strip() if len(line) >= 80 else ""
        except ValueError as exc:
            raise ValueError(
                f"Malformed {rec} record at PDB line {line_number}: "
                f"invalid numeric field"
            ) from exc
        missing_fields = [
            name
            for name, value in (
                ("serial", serial),
                ("atom_name", atom_name),
                ("resname", resname),
                ("resseq", resseq),
            )
            if value in (None, "")
        ]
        if missing_fields:
            raise ValueError(
                f"Malformed {rec} record at PDB line {line_number}: "
                f"missing required field(s): {', '.join(missing_fields)}"
            )
        records.append(PDBAtomRecord(
            record_name=rec,
            serial=serial,
            atom_name=atom_name,
            altloc=altloc,
            resname=resname,
            chain_id=chain_id,
            resseq=resseq,
            icode=icode,
            coord=(x, y, z),
            occupancy=occupancy,
            bfactor=bfactor,
            element=infer_element(atom_name, element),
            charge=charge,
        ))
    return records


def group_records_by_atom(records: Iterable[PDBAtomRecord]) -> Dict[AtomKey, List[PDBAtomRecord]]:
    out: Dict[AtomKey, List[PDBAtomRecord]] = {}
    for r in records:
        key = (r.chain_id, r.resseq, r.icode, r.atom_name)
        out.setdefault(key, []).append(r)
    return out


def group_records_by_residue(records: Iterable[PDBAtomRecord]) -> Dict[ResidueKey, List[PDBAtomRecord]]:
    out: Dict[ResidueKey, List[PDBAtomRecord]] = {}
    for r in records:
        key = (r.chain_id, r.resseq, r.icode)
        out.setdefault(key, []).append(r)
    return out


def choose_unique_record(records: Sequence[PDBAtomRecord], context: str = "") -> Optional[PDBAtomRecord]:
    """Return a unique coordinate record. Duplicate records indicate uncleaned input."""
    if not records:
        return None
    if len(records) == 1:
        return records[0]
    details = ", ".join(f"serial={r.serial}, altloc={r.altloc or '<blank>'}" for r in records)
    raise ValueError(
        f"Duplicate coordinate records for {context}. "
        "structure_json must be cleaned upstream before builder parsing. "
        f"Records: {details}"
    )


def format_atom_name(atom_name: str, element: Optional[str]) -> str:
    """Return a four-character PDB-ish atom name field."""
    name = atom_name[:4]
    elem = (element or infer_element(atom_name) or "").strip()
    # For one-letter elements, conventional PDB atom names are right-justified in columns 13-16
    # unless four-character names are needed. This is sufficient for simple ATOM output.
    if len(name) < 4 and len(elem) == 1 and not name[0].isdigit():
        return f" {name:<3s}"[:4]
    return f"{name:<4s}"[:4]


def format_pdb_atom_line(
    serial: int,
    record_name: str,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    icode: str,
    coord: Coord,
    occupancy: float = 1.00,
    bfactor: float = 0.00,
    element: Optional[str] = None,
    altloc: str = "",
) -> str:
    element = infer_element(atom_name, element) or ""
    atom_field = format_atom_name(atom_name, element)
    x, y, z = coord
    return (
        f"{record_name:<6s}{serial:5d} {atom_field}{altloc[:1]:1s}"
        f"{resname:>3s} {chain_id[:1]:1s}{resseq:4d}{icode[:1]:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element:>2s}"
    )


# -----------------------------------------------------------------------------
# Converting dictionary helpers
# -----------------------------------------------------------------------------

def load_json(path: str | Path) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def residue_exists_in_converting(converting: Mapping[str, Any], group: Optional[str], ff_resname: str) -> bool:
    return bool(group and group in converting and ff_resname in converting[group])


def expected_atom_order(converting: Mapping[str, Any], group: str, ff_resname: str) -> List[str]:
    return list(converting[group][ff_resname]["atom"].keys())


# -----------------------------------------------------------------------------
# Molecule construction
# -----------------------------------------------------------------------------

def build_molecule_from_forge_json(
    structure_json: Mapping[str, Any],
    converting: Mapping[str, Any],
    supported_mol_types: Sequence[str] = ("R", "D", "P"),
) -> Molecule:
    """Build Molecule object from a FORGE-style structure JSON.

    `supported_mol_types` controls which converting-dictionary mol_types become
    standardized builder residues. Other tokens are preserved as HETATM passthrough.
    """
    pdb_text = structure_json.get("pdb_text", "")
    records = parse_pdb_atoms(pdb_text)
    by_atom = group_records_by_atom(records)
    by_residue = group_records_by_residue(records)

    chains_json = structure_json.get("missing_atoms", {}).get("chains", {})
    mol = Molecule(chains={})
    consumed_ids: set[int] = set()


    for chain_id, chain_data in chains_json.items():
        chain = Chain(chain_id=chain_id)
        mol.chains[chain_id] = chain
        tokens = list(chain_data.get("tokens", []))
        for idx, token in enumerate(tokens):
            if token.get("is_gap"):
                continue
            group = token.get("group")
            ff_resname = token.get("ff_resname") or token.get("pdb_resname")
            pdb_resname = token.get("pdb_resname") or ff_resname
            res_chain = token.get("chain", chain_id) or chain_id
            resseq = int(token.get("resseq"))
            icode = token.get("icode") or ""
            is_supported = (
                group in supported_mol_types
                and residue_exists_in_converting(converting, group, ff_resname)
            )

            if is_supported:
                atom_order = expected_atom_order(converting, group, ff_resname)
                atoms: Dict[str, Atom] = {}
                for atom_name in atom_order:
                    recs = by_atom.get((res_chain, resseq, icode, atom_name), [])
                    rec = choose_unique_record(recs, context=f"{res_chain}:{resseq}{icode} {ff_resname} {atom_name}")
                    if rec is not None:
                        consumed_ids.add(id(rec))
                        atoms[atom_name] = Atom(
                            name=atom_name,
                            element=infer_element(atom_name, rec.element),
                            coord=rec.coord,
                            serial=rec.serial,
                            built=False,
                            build_source="input",
                            original_name=rec.atom_name,
                            altloc=rec.altloc,
                            occupancy=rec.occupancy,
                            bfactor=rec.bfactor,
                        )
                    else:
                        atoms[atom_name] = Atom(
                            name=atom_name,
                            element=infer_element(atom_name),
                            coord=None,
                            serial=None,
                            built=False,
                            build_source=None,
                        )

                observed_extra_atoms: Dict[str, Atom] = {}
                expected_names = set(atom_order)
                for rec in by_residue.get((res_chain, resseq, icode), []):
                    if rec.atom_name in expected_names:
                        continue
                    duplicate_records = [
                        candidate
                        for candidate in by_atom.get(
                            (res_chain, resseq, icode, rec.atom_name), []
                        )
                        if candidate.resname == rec.resname
                    ]
                    selected = choose_unique_record(
                        duplicate_records,
                        context=(
                            f"{res_chain}:{resseq}{icode} {ff_resname} "
                            f"extra atom {rec.atom_name}"
                        ),
                    )
                    if selected is None or rec is not selected:
                        continue
                    consumed_ids.add(id(selected))
                    observed_extra_atoms[rec.atom_name] = Atom(
                        name=rec.atom_name,
                        element=infer_element(rec.atom_name, rec.element),
                        coord=rec.coord,
                        serial=rec.serial,
                        built=False,
                        build_source="input",
                        original_name=rec.atom_name,
                        altloc=rec.altloc,
                        occupancy=rec.occupancy,
                        bfactor=rec.bfactor,
                    )

                json_atoms = set(token.get("atoms", []))
                expected = set(atom_order)
                extra_json_atoms = sorted(json_atoms - expected)
                if extra_json_atoms:
                    mol.warnings.append(
                        f"{res_chain}:{resseq}{icode} {ff_resname}: JSON token lists atoms not in converting: {extra_json_atoms}"
                    )
                chain.residues.append(Residue(
                    chain_id=res_chain,
                    resseq=resseq,
                    icode=icode,
                    ff_resname=ff_resname,
                    atoms=atoms,
                    index_in_chain=len(chain.residues),
                    original_resname=pdb_resname,
                    group=group,
                    is_gap=bool(token.get("is_gap", False)),
                    is_broken=bool(token.get("is_broken", False)),
                    terminus_reason=token.get("terminus_reason"),
                    connectivity_parts=list(token.get("connectivity_parts", [])),
                    observed_extra_atoms=observed_extra_atoms,
                ))
            else:
                # Unsupported/non-builder token: preserve actual selected PDB records as HETATM.
                res_recs = [
                    r for r in by_residue.get((res_chain, resseq, icode), [])
                    if (not pdb_resname or r.resname == pdb_resname)
                ]
                atom_to_recs: Dict[str, List[PDBAtomRecord]] = {}
                for r in res_recs:
                    atom_to_recs.setdefault(r.atom_name, []).append(r)
                for atom_name, recs in atom_to_recs.items():
                    rec = choose_unique_record(recs, context=f"{res_chain}:{resseq}{icode} {ff_resname} {atom_name}")
                    if rec is not None:
                        rec.group = group
                        rec.ff_resname = ff_resname
                        consumed_ids.add(id(rec))
                        mol.passthrough_atoms.append(rec)

    # Preserve any unconsumed HETATM records not represented by tokens.
    for rec in records:
        if id(rec) in consumed_ids:
            continue
        if rec.record_name == "HETATM":
            mol.unassigned_records.append(rec)
            mol.passthrough_atoms.append(rec)

    return mol


# -----------------------------------------------------------------------------
# Reporting and output
# -----------------------------------------------------------------------------

def molecule_summary(mol: Molecule) -> Dict[str, Any]:
    n_chains = len(mol.chains)
    n_res = sum(len(ch.residues) for ch in mol.chains.values())
    n_atoms = 0
    n_known = 0
    n_missing = 0
    per_chain = {}
    for cid, chain in mol.chains.items():
        cres = len(chain.residues)
        catoms = sum(len(r.atoms) for r in chain.residues)
        cknown = sum(1 for r in chain.residues for a in r.atoms.values() if a.coord is not None)
        cmissing = catoms - cknown
        per_chain[cid] = {"residues": cres, "atoms": catoms, "known": cknown, "missing": cmissing}
        n_atoms += catoms
        n_known += cknown
        n_missing += cmissing
    return {
        "chains": n_chains,
        "supported_residues": n_res,
        "supported_atoms_expected": n_atoms,
        "supported_atoms_known": n_known,
        "supported_atoms_missing": n_missing,
        "passthrough_atoms": len(mol.passthrough_atoms),
        "unassigned_records": len(mol.unassigned_records),
        "warnings": len(mol.warnings),
        "per_chain": per_chain,
    }


def write_molecule_report(mol: Molecule, path: str | Path) -> None:
    s = molecule_summary(mol)
    lines: List[str] = []
    lines.append("Molecule parser report")
    lines.append("======================")
    lines.append(f"Chains: {s['chains']}")
    lines.append(f"Supported residues: {s['supported_residues']}")
    lines.append(f"Supported expected atoms: {s['supported_atoms_expected']}")
    lines.append(f"Supported atoms with coordinates: {s['supported_atoms_known']}")
    lines.append(f"Supported missing atoms: {s['supported_atoms_missing']}")
    lines.append(f"Passthrough HETATM atoms: {s['passthrough_atoms']}")
    lines.append(f"Unassigned passthrough records: {s['unassigned_records']}")
    lines.append("")
    lines.append("Per chain:")
    for cid, c in s["per_chain"].items():
        lines.append(
            f"  {cid or '<blank>'}: residues={c['residues']} atoms={c['atoms']} known={c['known']} missing={c['missing']}"
        )
    lines.append("")
    lines.append("Supported residues:")
    for cid, chain in mol.chains.items():
        for res in chain.residues:
            known = sum(1 for a in res.atoms.values() if a.coord is not None)
            missing = [name for name, a in res.atoms.items() if a.coord is None]
            lines.append(
                f"  {res.chain_id}:{res.resseq}{res.icode or ''} {res.original_resname}->{res.ff_resname} "
                f"atoms={len(res.atoms)} known={known} missing={len(missing)}"
            )
            if missing:
                lines.append(f"    missing: {', '.join(missing)}")
    if mol.passthrough_atoms:
        lines.append("")
        lines.append("Passthrough residues/atoms:")
        # Compact by residue.
        grouped: Dict[Tuple[str, int, str, str], List[str]] = {}
        for r in mol.passthrough_atoms:
            grouped.setdefault((r.chain_id, r.resseq, r.icode, r.resname), []).append(r.atom_name)
        for (ch, seq, ic, rn), atoms in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1], x[0][2], x[0][3])):
            lines.append(f"  {ch}:{seq}{ic or ''} {rn}: {len(atoms)} atoms ({', '.join(atoms[:12])}{'...' if len(atoms)>12 else ''})")
    if mol.warnings:
        lines.append("")
        lines.append("Warnings:")
        for w in mol.warnings:
            lines.append(f"  - {w}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
