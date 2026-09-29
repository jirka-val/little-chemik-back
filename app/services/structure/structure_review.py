"""
Kontroly pro krok 3.5 "Structure Check" (Expert) nad vstupním PDB - věci,
o kterých builder sám nerozhoduje, ale výsledek simulace ovlivní:

- atomy s nulovou obsazeností (krystalograf je domodeloval, v datech nejsou)
  -> ponechat, nebo nechat builder postavit postranní řetězec znovu;
- ASN/GLN amid otočený o 180° (rentgen nerozliší O od N) -> zjednodušené
  skóre vodíkových vazeb obou orientací, jako Reduce/MolProbity.

Rozhodnutí uživatele se aplikují na PDB text před builderem
(apply_structure_edits). Protonace HIS a odstraněné heterogeny jsou ve
forge_service.py - tam žijí data, ze kterých vycházejí.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

ResidueKey = Tuple[str, int, str]

# Atomy páteře - ty builder znovu postavit neumí (kotví celý řetězec).
_PROTEIN_BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})
_NUCLEIC_BACKBONE = frozenset({
    "P", "OP1", "OP2", "OP3", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
})
_NUCLEIC_RESNAMES = frozenset({"A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI"})

# Amidová skupina: (kyslík, dusík, vodíky na dusíku).
_AMIDES = {
    "ASN": ("OD1", "ND2", ("HD21", "HD22", "1HD2", "2HD2")),
    "GLN": ("OE1", "NE2", ("HE21", "HE22", "1HE2", "2HE2")),
}

_WATER = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3"})
_METAL_ELEMENTS = frozenset({
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CO", "NI", "CU", "CD", "SR", "BA", "CS", "RB", "LI",
})

# Role polárních atomů jako partnerů vodíkové vazby.
_DONOR_ONLY = {
    "LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "TRP": {"NE1"},
    "ASN": {"ND2"}, "GLN": {"NE2"},
}
_ACCEPTOR_ONLY = {
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"}, "GLN": {"OE1"},
}
_EITHER = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "HIS": {"ND1", "NE2"},
}

_HBOND_MAX = 3.5    # donor-akceptor, Å
_CLASH_MAX = 3.2    # dva donory / dva akceptory proti sobě, Å
_METAL_MAX = 3.0


@dataclass(frozen=True)
class _Atom:
    record: str
    chain: str
    resseq: int
    icode: str
    resname: str
    name: str
    element: str
    occupancy: Optional[float]
    xyz: Tuple[float, float, float]

    @property
    def key(self) -> ResidueKey:
        return (self.chain, self.resseq, self.icode)

    @property
    def is_hydrogen(self) -> bool:
        return self.element in ("H", "D") or (not self.element and self.name.lstrip("0123456789").startswith("H"))


def _parse_atoms(pdb_text: str) -> List[_Atom]:
    atoms = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            resseq = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        try:
            occupancy: Optional[float] = float(line[54:60])
        except ValueError:
            occupancy = None
        atoms.append(_Atom(
            record=line[:6].strip(),
            chain=line[21].strip() or "?",
            resseq=resseq,
            icode=line[26].strip(),
            resname=line[17:20].strip(),
            name=line[12:16].strip(),
            element=line[76:78].strip().upper() if len(line) >= 78 else "",
            occupancy=occupancy,
            xyz=xyz,
        ))
    return atoms


def _dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def _label(key: ResidueKey) -> str:
    chain, resseq, icode = key
    return f"{chain}{resseq}{icode}"


def _ref(key: ResidueKey) -> Dict[str, Any]:
    chain, resseq, icode = key
    return {"residue": _label(key), "chain": chain, "resseq": resseq, "icode": icode}


# ---------------------------------------------------------------------------
# Nulová obsazenost
# ---------------------------------------------------------------------------

def find_zero_occupancy(pdb_text: str, decided: Iterable[ResidueKey] = ()) -> List[Dict[str, Any]]:
    """Polymerní rezidua s těžkými atomy o obsazenosti 0 (bez těch, o kterých už uživatel rozhodl)."""
    decided = set(decided)
    by_residue: Dict[ResidueKey, List[_Atom]] = {}
    for atom in _parse_atoms(pdb_text):
        if atom.record != "ATOM" or atom.is_hydrogen or atom.occupancy is None:
            continue
        if atom.occupancy == 0.0 and atom.key not in decided:
            by_residue.setdefault(atom.key, []).append(atom)

    out = []
    for key, atoms in by_residue.items():
        resname = atoms[0].resname
        backbone = _NUCLEIC_BACKBONE if resname in _NUCLEIC_RESNAMES else _PROTEIN_BACKBONE
        names = [a.name for a in atoms]
        out.append({
            **_ref(key),
            "resname": resname,
            "atoms": names,
            # Páteř builder neumí postavit znovu - takové reziduum jde jen ponechat.
            "rebuildable": not any(n in backbone for n in names),
        })
    return out


# ---------------------------------------------------------------------------
# Otočení amidu ASN/GLN
# ---------------------------------------------------------------------------

def _partner_role(atom: _Atom) -> Optional[str]:
    """donor | acceptor | either | metal | None (nepolární / nezajímavý)."""
    if atom.record == "HETATM":
        if atom.resname in _WATER:
            return "either"
        if atom.element in _METAL_ELEMENTS and atom.name.upper() == atom.element:
            return "metal"
        return None
    if atom.name == "N" and atom.resname != "PRO":
        return "donor"
    if atom.name in ("O", "OXT"):
        return "acceptor"
    for table, role in ((_DONOR_ONLY, "donor"), (_ACCEPTOR_ONLY, "acceptor"), (_EITHER, "either")):
        if atom.name in table.get(atom.resname, ()):
            return role
    return None


def _site_score(position, as_oxygen: bool, partners: List[Tuple[_Atom, str]]) -> Tuple[float, List[Dict[str, Any]]]:
    """Skóre jednoho amidového atomu (O nebo N) na dané pozici + kontakty, které k němu přispěly."""
    score = 0.0
    contacts = []
    for atom, role in partners:
        d = _dist(position, atom.xyz)
        delta = 0.0
        if role == "metal":
            if d <= _METAL_MAX:
                delta = 2.0 if as_oxygen else -2.0
        elif d <= _HBOND_MAX:
            good = "donor" if as_oxygen else "acceptor"
            bad = "acceptor" if as_oxygen else "donor"
            if role == good:
                delta = 1.0
            elif role == "either":
                delta = 0.5
            elif role == bad and d <= _CLASH_MAX:
                delta = -1.0
        if delta:
            score += delta
            contacts.append({
                "partner": f"{_label(atom.key)}:{atom.name}",
                "partner_resname": atom.resname,
                "role": role,
                "distance_angstrom": round(d, 2),
                "favorable": delta > 0,
            })
    return score, contacts


def find_amide_flips(pdb_text: str, decided: Iterable[ResidueKey] = ()) -> List[Dict[str, Any]]:
    """
    ASN/GLN, u kterých má orientace amidu v PDB nevýhodný kontakt a otočená
    vychází lépe (alespoň o jeden bod). Skóre: vodíková vazba se správným partnerem +1,
    s partnerem, který umí obojí (voda, OH, HIS) +0.5, dva akceptory nebo dva
    donory proti sobě -1, kov u kyslíku +2 / u dusíku -2.
    """
    decided = set(decided)
    atoms = _parse_atoms(pdb_text)
    polar = [(a, role) for a in atoms if not a.is_hydrogen for role in [_partner_role(a)] if role]

    residues: Dict[ResidueKey, Dict[str, _Atom]] = {}
    for a in atoms:
        if a.record == "ATOM" and a.resname in _AMIDES:
            residues.setdefault(a.key, {})[a.name] = a

    out = []
    for key, res in residues.items():
        if key in decided:
            continue
        resname = next(iter(res.values())).resname
        o_name, n_name, _h = _AMIDES[resname]
        o_atom, n_atom = res.get(o_name), res.get(n_name)
        if o_atom is None or n_atom is None:
            continue
        # Domodelované atomy řeší sekce nulové obsazenosti.
        if not o_atom.occupancy or not n_atom.occupancy:
            continue

        partners = [
            (a, role) for a, role in polar
            if a.key != key and min(_dist(a.xyz, o_atom.xyz), _dist(a.xyz, n_atom.xyz)) <= _HBOND_MAX
        ]
        if not partners:
            continue
        o_now, o_now_c = _site_score(o_atom.xyz, True, partners)
        n_now, n_now_c = _site_score(n_atom.xyz, False, partners)
        o_flip, o_flip_c = _site_score(n_atom.xyz, True, partners)
        n_flip, n_flip_c = _site_score(o_atom.xyz, False, partners)
        current, flipped = o_now + n_now, o_flip + n_flip
        # Jen když současná orientace má skutečný problém (dva akceptory /
        # donory proti sobě, kov u dusíku) - rozdíl daný jen vodami nebo
        # OH skupinami (umí obojí) není dost silný důkaz.
        has_bad_contact = any(not c["favorable"] for c in o_now_c + n_now_c)
        if not has_bad_contact or flipped - current < 1.0:
            continue
        out.append({
            **_ref(key),
            "resname": resname,
            "oxygen": o_name,
            "nitrogen": n_name,
            "score_current": current,
            "score_flipped": flipped,
            "contacts_current": [{"site": o_name, **c} for c in o_now_c] + [{"site": n_name, **c} for c in n_now_c],
            "contacts_flipped": [{"site": o_name, **c} for c in o_flip_c] + [{"site": n_name, **c} for c in n_flip_c],
        })
    return out


# ---------------------------------------------------------------------------
# Aplikace rozhodnutí
# ---------------------------------------------------------------------------

def _with_atom_name(line: str, name: str, element: str) -> str:
    # 4znakové jméno od sloupce 13, kratší od 14 (PDB konvence pro 1písmenné prvky).
    field = name.ljust(4) if len(name) == 4 else f" {name}".ljust(4)
    line = line[:12] + field + line[16:]
    if len(line) >= 78:
        line = line[:76] + element.rjust(2) + line[78:]
    return line


def apply_structure_edits(
    pdb_text: str,
    flip_residues: Iterable[ResidueKey] = (),
    rebuild_residues: Iterable[ResidueKey] = (),
) -> str:
    """
    - flip: prohodí jména O a N amidu (souřadnice zůstávají, takže se atomy
      vymění) a zahodí vodíky na dusíku - builder je postaví znovu;
    - rebuild: zahodí atomy s nulovou obsazeností (mimo páteř), builder je
      dostaví jako chybějící postranní řetězec.
    """
    flips: Set[ResidueKey] = set(flip_residues)
    rebuilds: Set[ResidueKey] = set(rebuild_residues)
    if not flips and not rebuilds:
        return pdb_text

    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                key = (line[21].strip() or "?", int(line[22:26]), line[26].strip())
            except ValueError:
                out.append(line)
                continue
            resname = line[17:20].strip()
            name = line[12:16].strip()

            if key in rebuilds:
                backbone = _NUCLEIC_BACKBONE if resname in _NUCLEIC_RESNAMES else _PROTEIN_BACKBONE
                try:
                    zero = float(line[54:60]) == 0.0
                except ValueError:
                    zero = False
                if zero and name not in backbone:
                    continue

            if key in flips and resname in _AMIDES:
                o_name, n_name, hydrogens = _AMIDES[resname]
                if name in hydrogens:
                    continue
                if name == o_name:
                    line = _with_atom_name(line, n_name, "N")
                elif name == n_name:
                    line = _with_atom_name(line, o_name, "O")
        out.append(line)
    return "\n".join(out)
