"""
Citace vstupní struktury pro report a popis postupu.

Primární článek ke struktuře je v hlavičce PDB v záznamech JRNL (AUTH, TITL,
REF, PMID, DOI). Čte se při načtení struktury - po přípravě už hlavička ve
structure.pdb není - a ukládá se do záznamu SOURCE. Kód PDB se bere z HEADER,
takže DOI záznamu v PDB jde uvést i u nahraného souboru.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# Citace samotné databáze, kterou RCSB žádá u každé stažené struktury.
RCSB_PDB_REFERENCE = ("Berman HM, Westbrook J, Feng Z, Gilliland G, Bhat TN, Weissig H, Shindyalov IN, Bourne PE. "
                      "The Protein Data Bank. Nucleic Acids Res. 2000;28:235-242. doi:10.1093/nar/28.1.235.")

_REF_RE = re.compile(r"^(?P<journal>.*?)\s+V\.\s*(?P<volume>\S*)\s+(?P<page>\S*)\s+(?P<year>\d{4})\s*$")


def _author(raw: str) -> str:
    """'D.J.KLEIN' -> 'Klein DJ'."""
    initials, _, surname = raw.strip().rpartition(".")
    if not surname:
        return raw.strip().title()
    return f"{surname.strip().title()} {initials.replace('.', '').replace('-', '')}".strip()


def parse_jrnl(pdb_text: str) -> Dict[str, Any]:
    """Primární citace a kód PDB z hlavičky. Co v souboru není, ve výsledku chybí."""
    fields: Dict[str, List[str]] = {}
    pdb_id = ""
    for line in pdb_text.splitlines():
        record = line[:6]
        if record == "HEADER":
            pdb_id = line[62:66].strip()
        elif record == "JRNL  ":
            fields.setdefault(line[12:16].strip(), []).append(line[19:79].strip())
        elif record in ("ATOM  ", "HETATM", "MODEL "):
            break

    out: Dict[str, Any] = {}
    if pdb_id:
        out["pdb_id"] = pdb_id.upper()
    if fields.get("AUTH"):
        # Pokračovací řádek navazuje za čárkou, jména se nedělí.
        out["authors"] = [_author(a) for a in "".join(fields["AUTH"]).split(",") if a.strip()]
    if fields.get("TITL"):
        out["title"] = " ".join(fields["TITL"]).rstrip(".")
    ref = " ".join(fields.get("REF") or [])
    match = _REF_RE.match(ref)
    if match:
        out.update({k: v for k, v in match.groupdict().items() if v})
    elif ref:
        out["journal"] = ref  # "TO BE PUBLISHED"
    if fields.get("PMID"):
        out["pmid"] = fields["PMID"][0]
    if fields.get("DOI"):
        out["doi"] = fields["DOI"][0].lower()
    return out


def pdb_entry_doi(pdb_id: str) -> str:
    return f"10.2210/pdb{pdb_id.lower()}/pdb"


def structure_pdb_id(source: Dict[str, Any]) -> str:
    return (source.get("pdb_code") or (source.get("citation") or {}).get("pdb_id") or "").upper()


def format_article(citation: Dict[str, Any]) -> Optional[str]:
    """Primární článek ve stylu Vancouver; None, když hlavička citaci nemá."""
    if not citation.get("title") and not citation.get("authors"):
        return None
    parts = []
    authors = citation.get("authors") or []
    if authors:
        parts.append(", ".join(authors[:6]) + (", et al" if len(authors) > 6 else ""))
    if citation.get("title"):
        parts.append(citation["title"])
    journal = citation.get("journal")
    if journal:
        details = journal.rstrip(".")
        if citation.get("year"):
            details += f" {citation['year']}"
        if citation.get("volume"):
            details += f";{citation['volume']}"
            if citation.get("page"):
                details += f":{citation['page']}"
        parts.append(details)
    if citation.get("doi"):
        parts.append(f"doi:{citation['doi']}")
    if citation.get("pmid"):
        parts.append(f"PMID:{citation['pmid']}")
    return ". ".join(p.rstrip(".") for p in parts) + "."


def structure_references(source: Optional[Dict[str, Any]]) -> List[str]:
    """Citace vstupu v pořadí: článek ke struktuře, záznam v PDB, databáze PDB (jen u stažených)."""
    if not source:
        return []
    entries = []
    article = format_article(source.get("citation") or {})
    if article:
        entries.append(article)
    pdb_id = structure_pdb_id(source)
    if pdb_id:
        entries.append(f"PDB entry {pdb_id}. doi:{pdb_entry_doi(pdb_id)}.")
    if source.get("kind") == "pdb":
        entries.append(RCSB_PDB_REFERENCE)
    return entries
