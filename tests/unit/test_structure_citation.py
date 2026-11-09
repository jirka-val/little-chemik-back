"""
Citace vstupní struktury z hlavičky PDB (JRNL) - app/services/report/citation.py.
"""

from pathlib import Path

from app.services.report.citation import (
    RCSB_PDB_REFERENCE, format_article, parse_jrnl, pdb_entry_doi, structure_references,
)

# Kopie z RCSB v repu (data/pdb_files je jen lokální cache stažených struktur, není v gitu).
PDB_FILES = Path(__file__).resolve().parents[1] / "fixtures" / "pdb" / "golden"


def test_parse_jrnl_joins_continuation_lines():
    citation = parse_jrnl((PDB_FILES / "1RNA.pdb").read_text())

    assert citation["pdb_id"] == "1RNA"
    assert citation["authors"] == [
        "Dock-Bregeon AC", "Chevrier B", "Podjarny A", "Johnson J",
        "De Bear JS", "Gough GR", "Gilham PT", "Moras D",
    ]
    assert citation["title"] == "CRYSTALLOGRAPHIC STRUCTURE OF AN RNA HELIX: [U(UA)6A]2"
    assert (citation["journal"], citation["volume"], citation["page"], citation["year"]) == \
        ("J.MOL.BIOL.", "209", "459", "1989")
    assert citation["doi"] == "10.1016/0022-2836(89)90010-7"
    assert citation["pmid"] == "2479753"


def test_format_article_vancouver_style():
    text = format_article(parse_jrnl((PDB_FILES / "1CRN.pdb").read_text()))
    assert text == ("Teeter MM. WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION: PENTAGON RINGS OF "
                    "WATER MOLECULES IN CRYSTALS OF CRAMBIN. PROC.NATL.ACAD.SCI.USA 1984;81:6014. "
                    "doi:10.1073/pnas.81.19.6014. PMID:16593516.")


def test_header_without_jrnl_still_gives_pdb_entry():
    pdb = "HEADER    RNA                                     01-FEB-90   1ABC              \nATOM      1  P"
    refs = structure_references({"kind": "upload", "filename": "x.pdb", "citation": parse_jrnl(pdb)})
    assert refs == [f"PDB entry 1ABC. doi:{pdb_entry_doi('1ABC')}."]


def test_downloaded_structure_also_cites_the_pdb():
    refs = structure_references({"kind": "pdb", "pdb_code": "1crn", "citation": {}})
    assert refs == ["PDB entry 1CRN. doi:10.2210/pdb1crn/pdb.", RCSB_PDB_REFERENCE]


def test_no_source_no_references():
    assert structure_references(None) == []
    assert structure_references({"kind": "upload", "filename": "x.pdb"}) == []
