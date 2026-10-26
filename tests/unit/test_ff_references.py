"""
Citace force fieldů - app/services/report/ff_references.py.
"""

from app.services.report.ff_references import ARTICLES, FF_REFERENCES, ff_references


def test_every_ff_cites_known_articles():
    for name, keys in FF_REFERENCES.items():
        assert keys, name
        assert all(k in ARTICLES for k in keys), name


def test_table_wins_over_catalog():
    refs = ff_references({"ff_name": "ol3-ljbb", "reference_article_doi": "10.1/wrong", "creators": ["X"]})

    assert [r.split(".")[0] for r in refs] == [
        "Cornell et al", "Wang et al", "Pérez et al", "Zgarbová et al", "Steinbrecher et al", "Bergonzo and Cheatham (2015)",
    ]
    assert refs[-1].endswith("doi:10.1021/acs.jctc.5b00444.")
    assert not any("10.1/wrong" in r for r in refs)


def test_desres_has_full_doi():
    assert ff_references({"display_name": "DESRES"}) == [
        "Tan et al. (2018). RNA Force Field with Accuracy Comparable to State-of-the-Art Protein Force Fields. "
        "doi:10.1073/pnas.1713027115.",
    ]


def test_other_ff_falls_back_to_catalog():
    refs = ff_references({"ff_name": "FF14SB", "reference_article_doi": "10.1021/acs.jctc.5b00255",
                          "creators": ["Maier, James A.", "Simmerling, Carlos"]})

    assert refs == ["Maier, James A.; Simmerling, Carlos. doi:10.1021/acs.jctc.5b00255."]
    assert ff_references({"ff_name": "X"}) == []
