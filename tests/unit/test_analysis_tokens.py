"""
Unit testy pro app/services/analysis_service.py::build_sequence_tokens.

Pokrývá:
- detekci sekvenčních děr (is_gap tokeny) a rozlišení "opravdové" díry od
  pouhého přejmenovaného rezidua
- terminalitu na okraji díry (gap -> C-terminál/3'-terminál / N-terminál/5'-terminál)
  - toto je regresní sada pro opravu popsanou v konverzaci: GLU83/PHE89 a
    U125/A128 příklady vycházejí z reálného PDB 1JJ2
- terminalitu na skutečném konci řetězce (beze změny oproti původnímu chování)
- detekci HIS variant podle přítomných vodíků

Nahrazuje původní tests/test_gap_fix.py a tests/test_real_gap.py, které byly
jen debug skripty bez assercí.
"""

import pytest

from app.services.analysis_service import build_sequence_tokens

pytestmark = pytest.mark.unit


def _tokens(pdb_text: str, chain: str, fill_gaps: bool = True):
    result = build_sequence_tokens(pdb_text, chain=chain, fill_gaps=fill_gaps)
    return result["chains"][chain]["tokens"]


def _by_resseq(tokens, resseq):
    return next(t for t in tokens if not t["is_gap"] and t["resseq"] == resseq)


class TestGapDetection:
    def test_no_gap_when_residue_only_renamed_not_missing(self):
        """
        Pokud je číslování souvislé (i když je prostřední reziduum neznámého
        jména), nejde o sekvenční díru - žádný is_gap token se nevloží.
        """
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   GAL A   2       2.000   2.000   2.000  1.00  0.00           N\n"
            "ATOM      4  CA  GAL A   2       3.000   3.000   3.000  1.00  0.00           C\n"
            "ATOM      5  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      6  CA  ALA A   3       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert sum(1 for t in tokens if t["is_gap"]) == 0

    def test_gap_token_inserted_when_residue_truly_missing(self):
        """Reziduum 2 v PDB vůbec neexistuje -> mezi 1 a 3 se vloží is_gap token."""
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      4  CA  ALA A   3       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert sum(1 for t in tokens if t["is_gap"]) == 1

    def test_fill_gaps_false_disables_detection(self):
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N\n"
        )
        tokens = _tokens(pdb, "A", fill_gaps=False)
        assert sum(1 for t in tokens if t["is_gap"]) == 0


class TestGapTerminality:
    """Regresní sada pro terminalitu na okraji sekvenční díry (viz PDB 1JJ2)."""

    def test_protein_gap_produces_c_and_n_terminal_variants(self, pdb_protein_gap):
        tokens = _tokens(pdb_protein_gap, "K")

        glu83 = _by_resseq(tokens, 83)
        phe89 = _by_resseq(tokens, 89)

        assert glu83["ff_resname"] == "CGLU"
        assert glu83["terminus_reason"] == "gap"
        assert phe89["ff_resname"] == "NPHE"
        assert phe89["terminus_reason"] == "gap"

        # mezi 84 a 88 chybí 5 reziduí -> 5 gap tokenů
        assert sum(1 for t in tokens if t["is_gap"]) == 5

    def test_rna_gap_produces_3_and_5_prime_variants(self, pdb_rna_gap):
        # Fixtura je výřez z reálné struktury 1RNA (chain A, residua 1-3 a 6-8,
        # 4-5 vynechána) - reálná geometrie, ne syntetické souřadnice, protože
        # FORGE builder na degenerovaných/kolineárních atomech spadne na
        # "Cannot normalize near-zero central dihedral bond".
        tokens = _tokens(pdb_rna_gap, "A")

        a3 = _by_resseq(tokens, 3)
        u6 = _by_resseq(tokens, 6)

        assert a3["ff_resname"] == "RA3"
        assert a3["terminus_reason"] == "gap"
        assert u6["ff_resname"] == "RU5"
        assert u6["terminus_reason"] == "gap"
        assert sum(1 for t in tokens if t["is_gap"]) == 2

    def test_gap_terminus_warning_is_emitted(self, pdb_protein_gap):
        result = build_sequence_tokens(pdb_protein_gap, chain="K", fill_gaps=True)
        warnings = result["chains"]["K"]["warnings"]
        assert any("artificial C-terminus" in w for w in warnings)
        assert any("artificial N-terminus" in w for w in warnings)

    def test_regular_chain_ends_unaffected_by_gap_logic(self):
        """
        Regrese: obyčejný 3reziduový řetězec beze změny - kraje dostanou
        chain_end terminalitu, prostředek žádnou.
        """
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   GLY A   2       2.000   2.000   2.000  1.00  0.00           N\n"
            "ATOM      4  CA  GLY A   2       3.000   3.000   3.000  1.00  0.00           C\n"
            "ATOM      5  N   SER A   3       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      6  CA  SER A   3       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        ala1, gly2, ser3 = (_by_resseq(tokens, i) for i in (1, 2, 3))

        assert ala1["ff_resname"] == "NALA"
        assert ala1["terminus_reason"] == "chain_end"
        assert gly2["ff_resname"] == "GLY"
        assert gly2["terminus_reason"] is None
        assert ser3["ff_resname"] == "CSER"
        assert ser3["terminus_reason"] == "chain_end"


class TestGapBoundaryExclusion:
    """
    GLU83 v reálném 1JJ2 zůstává v ATOM záznamech jen s N/CA/C/O/CB -
    CG/CD/OE1/OE2 chybí i po přeznačení na CGLU (ověřeno přímo na staženém
    1JJ2.pdb z RCSB, ne jen odhadem). Backbone kotva (C) ale přítomná je, takže
    reziduum se v modelu ponechá - nový builder (v1.0.0) umí tenhle typ mezery
    (jeden kotevní bod, chybí jen distální postranní řetězec) bezpečně
    dostavit interaktivně přes GUI (viz app/builder/INTEGRATION_CONTRACT.md,
    "residue_local_open_branch"), takže tichá exkluze celého rezidua by teď
    jen zahodila data, která builder umí zpracovat. Vylučuje se ONLY reziduum,
    kterému chybí i tahle backbone kotva sama (viz
    _GAP_BOUNDARY_ANCHOR_ATOM v analysis_service.py) - to je jediný případ,
    kdy builder opravdu nemá nic, k čemu by cokoliv připojil.
    """

    def test_incomplete_boundary_residue_with_anchor_is_kept(
        self, pdb_protein_gap_incomplete_boundary
    ):
        tokens = _tokens(pdb_protein_gap_incomplete_boundary, "K")

        glu83 = _by_resseq(tokens, 83)
        assert glu83["ff_resname"] == "CGLU"
        assert glu83["terminus_reason"] == "gap"
        # Backbone kotva (C) přítomná, distální postranní řetězec (CG/CD/OE1/
        # OE2) chybí - přesně ten případ, který teď řeší interaktivní GUI
        # místo tiché exkluze.
        assert glu83["gap_boundary_incomplete"] is False
        assert "CG" in glu83["missing_atoms"]

        # ALA82 (první reziduum tohoto fragmentu, tedy přirozený N-konec) se
        # terminalitou kvůli gapu NEPŘEZNAČUJE - ta zůstává na GLU83.
        ala82 = _by_resseq(tokens, 82)
        assert ala82["terminus_reason"] != "gap"

        phe89 = _by_resseq(tokens, 89)
        assert phe89["ff_resname"] == "NPHE"

    def test_incomplete_boundary_residue_warns_instead_of_excluding(
        self, pdb_protein_gap_incomplete_boundary
    ):
        result = build_sequence_tokens(pdb_protein_gap_incomplete_boundary, chain="K", fill_gaps=True)
        warnings = result["chains"]["K"]["warnings"]
        assert not any("Excluded" in w for w in warnings)
        assert any("K:83" in w and "interactive side-chain completion" in w for w in warnings)

    def test_boundary_residue_missing_its_own_anchor_is_still_excluded(self):
        # ASP83 tady nemá ani C - i s novým builderem nemá vůbec nic, k čemu
        # by se cokoliv připojilo, takže musí zůstat vyloučené jako dřív.
        pdb = (
            "ATOM      1  N   ALA K  82       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA K  82       1.500   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA K  82       3.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA K  82       3.500   1.200   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA K  82      -0.200   0.800  -1.300  1.00  0.00           C\n"
            "ATOM      6  N   ASP K  83       4.500  -1.000   0.000  1.00  0.00           N\n"
            "ATOM      7  CA  ASP K  83       6.000  -1.000   0.000  1.00  0.00           C\n"
            "ATOM      8  N   PHE K  89       8.000  -3.000   0.000  1.00  0.00           N\n"
            "ATOM      9  CA  PHE K  89       9.500  -3.000   0.000  1.00  0.00           C\n"
            "ATOM     10  C   PHE K  89      11.000  -3.000   0.000  1.00  0.00           C\n"
            "ATOM     11  O   PHE K  89      11.500  -1.800   0.000  1.00  0.00           O\n"
        )
        tokens = _tokens(pdb, "K")
        assert all(t.get("resseq") != 83 for t in tokens if not t["is_gap"])

        ala82 = _by_resseq(tokens, 82)
        assert ala82["ff_resname"] == "CALA"
        assert ala82["terminus_reason"] == "gap"

    def test_complete_boundary_residue_is_not_excluded(self, pdb_protein_gap):
        # protein_gap.pdb má GLU83 s kompletním postranním řetězcem - žádné
        # vylučování se nesmí spustit, chování zůstává jako předtím.
        tokens = _tokens(pdb_protein_gap, "K")
        glu83 = _by_resseq(tokens, 83)
        assert glu83["ff_resname"] == "CGLU"
        assert glu83["gap_boundary_incomplete"] is False

    def test_boundary_residue_with_only_hydrogens_missing_is_not_flagged(self, pdb_rna_gap):
        tokens = _tokens(pdb_rna_gap, "A")
        a3 = _by_resseq(tokens, 3)
        assert a3["gap_boundary_incomplete"] is False


class TestNonNumericChainBreaks:
    """
    Sekvenční díra odvozená z čísel reziduí není jediný důkaz přerušení
    polymeru - viz INTEGRATION_CONTRACT.md. Explicitní TER uprostřed řetězce
    nebo chemicky nemožná meziresiduová vzdálenost musí vést ke stejné
    terminální úpravě, i když číslování zůstává souvislé.
    """

    # ALA1 musí mít kompletní těžké atomy (N/CA/C/O/CB) - jinak by ho nové
    # kaskádové vylučování neúplných okrajových reziduí (viz
    # TestGapBoundaryExclusion) samo vyřadilo z modelu ještě předtím, než se
    # vůbec dostaneme k ověřování TER/geometrického přerušení.
    _PROTEIN_TER_BREAK = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       3.500   1.200   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.200  -1.200   0.000  1.00  0.00           C\n"
        "TER       6      ALA A   1\n"
        "ATOM      7  N   GLY A   2       4.330   0.000   0.000  1.00  0.00           N\n"
        "ATOM      8  CA  GLY A   2       5.830   0.000   0.000  1.00  0.00           C\n"
    )

    _PROTEIN_GEOMETRY_BREAK = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       3.500   1.200   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       2.200  -1.200   0.000  1.00  0.00           C\n"
        "ATOM      6  N   GLY A   2      50.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      7  CA  GLY A   2      51.500   0.000   0.000  1.00  0.00           C\n"
    )

    def test_mid_chain_ter_forces_termini_without_numbering_gap(self):
        tokens = _tokens(self._PROTEIN_TER_BREAK, "A")
        ala1, gly2 = _by_resseq(tokens, 1), _by_resseq(tokens, 2)
        assert sum(1 for t in tokens if t["is_gap"]) == 0
        assert ala1["ff_resname"] == "CALA" and ala1["terminus_reason"] == "ter"
        assert gly2["ff_resname"] == "NGLY" and gly2["terminus_reason"] == "ter"

    def test_implausible_bond_distance_forces_termini_without_numbering_gap(self):
        tokens = _tokens(self._PROTEIN_GEOMETRY_BREAK, "A")
        ala1, gly2 = _by_resseq(tokens, 1), _by_resseq(tokens, 2)
        assert sum(1 for t in tokens if t["is_gap"]) == 0
        assert ala1["ff_resname"] == "CALA" and ala1["terminus_reason"] == "geometry"
        assert gly2["ff_resname"] == "NGLY" and gly2["terminus_reason"] == "geometry"

    def test_plausible_bond_distance_without_ter_is_not_flagged(self):
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      4  N   GLY A   2       4.330   0.000   0.000  1.00  0.00           N\n"
            "ATOM      5  CA  GLY A   2       5.830   0.000   0.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert _by_resseq(tokens, 2)["terminus_reason"] is None


class TestZeroResidueNumberSkip:
    """
    Regresní sada pro opravu popsanou v konverzaci: reálné PDB 2OUE má
    chain A číslovaný -5..-2 (uměle přidaný 5' leader), pak podle běžné PDB
    konvence PŘESKOČÍ sekvenční číslo 0 a pokračuje 1, 2, ... Sekvenční
    číslo 0 v PDB záznamech principiálně nikdy neexistuje, takže tenhle
    skok není skutečná mezera - kód dřív "reziduum 0" bral jako chybějící a
    vytvořil falešný GAP placeholder token (pdb_resname="0", resseq=None,
    zobrazený v UI jako "UNK ??").
    """

    def test_no_gap_when_numbering_skips_zero_by_convention(self):
        pdb = (
            "ATOM      1  N   ALA A  -1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A  -1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   ALA A   1       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      4  CA  ALA A   1       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert sum(1 for t in tokens if t["is_gap"]) == 0

    def test_real_gap_still_detected_when_crossing_zero_boundary(self):
        """
        Reziduum 1 opravdu v PDB chybí (skok rovnou z -1 na 2). Přeskočení
        neexistujícího čísla 0 nesmí "spolknout" i tenhle skutečně chybějící
        úsek - musí vzniknout přesně jeden GAP token, ne nula.
        """
        pdb = (
            "ATOM      1  N   ALA A  -1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A  -1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   ALA A   2       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      4  CA  ALA A   2       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert sum(1 for t in tokens if t["is_gap"]) == 1


class TestExtraAtomDetection:
    """
    _check_extra_atoms je zrcadlová funkce k _check_missing_atoms: místo
    "co v šabloně je, ale v PDB chybí" hlídá "co je v PDB navíc oproti
    šabloně". Builder tenhle případ interně eviduje jako
    observed_extra_atoms (viz forge_molecule_parser.py) - tenhle token je
    analytická vrstva nad stejným konceptem.
    """

    _ALA_WITH_EXTRA = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       2.500   1.000   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1      -0.500   1.000  -1.000  1.00  0.00           C\n"
        "ATOM      6  CX  ALA A   1       9.000   9.000   9.000  1.00  0.00           C\n"
        "ATOM      7  HZ9 ALA A   1       9.500   9.500   9.500  1.00  0.00           H\n"
    )

    def test_no_extra_atoms_for_clean_residue(self):
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       2.500   1.000   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA A   1      -0.500   1.000  -1.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert _by_resseq(tokens, 1)["extra_atoms"] == []

    def test_heavy_and_hydrogen_extra_atoms_are_both_reported(self):
        tokens = _tokens(self._ALA_WITH_EXTRA, "A")
        ala1 = _by_resseq(tokens, 1)
        assert ala1["extra_atoms"] == ["CX", "HZ9"]

    def test_extra_atom_warning_is_emitted(self):
        result = build_sequence_tokens(self._ALA_WITH_EXTRA, chain="A", fill_gaps=True)
        warnings = result["chains"]["A"]["warnings"]
        assert any("Extra atoms not in template" in w and "CX" in w for w in warnings)

    def test_5prime_gap_phosphate_is_not_flagged_as_extra(self, pdb_rna_gap):
        """
        Regrese: reziduum hned za mezerou dostane 5'-terminální variantu
        (RU5/RA5/...), jejíž FF šablona fosfátovou skupinu neobsahuje - ale
        reálná struktura ji stejně nese (fosfát, který by se normálně
        napojoval na chybějící předchozí reziduum). To NENÍ neznámý/extra
        atom, je to očekávaný stav umělého 5' konce - viz
        _TERMINAL_ALLOWED_EXTRA_HEAVY_ATOMS.
        """
        tokens = _tokens(pdb_rna_gap, "A")
        u6 = _by_resseq(tokens, 6)
        assert u6["ff_resname"] == "RU5"
        assert u6["terminus_reason"] == "gap"
        assert "P" in u6["atoms"]
        assert u6["extra_atoms"] == []

    def test_unknown_residue_reports_no_extra_atoms(self):
        """Bez šablony (neznámé reziduum) nelze nic porovnávat - stejně jako u missing_atoms."""
        pdb = (
            "ATOM      1  N   XYZ A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  XYZ A   1       1.000   1.000   1.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        xyz1 = _by_resseq(tokens, 1)
        assert xyz1["known"] is False
        assert xyz1["extra_atoms"] == []


class TestHisVariantDetection:
    # Bare "HIS" není samo o sobě klíčem v converting_dictionary.json (jen
    # jeho vyřešené HID/HIE/HIP varianty jsou) - _infer_group má proto pro
    # "HIS" explicitní výjimku vracející "P" (viz komentář u _infer_group).
    # Bez ní _pick_variant/_get_res_def sice i tak správně dohledaly
    # HID/HIE/HIP podle přítomných vodíků, ale reziduum s group=None
    # nespadlo do main_chain a builder ho dostal jako nesouvisející HETATM
    # bez napojení na sousední rezidua - potvrzeno pádem na reálném 1JJ2
    # ("Passthrough atom HIS:N ... lacks converting identity" při
    # solvataci). To porušovalo invariantu #5 z INTEGRATION_CONTRACT.md
    # (generické HIS musí být vyřešené na HID/HIE/HIP PŘI ZACHOVÁNÍ původní
    # pozice v řetězci, ne jen se správným ff_resname).

    def _his_pdb(self, atoms):
        # HIS je zasazené mezi dvě ALA (resseq 1 a 3), ať je jednoznačně
        # mid-chain (group="P" ho teď od opravy _infer_group správně řadí
        # do main_chain) - jinak by jako jediné/první reziduum svého řetězce
        # samo dostalo N-terminální variantu (NHID/NHIE/NHIP) a testovalo by
        # se tím terminus_reason, ne rozpoznání tautomeru podle vodíků.
        # Páteřní C/N atomy musí být v bonding distance mezi rezidui - jinak
        # by je nový geometrický break-detektor (viz TestNonNumericChainBreaks)
        # falešně vyhodnotil jako přerušení i v tomhle obyčejném souvislém
        # řetězci. Postranní řetězec HIS na přesné souřadnici nezáleží (do
        # geometrické kontroly vstupují jen atomy C/N/O3'/P), proto zůstává
        # na degenerovaném společném bodě jako v ostatních fixturách v tomhle
        # souboru.
        lines = [
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
            "ATOM      2  CA  ALA A   1       0.500   0.000   0.000  1.00  0.00           C",
            "ATOM      3  C   ALA A   1       1.000   0.000   0.000  1.00  0.00           C",
            "ATOM      4  O   ALA A   1       1.200   0.500   0.000  1.00  0.00           O",
            "ATOM      5  CB  ALA A   1       0.700  -0.500   0.000  1.00  0.00           C",
        ]
        backbone = {
            "N": (2.330, 0.000, 0.000),
            "CA": (2.830, 0.000, 0.000),
            "C": (3.330, 0.000, 0.000),
            "O": (3.530, 0.500, 0.000),
            "CB": (3.030, -0.500, 0.000),
        }
        for i, name in enumerate(["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"] + atoms, start=6):
            coord = backbone.get(name, backbone["CB"])
            lines.append(
                f"ATOM  {i:>5} {name:<4} HIS A   2      {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                f"  1.00  0.00           {name[0]}"
            )
        next_i = 6 + 10 + len(atoms)
        lines.append(f"ATOM  {next_i:>5} N    ALA A   3       4.660   0.000   0.000  1.00  0.00           N")
        lines.append(f"ATOM  {next_i + 1:>5} CA   ALA A   3       5.160   0.000   0.000  1.00  0.00           C")
        return "\n".join(lines)

    def _his_token(self, atoms):
        return next(t for t in _tokens(self._his_pdb(atoms), "A") if t["pdb_resname"] == "HIS")

    def test_hip_when_both_hydrogens_present(self):
        his = self._his_token(["HD1", "HE2"])
        assert his["ff_resname"] == "HIP"
        assert his["group"] == "P"
        assert his["terminus_reason"] is None

    def test_hie_when_only_epsilon_hydrogen_present(self):
        his = self._his_token(["HE2"])
        assert his["ff_resname"] == "HIE"
        assert his["group"] == "P"

    def test_hid_default_when_no_extra_hydrogens(self):
        his = self._his_token([])
        assert his["ff_resname"] == "HID"
        assert his["group"] == "P"

    def test_his_in_middle_of_chain_stays_in_main_chain(self):
        """
        HIS obklopené jinými rezidui musí zůstat součástí main_chain, ať se
        na něm nepřeruší kontinuita řetězce - flankující rezidua smí dostat
        terminální variantu jen na skutečném konci řetězce, nikde jinde.
        """
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C\n"
            "ATOM      3  N   HIS A   2       2.000   2.000   2.000  1.00  0.00           N\n"
            "ATOM      4  CA  HIS A   2       3.000   3.000   3.000  1.00  0.00           C\n"
            "ATOM      5  N   ALA A   3       4.000   4.000   4.000  1.00  0.00           N\n"
            "ATOM      6  CA  ALA A   3       5.000   5.000   5.000  1.00  0.00           C\n"
        )
        tokens = _tokens(pdb, "A")
        assert sum(1 for t in tokens if t["is_gap"]) == 0
        his2 = _by_resseq(tokens, 2)
        assert his2["group"] == "P" and his2["terminus_reason"] is None
        assert _by_resseq(tokens, 1)["terminus_reason"] == "chain_end"
        assert _by_resseq(tokens, 3)["terminus_reason"] == "chain_end"
