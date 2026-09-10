"""
Unit testy pro app/services/forcefield_service.py.

Běží izolovaně přes monkeypatch CACHE_DIR/FORGE_CACHE_DIR do tmp_path -
NIKDY nesahají na skutečnou data/ff_cache/ (reálná, draze stažená data).
"""

import base64

import pytest

from app.services.forcefield_service import ForceFieldService

pytestmark = pytest.mark.unit


@pytest.fixture
def service(tmp_path, monkeypatch):
    monkeypatch.setattr(ForceFieldService, "CACHE_DIR", tmp_path / "ff_cache")
    monkeypatch.setattr(ForceFieldService, "FORGE_CACHE_DIR", tmp_path / "ff_cache_forge")
    return ForceFieldService()


def _fake_ff_data(name="TESTFF"):
    return {
        "display_name": name,
        "force_field_file": "[ defaults ]\n1 2 yes 0.5 0.8333\n#include nonbonded\n",
        "residue_lib_ff_file": "[ RES ]\n [ atoms ]\n  A1  T1  0.5  1\n",
        "nonbonded_ff_file": "[ atomtypes ]\nT1 6 12.0 0.0 A 0.35 0.4\n",
        "bonded_ff_file": "[ bondtypes ]\nT1 T1 1 0.15 1000\n",
        "atom_type_ff_file": "T1 12.0\n",
    }


class TestDecodeContent:
    def test_decodes_valid_base64(self, service):
        original = "hello world"
        encoded = base64.b64encode(original.encode()).decode()
        assert service._decode_content(encoded) == original

    def test_strips_data_uri_prefix_before_decoding(self, service):
        original = "hello world"
        encoded = "data:text/plain;base64," + base64.b64encode(original.encode()).decode()
        assert service._decode_content(encoded) == original

    def test_falls_back_to_plain_text_when_not_valid_base64(self, service):
        plain = "[ defaults ]!!not-base64-length"
        assert service._decode_content(plain) == plain

    def test_empty_or_missing_content_returns_empty_string(self, service):
        assert service._decode_content(None) == ""
        assert service._decode_content("") == ""


class TestPrepareForcefieldFiles:
    def test_writes_expected_files_to_cache_dir(self, service):
        target_dir = service.prepare_forcefield_files(_fake_ff_data("TESTFF"))
        assert target_dir == service.CACHE_DIR / "TESTFF"
        for fname in ["TESTFF.rtp", "nonbonded_TESTFF.itp", "bonded_TESTFF.itp", "TESTFF.atp"]:
            assert (target_dir / fname).exists()

    def test_flattens_include_directives(self, service):
        target_dir = service.prepare_forcefield_files(_fake_ff_data("TESTFF"))
        rtp_content = (target_dir / "TESTFF.rtp").read_text(encoding="utf-8")
        assert "#include" not in rtp_content
        assert "[ atomtypes ]" in rtp_content  # obsah z nonbonded_ff_file byl vlepen misto #include

    def test_display_name_spaces_become_underscores(self, service):
        target_dir = service.prepare_forcefield_files(_fake_ff_data("My Force Field"))
        assert target_dir.name == "My_Force_Field"


class TestPrepareForgeForceFieldDirectory:
    def test_creates_mol_type_suffixed_directory(self, service):
        target_dir = service.prepare_forge_force_field_directory(_fake_ff_data("TESTFF"), "R")
        assert target_dir == service.FORGE_CACHE_DIR / "TESTFF_R"
        assert target_dir.exists()

    def test_contains_files_matching_builder_naming_contract(self, service):
        """
        app/builder/forge_molecule_solvation.py::_one_matching_file hledá soubory
        podle substringu "residue-lib"/"residue_lib", "nonbonded", "bonded"
        (bez "nonbonded"), "forcefield"/"force_field".
        """
        target_dir = service.prepare_forge_force_field_directory(_fake_ff_data("TESTFF"), "R")
        names = {p.name for p in target_dir.iterdir()}
        assert any("residue_lib" in n for n in names)
        assert any("forcefield" in n for n in names)
        assert any(n.startswith("nonbonded_") for n in names)
        assert any(n.startswith("bonded_") and not n.startswith("nonbonded_") for n in names)

    def test_two_different_mol_types_get_separate_directories(self, service):
        dir_r = service.prepare_forge_force_field_directory(_fake_ff_data("TESTFF"), "R")
        dir_d = service.prepare_forge_force_field_directory(_fake_ff_data("TESTFF"), "D")
        assert dir_r != dir_d
        assert dir_r.exists() and dir_d.exists()

    def test_does_not_touch_unrelated_cache_entries(self, service):
        """Volání pro jeden FF nesmí nijak ovlivnit jiný, dřív připravený FF."""
        service.prepare_forge_force_field_directory(_fake_ff_data("FIRST"), "R")
        service.prepare_forge_force_field_directory(_fake_ff_data("SECOND"), "R")
        assert (service.FORGE_CACHE_DIR / "FIRST_R").exists()
        assert (service.FORGE_CACHE_DIR / "SECOND_R").exists()


class TestClearCache:
    def test_clears_both_cache_directories(self, service):
        service.prepare_forge_force_field_directory(_fake_ff_data("TESTFF"), "R")
        assert any(service.CACHE_DIR.iterdir())
        assert any(service.FORGE_CACHE_DIR.iterdir())

        service.clear_cache()

        assert list(service.CACHE_DIR.iterdir()) == []
        assert list(service.FORGE_CACHE_DIR.iterdir()) == []
