"""
Unit testy pro app/services/ff_classification_service.py.

Běží nad izolovanou kopií data/force_fields.json v tmp_path - NIKDY nesahají
na skutečný vendorovaný soubor, aby test nemohl omylem přepsat produkční
klasifikaci.
"""

import json
import shutil
from pathlib import Path

import pytest

from app.services.ff_classification_service import ForceFieldClassificationService

pytestmark = pytest.mark.unit

_SOURCE_FILE = Path(__file__).resolve().parents[2] / "data" / "force_fields.json"


@pytest.fixture
def service(tmp_path) -> ForceFieldClassificationService:
    target = tmp_path / "force_fields.json"
    shutil.copy(_SOURCE_FILE, target)
    return ForceFieldClassificationService(path=target)


def _name(display_name: str) -> str:
    return display_name.replace(" ", "_")


class TestClassifySolute:
    def test_default_is_recommended_and_flagged(self, service):
        result = service.classify_solute("P", "FF14SB")
        assert result == {"tier": "recommended", "is_default": True}

    def test_supported_is_not_default(self, service):
        result = service.classify_solute("P", "FF12SB")
        assert result == {"tier": "supported", "is_default": False}

    def test_obsolete(self, service):
        assert service.classify_solute("P", "FF99SB")["tier"] == "obsolete"

    def test_unknown_ff_is_new_unclassified(self, service):
        assert service.classify_solute("P", "SomeBrandNewFF") == {
            "tier": "new_unclassified",
            "is_default": False,
        }


class TestIonTierMap:
    def test_default_water_type_ion_is_flagged_default(self, service):
        # solvent.default_water_type = W3, jehož default_profile (SPCE_JC_LM)
        # používá pro I1 "JC-SPCE_I1" - viz data/force_fields.json.
        entry = service.classify_ion("I1", "JC-SPCE_I1")
        assert entry == {"tier": "recommended", "is_default": True}

    def test_ion_used_only_by_non_default_profile_is_not_default(self, service):
        entry = service.classify_ion("I1", "LM-SPCE-I1")
        assert entry["is_default"] is False

    def test_unknown_ion_is_new_unclassified(self, service):
        assert service.classify_ion("I1", "TotallyUnknownIonFF") == {
            "tier": "new_unclassified",
            "is_default": False,
        }


class TestWaterProfiles:
    def test_exactly_one_profile_is_global_default(self, service):
        defaults = [p for p in service.all_water_profiles() if p["is_default"]]
        assert len(defaults) == 1
        assert defaults[0]["profile_id"] == "SPCE_JC_LM"

    def test_each_water_type_still_exposes_its_own_local_default(self, service):
        w4_profiles = service.water_profiles_for("W4")
        local_defaults = [p for p in w4_profiles if p["is_default_for_water_type"]]
        assert len(local_defaults) == 1
        assert local_defaults[0]["profile_id"] == "OPC_LM"
        # ale globální is_default u W4 profilů je vždy False (default vody je W3)
        assert all(not p["is_default"] for p in w4_profiles)


class TestReconcile:
    def _fake_catalog(self):
        return [
            {"display_name": "FF14SB", "molecule_type": ["P"]},  # už známý
            {"display_name": "Brand New Protein FF", "molecule_type": ["P"]},
            {"display_name": "SomeNewIon", "molecule_type": ["I1"]},
            {"display_name": "Brand New Water", "molecule_type": ["W3"]},
        ]

    def test_adds_only_unknown_names_to_new_unclassified(self, service):
        added = service.reconcile(self._fake_catalog(), lambda ff: _name(ff["display_name"]))
        assert added == 3

        data = json.loads(service.path.read_text(encoding="utf-8"))
        assert data["solute_force_fields"]["P"]["new_unclassified"] == ["Brand_New_Protein_FF"]
        assert data["solvent"]["new_unclassified"]["I1"] == ["SomeNewIon"]
        assert data["solvent"]["new_unclassified"]["W3"] == ["Brand_New_Water"]

    def test_is_idempotent(self, service):
        service.reconcile(self._fake_catalog(), lambda ff: _name(ff["display_name"]))
        added_again = service.reconcile(self._fake_catalog(), lambda ff: _name(ff["display_name"]))
        assert added_again == 0

    def test_leaves_file_untouched_when_nothing_new(self, service):
        catalog = [{"display_name": "FF14SB", "molecule_type": ["P"]}]
        mtime_before = service.path.stat().st_mtime_ns
        added = service.reconcile(catalog, lambda ff: _name(ff["display_name"]))
        assert added == 0
        assert service.path.stat().st_mtime_ns == mtime_before


class TestSetSoluteTier:
    def test_moves_ff_between_tiers(self, service):
        service.set_solute_tier("P", "FF12SB", "obsolete")
        data = json.loads(service.path.read_text(encoding="utf-8"))
        assert "FF12SB" not in data["solute_force_fields"]["P"]["supported"]
        assert "FF12SB" in data["solute_force_fields"]["P"]["obsolete"]

    def test_promotes_new_unclassified_ff(self, service):
        service.reconcile(
            [{"display_name": "Brand New Protein FF", "molecule_type": ["P"]}],
            lambda ff: _name(ff["display_name"]),
        )
        service.set_solute_tier("P", "Brand_New_Protein_FF", "recommended")

        data = json.loads(service.path.read_text(encoding="utf-8"))
        assert "Brand_New_Protein_FF" in data["solute_force_fields"]["P"]["recommended"]
        assert "Brand_New_Protein_FF" not in data["solute_force_fields"]["P"]["new_unclassified"]

    def test_rejects_unknown_tier(self, service):
        with pytest.raises(ValueError):
            service.set_solute_tier("P", "FF14SB", "not_a_real_tier")

    def test_rejects_unknown_group(self, service):
        with pytest.raises(ValueError):
            service.set_solute_tier("X", "FF14SB", "recommended")
