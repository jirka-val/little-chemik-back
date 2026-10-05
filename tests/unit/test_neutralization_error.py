"""
Náboj, který neutralizační iont první soli nevyrovná přesně (lichý náboj a
Mg2+), je chyba nastavení pro uživatele (400), ne pád aplikace (500).
"""

import pytest

from app.services.structure.forge_service import ForgeNeutralizationError, neutralization_error_from

pytestmark = pytest.mark.unit


def test_builder_message_becomes_user_error():
    error = neutralization_error_from(ValueError("System charge -57 cannot be exactly neutralized by Mg2+ (+2)"))
    assert isinstance(error, ForgeNeutralizationError)
    assert error.status_code == 400
    assert error.code == "neutralization_impossible"
    assert error.payload == {"system_charge": -57, "ion": "Mg2+", "ion_charge": 2}
    assert "first salt" in error.message


def test_other_value_errors_are_left_alone():
    assert neutralization_error_from(ValueError("something else went wrong")) is None
