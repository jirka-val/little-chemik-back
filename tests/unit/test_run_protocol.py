"""
Protokol celé přípravy (app/services/run_protocol.py).
"""

import pytest

from app.core.exceptions import BadRequestError
from app.services.run_protocol import (
    FORMAT_VERSION, RunProtocol, parse_protocol, protocol_filename, protocol_fingerprint,
)

pytestmark = pytest.mark.unit


def _raw(**extra):
    return {"format": "forge-protocol", "format_version": FORMAT_VERSION, "name": "RNA standard", **extra}


def test_minimal_protocol_gets_defaults():
    protocol, warnings = parse_protocol(_raw())
    assert warnings == []
    assert protocol.mode == "standard"
    assert protocol.preparation.ph == 7.0 and protocol.preparation.existing_solvent == "keep"
    assert protocol.simulation.temp0 == 298.16
    assert protocol.export.top is True


def test_full_protocol_round_trips():
    raw = _raw(
        mode="expert",
        forcefields={"R": {"id": 12, "name": "OL3"}, "W": {"id": 3, "name": "OPC"}},
        preparation={"ph": 6.5, "add_solvent": True, "box_shape": "octahedron",
                     "additional_salts": [{"positive_ion": "K+", "negative_ion": "Cl-", "ionic_strength": 0.1}]},
        protonation=[{"chain": "A", "resseq": 12, "state": "HIE", "resname": "HIS"}],
        simulation={"protocol": True, "steps": [{"name": "heat1", "kind": "heating", "temp_target": 100}]},
    )
    protocol, _ = parse_protocol(raw)
    again, _ = parse_protocol(protocol.model_dump(mode="json"))
    assert again == protocol
    assert again.simulation.steps[0].temp_target == 100


@pytest.mark.parametrize("raw, code", [
    ([], "protocol_invalid"),
    ({"format": "something-else"}, "protocol_invalid"),
    ({"format": "forge-protocol", "name": "x"}, "protocol_invalid"),  # bez verze
    (_raw(format_version=FORMAT_VERSION + 1), "protocol_too_new"),
])
def test_foreign_or_newer_files_are_rejected(raw, code):
    with pytest.raises(BadRequestError) as e:
        parse_protocol(raw)
    assert e.value.code == code


@pytest.mark.parametrize("extra, where", [
    ({"preparation": {"ph": 20}}, "preparation.ph"),
    ({"simulation": {"temp0": -5}}, "simulation.temp0"),
    ({"simulation": {"steps": [{"name": "a", "kind": "minimization"}]}}, "simulation.steps"),
    ({"name": ""}, "name"),
    ({"mode": "wizard"}, "mode"),
])
def test_invalid_values_name_the_setting(extra, where):
    with pytest.raises(BadRequestError, match=f"^The protocol is not valid - {where}"):
        parse_protocol(_raw(**extra))


def test_unknown_settings_are_only_warnings():
    _, warnings = parse_protocol(_raw(color="blue", preparation={"ph": 7, "magic": 1}))
    assert warnings == ["Unknown setting 'color' was ignored.", "Unknown setting 'preparation.magic' was ignored."]


def test_fingerprint_ignores_name_but_not_settings():
    a = RunProtocol(name="A")
    assert protocol_fingerprint(a) == protocol_fingerprint(RunProtocol(name="B", notes="x"))
    changed = RunProtocol(name="A", preparation={"ph": 6.0})
    assert protocol_fingerprint(a) != protocol_fingerprint(changed)


def test_filename_is_safe():
    assert protocol_filename(RunProtocol(name="RNA / 1 µs, OL3")) == "RNA___1__s__OL3.forge-protocol.json"
    assert protocol_filename(RunProtocol(name="***")) == "protocol.forge-protocol.json"
