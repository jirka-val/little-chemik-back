"""
Segmenty produkce a celý MD běh pro Export (app/services/md_protocol.py).
"""

import pytest
from pydantic import ValidationError

from app.api.v1.endpoints.simulation import AmberMdinRequest, _mdin_filename, mdin_overview, render_amber_mdin
from app.core.exceptions import BadRequestError
from app.services.md_protocol import render_protocol, solute_residue_count
from app.services.md_steps import ProtocolStep, default_steps

pytestmark = pytest.mark.unit


def _atom(serial, name, resn, chain, resi):
    return f"ATOM  {serial:5d} {name:<4} {resn:>3} {chain}{resi:4d}    {0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00"


_PDB = "\n".join([
    _atom(1, "N", "ALA", "A", 1), _atom(2, "CA", "ALA", "A", 1),
    _atom(3, "N", "GLY", "A", 2), _atom(4, "CA", "GLY", "A", 2),
    _atom(5, "NA", "Na+", "B", 1),
    _atom(6, "O", "WAT", "C", 1), _atom(7, "H1", "WAT", "C", 1),
])


def test_segments_split_total_length():
    req = AmberMdinRequest(duration_ns=100.0, segments=10, dt_fs=2.0)
    mdin = render_amber_mdin(req)
    assert "nstlim=5000000," in mdin  # 10 ns / 2 fs
    assert "10 x 10 ns = 100 ns" in mdin.splitlines()[0]
    assert mdin_overview(req)[1] == 5_000_000
    assert _mdin_filename(req) == "amber_nvt_production_10ns.mdin"


def test_single_segment_is_unchanged():
    req = AmberMdinRequest(duration_ns=1.0, dt_fs=2.0)
    assert "nstlim=500000," in render_amber_mdin(req)
    assert "irest=0, ntx=1," in render_amber_mdin(req)


def test_solute_ends_at_first_ion_or_water():
    assert solute_residue_count(_PDB) == 2
    assert solute_residue_count(_atom(1, "O", "WAT", "A", 1)) is None


def test_protocol_files_and_restraint_mask():
    req = AmberMdinRequest(duration_ns=100.0, segments=4, protocol=True, ensemble="NVT")
    files = dict(render_protocol(req, _PDB))
    assert list(files) == ["mini0.in", "mini1.in", "heat.in", "equil.in", "prod.in", "run_md.sh"]
    assert "restraintmask=':1-2 & !@H='" in files["mini0.in"]
    assert "ntr=" not in files["mini1.in"]
    assert "type='TEMP0'" in files["heat.in"] and "ntr=1" in files["heat.in"]
    assert "ntb=2, ntp=1" in files["equil.in"]  # NPT i při NVT produkci
    # produkce v protokolu vždy navazuje na restart
    assert "irest=1, ntx=5," in files["prod.in"]
    assert "SEGMENTS=4" in files["run_md.sh"]


def test_default_steps_are_the_former_fixed_protocol():
    assert [(s.name, s.kind) for s in default_steps()] == [
        ("mini0", "minimization"), ("mini1", "minimization"), ("heat", "heating"), ("equil", "equilibration"),
    ]


def test_custom_steps_build_files_and_run_script():
    steps = [
        ProtocolStep(name="min_h", kind="minimization", restraint="heavy", restraint_wt=1000.0,
                     maxcyc=2000, ncyc=250, extra="drms=0.01, dx0=0.01"),
        ProtocolStep(name="min_solv", kind="minimization", restraint="solute", maxcyc=3000),
        ProtocolStep(name="min_all", kind="minimization", maxcyc=4000),
        ProtocolStep(name="eq_nvt", kind="equilibration", restraint="custom", restraint_mask=":1-2@P",
                     length_ps=50.0, pressure=False),
        ProtocolStep(name="eq_npt", kind="equilibration", length_ps=500.0),
    ]
    files = dict(render_protocol(AmberMdinRequest(protocol=True, steps=steps), _PDB))

    assert list(files) == ["min_h.in", "min_solv.in", "min_all.in", "eq_nvt.in", "eq_npt.in", "prod.in", "run_md.sh"]
    assert "maxcyc=2000, ncyc=250," in files["min_h.in"]
    assert "restraint_wt=1000.0, restraintmask='!@H='," in files["min_h.in"]
    assert "drms=0.01, dx0=0.01," in files["min_h.in"]
    assert "restraintmask=':1-2 & !@H='" in files["min_solv.in"]
    # ekvilibrace po minimalizaci začíná znovu, po MD kroku navazuje
    assert "irest=0, ntx=1, tempi=298.16," in files["eq_nvt.in"] and "ntb=1, ntp=0" in files["eq_nvt.in"]
    assert "restraintmask=':1-2@P'" in files["eq_nvt.in"]
    assert "irest=1, ntx=5," in files["eq_npt.in"] and "ntb=2, ntp=1" in files["eq_npt.in"]
    assert "nstlim=250000," in files["eq_npt.in"]  # 500 ps / 2 fs
    script = files["run_md.sh"]
    assert 'run_step min_h min_h.in "$CRD" "$CRD"' in script
    assert "run_step min_all min_all.in min_solv.rst7\n" in script
    assert "run_step eq_npt eq_npt.in eq_nvt.rst7\n" in script
    assert "prev=eq_npt" in script


def test_heating_ramps_from_start_temperature():
    steps = [ProtocolStep(name="heat_a", kind="heating", temp_start=100.0, length_ps=20.0)]
    heat = dict(render_protocol(AmberMdinRequest(protocol=True, steps=steps, temp0=300.0), None))["heat_a.in"]
    assert "tempi=100.0, temp0=300.0" in heat
    assert "value1=100.0, value2=300.0," in heat
    assert "ntr=" not in heat


def test_chained_heating_to_own_target_temperatures():
    steps = [
        ProtocolStep(name="heat1", kind="heating", temp_start=0.0, temp_target=100.0, length_ps=50.0),
        ProtocolStep(name="heat2", kind="heating", temp_start=100.0, temp_target=250.0, length_ps=200.0),
        ProtocolStep(name="eq_cold", kind="equilibration", temp_target=250.0, length_ps=100.0, pressure=False),
        ProtocolStep(name="heat3", kind="heating", temp_start=250.0, length_ps=100.0),
    ]
    files = dict(render_protocol(AmberMdinRequest(protocol=True, steps=steps, temp0=300.0), None))

    heat1, heat2, heat3 = files["heat1.in"], files["heat2.in"], files["heat3.in"]
    assert "0 -> 100 K" in heat1.splitlines()[0]
    assert "irest=0, ntx=1," in heat1 and "tempi=0.0, temp0=100.0" in heat1
    assert "value1=0.0, value2=100.0," in heat1
    # druhý ohřev navazuje na rychlosti prvního, tempi nemá
    assert "irest=1, ntx=5," in heat2 and "tempi=" not in heat2
    assert "value1=100.0, value2=250.0," in heat2 and "temp0=250.0" in heat2
    eq = files["eq_cold.in"]
    assert "temp0=250.0" in eq and "irest=1, ntx=5," in eq
    # bez cílové teploty ohřev končí na teplotě produkce
    assert "value1=250.0, value2=300.0," in heat3
    assert "temp0=300.0" in files["prod.in"]


def test_report_describes_step_temperatures():
    from app.services.md_steps import describe_step
    heat = ProtocolStep(name="h", kind="heating", temp_start=10.0, temp_target=100.0, length_ps=50.0)
    assert describe_step(heat, 300.0).startswith("heating from 10 K to 100 K over 50 ps")
    eq = ProtocolStep(name="e", kind="equilibration", length_ps=2000.0)
    assert describe_step(eq, 300.0).startswith("NPT equilibration at 300 K for 2 ns")


@pytest.mark.parametrize("steps, message", [
    ([], "at least one"),
    ([{"name": "a", "kind": "minimization"}], "last preparation step"),
    ([{"name": "a", "kind": "heating"}, {"name": "A", "kind": "equilibration"}], "used twice"),
    ([{"name": "prod_x", "kind": "heating"}], "reserved"),
    ([{"name": "a b", "kind": "heating"}], "pattern"),
    ([{"name": "a", "kind": "minimization", "maxcyc": 10, "ncyc": 20}, {"name": "b", "kind": "heating"}], "ncyc"),
    ([{"name": "a", "kind": "heating", "restraint": "custom"}], "needs a mask"),
    ([{"name": "a", "kind": "heating", "extra": "nonsense"}], "key=value"),
])
def test_invalid_step_lists_are_rejected(steps, message):
    with pytest.raises(ValidationError, match=message):
        AmberMdinRequest(protocol=True, steps=steps)


@pytest.mark.parametrize("field, value", [
    ("temp_target", -5.0), ("temp_target", 5000.0), ("temp_target", float("inf")), ("temp_target", float("nan")),
    ("temp_start", float("inf")), ("temp_start", 2000.0),
    ("length_ps", 0.004),  # pár kroků integrace - dřív rozbitý &wt blok
    ("restraint_wt", float("inf")),
])
def test_invalid_step_values_are_rejected(field, value):
    with pytest.raises(ValidationError, match=field):
        AmberMdinRequest(protocol=True, steps=[{"name": "h", "kind": "heating", field: value}])


@pytest.mark.parametrize("payload", [
    {"temp0": float("nan")}, {"temp0": -1.0}, {"temp0": 5000.0}, {"duration_ns": float("inf")},
    {"duration_ns": 1e-6},  # nstlim=0
    {"duration_ns": 0.002, "segments": 10},
    {"cut": 0.0}, {"tol": -1.0}, {"pres0": 0.0}, {"dt_fs": 50.0},
])
def test_invalid_run_values_are_rejected(payload):
    with pytest.raises(ValidationError):
        AmberMdinRequest(**payload)


def test_report_survives_steps_saved_before_current_limits():
    from app.services.md_steps import steps_from_settings
    steps = steps_from_settings({"steps": [{"name": "h", "kind": "heating", "length_ps": 0.5}]})
    assert steps[0].length_ps == 0.5 and steps[0].temp_target is None


def test_extra_cannot_override_step_settings():
    steps = [ProtocolStep(name="h", kind="heating", extra="nstlim=5")]
    with pytest.raises(BadRequestError, match="nstlim"):
        render_protocol(AmberMdinRequest(protocol=True, steps=steps), None)


def test_protocol_ghbfix_in_md_steps_only():
    req = AmberMdinRequest(protocol=True, ghbfix=True)
    files = dict(render_protocol(req, None))
    for name in ("heat.in", "equil.in", "prod.in"):
        assert "DISANG=ghbfix.f" in files[name]
    assert "DISANG" not in files["mini0.in"]
    # bez PDB aspoň všechno kromě vody
    assert "restraintmask='!:WAT,HOH,SOL & !@H='" in files["mini0.in"]
