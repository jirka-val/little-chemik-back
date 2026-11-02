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


# ---------------------------------------------------------------------------
# Teplotní průběh ohřevu (temp_points)
# ---------------------------------------------------------------------------

def _schedule(*points):
    return [{"time_ps": t, "temp": temp} for t, temp in points]


def test_heating_schedule_renders_one_wt_block_per_segment():
    # průběh z referenčního vstupu skupiny: 10 -> 50 K, 50 -> 298.16 K, držení; 100 ps po 2 fs
    step = ProtocolStep(name="heat1", kind="heating", length_ps=100.0,
                        temp_points=_schedule((0, 10.0), (20, 50.0), (90, 298.16), (100, 298.16)))
    assert (step.temp_start, step.temp_target) == (10.0, 298.16)
    heat = dict(render_protocol(AmberMdinRequest(protocol=True, steps=[step], temp0=298.16), None))["heat1.in"]
    assert "nstlim=50000," in heat
    assert "tempi=10.0, temp0=298.16" in heat
    assert "type='TEMP0', istep1=0, istep2=10000, value1=10.0, value2=50.0," in heat
    assert "type='TEMP0', istep1=10001, istep2=45000, value1=50.0, value2=298.16," in heat
    assert "type='TEMP0', istep1=45001, istep2=50000, value1=298.16, value2=298.16," in heat
    assert heat.count("type='TEMP0'") == 3
    assert "3 temperature segments" in heat.splitlines()[0]


def test_heating_without_schedule_is_unchanged():
    step = ProtocolStep(name="heat", kind="heating", temp_start=0.0, length_ps=100.0)
    heat = dict(render_protocol(AmberMdinRequest(protocol=True, steps=[step], temp0=300.0), None))["heat.in"]
    assert "istep1=0, istep2=40000, value1=0.0, value2=300.0," in heat
    assert "istep1=40001, istep2=50000, value1=300.0, value2=300.0," in heat


def test_heating_schedule_with_ghbfix_keeps_the_rest_block():
    step = ProtocolStep(name="heat1", kind="heating", length_ps=10.0,
                        temp_points=_schedule((0, 0.0), (5, 300.0), (10, 300.0)))
    heat = dict(render_protocol(AmberMdinRequest(protocol=True, steps=[step], temp0=300.0, ghbfix=True), None))["heat1.in"]
    assert heat.count("type='TEMP0'") == 2 and "type='REST'" in heat and "DISANG=" in heat


@pytest.mark.parametrize("points, length, message", [
    (((0, 10.0),), 100.0, "at least 2 points"),
    (((5, 10.0), (100, 300.0)), 100.0, "start at 0 ps"),
    (((0, 10.0), (50, 100.0), (40, 200.0), (100, 300.0)), 100.0, "time order"),
    (((0, 10.0), (50, 300.0)), 100.0, "end at the step length"),
])
def test_invalid_heating_schedules(points, length, message):
    with pytest.raises(ValidationError, match=message):
        ProtocolStep(name="heat1", kind="heating", length_ps=length, temp_points=_schedule(*points))


def test_schedule_only_for_heating():
    with pytest.raises(ValidationError, match="only for heating"):
        ProtocolStep(name="eq", kind="equilibration", length_ps=10.0, temp_points=_schedule((0, 300.0), (10, 300.0)))


def test_schedule_is_described_in_the_report():
    from app.services.md_steps import describe_step
    step = ProtocolStep(name="heat1", kind="heating", length_ps=100.0,
                        temp_points=_schedule((0, 10.0), (20, 50.0), (90, 298.16), (100, 298.16)))
    text = describe_step(step, 298.16)
    assert "from 10 K to 298.16 K over 100 ps" in text and "50 K at 20 ps" in text


# ---------------------------------------------------------------------------
# Ohřev v NPT a intervaly výpisu/trajektorie v ps (referenční vstupy skupiny)
# ---------------------------------------------------------------------------

def _files(steps, **req):
    return dict(render_protocol(AmberMdinRequest(protocol=True, steps=steps, **req), None))


def test_heating_in_npt_uses_the_production_barostat():
    step = ProtocolStep(name="heat", kind="heating", heat_pressure=True, length_ps=500.0,
                        temp_points=_schedule((0, 10.0), (400, 298.16), (500, 298.16)))
    heat = _files([step], temp0=298.16, ensemble="NPT", barostat="monte-carlo", pres0=1.0, taup=2.0)["heat.in"]
    assert "ntb=2, ntp=1, cut=10.0, barostat=2, pres0=1.0, taup=2.0," in heat
    assert "NPT" in heat.splitlines()[0] and "nstlim=250000," in heat
    # bez barostatu v produkci Monte Carlo (jako u ekvilibrace)
    heat = _files([step], temp0=298.16, ensemble="NVT")["heat.in"]
    assert "barostat=2," in heat


def test_heating_is_nvt_by_default():
    heat = _files([ProtocolStep(name="heat", kind="heating", length_ps=100.0)], temp0=300.0)["heat.in"]
    assert "ntb=1, ntp=0," in heat and "barostat" not in heat


def test_output_intervals_in_ps():
    steps = [
        ProtocolStep(name="heat", kind="heating", length_ps=100.0, log_interval_ps=1.0, traj_interval_ps=5.0),
        ProtocolStep(name="eq", kind="equilibration", length_ps=100.0, pressure=False,
                     log_interval_ps=2.0, traj_interval_ps=10.0),
    ]
    files = _files(steps, temp0=300.0, dt_fs=2.0)
    assert "ntpr=500, ntwx=2500," in files["heat.in"]          # 1 ps a 5 ps po 2 fs
    assert "ntpr=1000, ntwx=5000," in files["eq.in"]
    # bez vlastních intervalů: ohřev 1 % kroku, ekvilibrace podle produkce
    files = _files([ProtocolStep(name="heat", kind="heating", length_ps=100.0),
                    ProtocolStep(name="eq", kind="equilibration", length_ps=100.0)],
                   temp0=300.0, log_interval_ps=4.0, traj_interval_ps=8.0)
    assert "ntpr=500, ntwx=500," in files["heat.in"]
    assert "ntpr=2000, ntwx=4000," in files["eq.in"]


def test_output_interval_longer_than_the_step_is_rejected():
    with pytest.raises(ValidationError, match="log interval"):
        ProtocolStep(name="eq", kind="equilibration", length_ps=10.0, log_interval_ps=20.0)


def test_minimization_prints_every_10_cycles():
    files = _files([ProtocolStep(name="mini", kind="minimization"),
                    ProtocolStep(name="heat", kind="heating", length_ps=10.0)], temp0=300.0)
    assert "ntpr=10," in files["mini.in"]
