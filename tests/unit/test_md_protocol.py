"""
Segmenty produkce a celý MD běh pro Export (app/services/md_protocol.py).
"""

import pytest

from app.api.v1.endpoints.simulation import AmberMdinRequest, _mdin_filename, mdin_overview, render_amber_mdin
from app.services.md_protocol import render_protocol, solute_residue_count

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


def test_protocol_ghbfix_in_md_steps_only():
    req = AmberMdinRequest(protocol=True, ghbfix=True)
    files = dict(render_protocol(req, None))
    for name in ("heat.in", "equil.in", "prod.in"):
        assert "DISANG=ghbfix.f" in files[name]
    assert "DISANG" not in files["mini0.in"]
    # bez PDB aspoň všechno kromě vody
    assert "restraintmask='!:WAT,HOH,SOL & !@H='" in files["mini0.in"]
