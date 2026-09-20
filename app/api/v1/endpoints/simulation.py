# app/api/v1/endpoints/simulation.py
"""
Generuje AMBER `&cntrl` mdin soubor pro jeden produkční MD segment (sander/
pmemd). Vychází z FORGE_general_design_v5.xlsx (list "MD input mapping",
sekce W4.1-W4.7) a z Amber26.pdf sekce 23.6 (General minimization and
dynamics parameters).

ROZSAH (v1): jeden segment na jeden request - žádné automatické dělení
dlouhého běhu na N+1 segmentů s restart orchestrací (to je podle designového
listu samo o sobě označené jako "future extension", ne součást dnešního
rozsahu). Uživatel si segmenty spouští ručně, po jednom, s run_type
"new"/"continue".

Generovaný text je čistě funkcí vstupních parametrů formuláře - nezávisí na
konkrétní struktuře/topologii daného workspace (mdin soubor neobsahuje žádné
atomové ani rezidové informace). workspace_id v cestě je jen kvůli
konzistenci s zbytkem API (a ověření, že workspace existuje), ne proto, že by
se z něj něco čtlo.
"""
import logging
from typing import Literal

from fastapi import APIRouter
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)
router = APIRouter()


class AmberMdinRequest(BaseModel):
    # W4.2 Run structure
    run_type: Literal["new", "continue"] = "new"
    duration_ns: float = Field(100.0, gt=0)
    dt_fs: float = Field(2.0, gt=0)
    hmr: bool = False

    # W4.3 Ensemble / temperature
    ensemble: Literal["NVE", "NVT", "NPT"] = "NVT"
    temp0: float = 298.16
    gamma_ln: float = Field(1.0, ge=0)

    # W4.4 Pressure (NPT only)
    pres0: float = 1.0
    taup: float = 2.0

    # W4.5 Constraints
    constraints: Literal["none", "h-bonds", "all-bonds"] = "h-bonds"
    tol: float = 0.00001

    # W4.6 PBC / nonbonded
    cut: float = 10.0

    # W4.2 Output (physical intervals in ps - converted to step counts below)
    log_interval_ps: float = Field(10.0, gt=0)
    traj_interval_ps: float = Field(10.0, gt=0)
    restart_interval_ps: float = Field(1000.0, gt=0)
    write_energy_file: bool = False

    # W4.6/4.7 format & expert
    ioutfm: Literal[0, 1] = 1
    ntxo: Literal[1, 2] = 2
    iwrap: Literal[0, 1] = 1


_NTC_NTF = {"none": 1, "h-bonds": 2, "all-bonds": 3}


def _fmt(value: float) -> str:
    """
    Fixed-point formatting matching the AMBER mdin style in the reference
    file (e.g. "cut=10.0", "tol=0.00001") - never scientific notation
    (Fortran namelist parsers accept it, but it doesn't match house style),
    trailing zeros trimmed down to one decimal digit minimum.
    """
    text = f"{value:.8f}".rstrip("0")
    return text if not text.endswith(".") else text + "0"


def _build_title(req: AmberMdinRequest) -> str:
    parts = f"Production {req.ensemble}, {req.duration_ns:g} ns segment, {req.dt_fs:g} fs"
    if req.hmr:
        parts += " with HMR in the topology"
    return parts


def render_amber_mdin(req: AmberMdinRequest) -> str:
    dt_ps = req.dt_fs / 1000.0
    nstlim = round(req.duration_ns * 1000.0 / dt_ps)
    ntpr = max(1, round(req.log_interval_ps / dt_ps))
    ntwx = max(1, round(req.traj_interval_ps / dt_ps))
    ntwr = max(1, round(req.restart_interval_ps / dt_ps))
    ntc = ntf = _NTC_NTF[req.constraints]

    is_npt = req.ensemble == "NPT"
    is_nve = req.ensemble == "NVE"
    irest, ntx = (0, 1) if req.run_type == "new" else (1, 5)
    ntb = 2 if is_npt else 1
    ntp = 1 if is_npt else 0
    ntt = 0 if is_nve else 3

    lines = [_build_title(req), "&cntrl"]
    lines.append("  imin=0,")
    lines.append(f"  irest={irest}, ntx={ntx},")
    lines.append(f"  nstlim={nstlim}, dt={_fmt(dt_ps)},")
    lines.append(f"  ntb={ntb}, ntp={ntp}, cut={_fmt(req.cut)},")
    lines.append(f"  ntc={ntc}, ntf={ntf}, tol={_fmt(req.tol)},")

    if not is_nve:
        # ig=-1 (auto-seed from date/time) on every segment, including restarts -
        # Amber26.pdf 23.6.7 explicitly warns against reusing a fixed seed across
        # Langevin restarts ("synchronization" artifacts).
        lines.append(f"  ntt={ntt}, temp0={_fmt(req.temp0)}, gamma_ln={_fmt(req.gamma_ln)}, ig=-1,")
    else:
        lines.append(f"  ntt={ntt},")

    if is_npt:
        lines.append(f"  barostat=2, pres0={_fmt(req.pres0)}, taup={_fmt(req.taup)},")

    output_line = f"  ntpr={ntpr}, ntwx={ntwx}, ntwr={ntwr},"
    if req.write_energy_file:
        output_line += f" ntwe={ntpr},"
    lines.append(output_line)

    lines.append(f"  ioutfm={req.ioutfm}, ntxo={req.ntxo},")
    lines.append(f"  iwrap={req.iwrap},")
    lines.append("/")
    return "\n".join(lines) + "\n"


def _mdin_filename(req: AmberMdinRequest) -> str:
    hmr_suffix = "_hmr" if req.hmr else ""
    duration_label = f"{req.duration_ns:g}".replace(".", "p")
    return f"amber_{req.ensemble.lower()}_production_{duration_label}ns{hmr_suffix}.mdin"


@router.post("/{workspace_id}/amber-mdin", summary="Generuje AMBER produkční mdin (&cntrl) soubor")
async def generate_amber_mdin(workspace_id: str, req: AmberMdinRequest):
    workspace_manager.require_workspace(workspace_id)

    content = render_amber_mdin(req)
    filename = _mdin_filename(req)

    logger.info(
        f"Generated AMBER mdin for workspace {workspace_id}: "
        f"{req.ensemble}, {req.duration_ns:g} ns, dt={req.dt_fs:g} fs, hmr={req.hmr}"
    )

    return PlainTextResponse(
        content=content,
        media_type="text/plain",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )
