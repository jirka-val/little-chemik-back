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

from app.core.exceptions import BadRequestError
from app.services.ghbfix import GHBFIX_FILENAME, GHBFIX_LISTOUT
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
    #
    # `thermostat`/`barostat` jsou primární volba (frontend od redesignu
    # vybírá termostat a barostat zvlášť, ensemble z nich vyplývá). Pokud
    # nepřijdou (starší klient), odvodí se ze `ensemble` přesně jako dřív:
    # NVE = žádný, NVT = Langevin, NPT = Langevin + Monte Carlo.
    ensemble: Literal["NVE", "NVT", "NPT"] = "NVT"
    thermostat: Literal["none", "langevin", "berendsen", "andersen"] | None = None
    barostat: Literal["none", "monte-carlo", "berendsen"] | None = None
    temp0: float = 298.16
    gamma_ln: float = Field(1.0, ge=0)      # Langevin (ntt=3)
    tautp: float = Field(1.0, gt=0)         # Berendsen weak coupling (ntt=1)
    vrand: int = Field(1000, gt=0)          # Andersen - kroky mezi randomizacemi (ntt=2)

    # W4.4 Pressure (jen s barostatem)
    pres0: float = 1.0
    taup: float = 2.0                       # Berendsen barostat (barostat=1)

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

    # Vybraný FF má gHBfix korekce (has_ghbfix z /api/forcefields) - mdin je
    # načte jako restrainty z ghbfix.f, viz app/services/ghbfix.py.
    ghbfix: bool = False


_NTC_NTF = {"none": 1, "h-bonds": 2, "all-bonds": 3}
_NTT = {"none": 0, "berendsen": 1, "andersen": 2, "langevin": 3}
_BAROSTAT = {"berendsen": 1, "monte-carlo": 2}


def _resolve_coupling(req: AmberMdinRequest) -> tuple[str, str]:
    """(thermostat, barostat) - explicitní volba, jinak odvozené z `ensemble`."""
    thermostat = req.thermostat
    barostat = req.barostat
    if thermostat is None:
        thermostat = "none" if req.ensemble == "NVE" else "langevin"
    if barostat is None:
        barostat = "monte-carlo" if req.ensemble == "NPT" else "none"
    if barostat != "none" and thermostat == "none":
        # Monte Carlo barostat potřebuje cílovou teplotu (temp0) pro
        # akceptační kritérium a NPH bez termostatu není podporovaný
        # produkční režim tohohle generátoru.
        raise BadRequestError(
            "A barostat requires a thermostat - pressure coupling without temperature control (NPH) is not supported.",
            code="barostat_requires_thermostat",
        )
    return thermostat, barostat


def _ensemble_label(thermostat: str, barostat: str) -> str:
    if barostat != "none":
        return "NPT"
    return "NVE" if thermostat == "none" else "NVT"


def _fmt(value: float) -> str:
    """
    Fixed-point formatting matching the AMBER mdin style in the reference
    file (e.g. "cut=10.0", "tol=0.00001") - never scientific notation
    (Fortran namelist parsers accept it, but it doesn't match house style),
    trailing zeros trimmed down to one decimal digit minimum.
    """
    text = f"{value:.8f}".rstrip("0")
    return text if not text.endswith(".") else text + "0"


def _build_title(req: AmberMdinRequest, ensemble: str) -> str:
    parts = f"Production {ensemble}, {req.duration_ns:g} ns segment, {req.dt_fs:g} fs"
    if req.hmr:
        parts += " with HMR in the topology"
    return parts


def render_amber_mdin(req: AmberMdinRequest) -> str:
    thermostat, barostat = _resolve_coupling(req)
    ensemble = _ensemble_label(thermostat, barostat)

    dt_ps = req.dt_fs / 1000.0
    nstlim = round(req.duration_ns * 1000.0 / dt_ps)
    ntpr = max(1, round(req.log_interval_ps / dt_ps))
    ntwx = max(1, round(req.traj_interval_ps / dt_ps))
    ntwr = max(1, round(req.restart_interval_ps / dt_ps))
    ntc = ntf = _NTC_NTF[req.constraints]

    has_barostat = barostat != "none"
    irest, ntx = (0, 1) if req.run_type == "new" else (1, 5)
    ntb = 2 if has_barostat else 1
    ntp = 1 if has_barostat else 0
    ntt = _NTT[thermostat]

    lines = [_build_title(req, ensemble), "&cntrl"]
    lines.append("  imin=0,")
    lines.append(f"  irest={irest}, ntx={ntx},")
    lines.append(f"  nstlim={nstlim}, dt={_fmt(dt_ps)},")
    lines.append(f"  ntb={ntb}, ntp={ntp}, cut={_fmt(req.cut)},")
    lines.append(f"  ntc={ntc}, ntf={ntf}, tol={_fmt(req.tol)},")

    # Parametry termostatu - vždy jen ty, které zvolený ntt opravdu čte
    # (design list W4.3: "gamma_ln only for ntt=3; tautp only for ntt=1.
    # Never emit both blindly.").
    if thermostat == "langevin":
        # ig=-1 (auto-seed from date/time) on every segment, including restarts -
        # Amber26.pdf 23.6.7 explicitly warns against reusing a fixed seed across
        # Langevin restarts ("synchronization" artifacts).
        lines.append(f"  ntt={ntt}, temp0={_fmt(req.temp0)}, gamma_ln={_fmt(req.gamma_ln)}, ig=-1,")
    elif thermostat == "andersen":
        # Andersen je taky stochastický - stejný důvod pro ig=-1 jako u Langevinu.
        lines.append(f"  ntt={ntt}, temp0={_fmt(req.temp0)}, vrand={req.vrand}, ig=-1,")
    elif thermostat == "berendsen":
        lines.append(f"  ntt={ntt}, temp0={_fmt(req.temp0)}, tautp={_fmt(req.tautp)},")
    else:
        lines.append(f"  ntt={ntt},")

    if has_barostat:
        # Řádek pro Monte Carlo je beze změny proti verzi před redesignem
        # (včetně taup). Design list W4.4 říká, že taup se pro MC nemusí
        # uvádět - odstranění je samostatné rozhodnutí, ne vedlejší efekt.
        lines.append(f"  barostat={_BAROSTAT[barostat]}, pres0={_fmt(req.pres0)}, taup={_fmt(req.taup)},")

    output_line = f"  ntpr={ntpr}, ntwx={ntwx}, ntwr={ntwr},"
    if req.write_energy_file:
        output_line += f" ntwe={ntpr},"
    lines.append(output_line)

    lines.append(f"  ioutfm={req.ioutfm}, ntxo={req.ntxo},")
    lines.append(f"  iwrap={req.iwrap},")
    if req.ghbfix:
        lines.append("  nmropt=1,")
    lines.append("/")
    if req.ghbfix:
        # Váha restraintů 1.0 po celý segment. Při výměnách replik (numexchg
        # > 1) by istep2 bylo nstlim*numexchg - tenhle generátor REMD nemá.
        lines += [
            "&wt",
            f"  type='REST', istep1=1, istep2={nstlim}, value1=1.0, value2=1.0,",
            "/",
            "&wt",
            "  type='END',",
            "/",
            f"LISTOUT={GHBFIX_LISTOUT}",
            f"DISANG={GHBFIX_FILENAME}",
        ]
    return "\n".join(lines) + "\n"


def mdin_overview(req: AmberMdinRequest) -> tuple[str, int]:
    """(ensemble, nstlim) pro souhrnný report - stejný výpočet jako render_amber_mdin."""
    ensemble = _ensemble_label(*_resolve_coupling(req))
    return ensemble, round(req.duration_ns * 1000.0 / (req.dt_fs / 1000.0))


def _mdin_filename(req: AmberMdinRequest) -> str:
    ensemble = _ensemble_label(*_resolve_coupling(req))
    hmr_suffix = "_hmr" if req.hmr else ""
    duration_label = f"{req.duration_ns:g}".replace(".", "p")
    return f"amber_{ensemble.lower()}_production_{duration_label}ns{hmr_suffix}.mdin"


@router.post("/{workspace_id}/amber-mdin", summary="Generuje AMBER produkční mdin (&cntrl) soubor")
async def generate_amber_mdin(workspace_id: str, req: AmberMdinRequest, preview: bool = False):
    """
    `preview=true` = živý náhled v Simulation panelu (volá se debounced při
    každé změně formuláře) - stejný výstup, jen se neloguje, ať náhled
    nezaplaví log při každém stisku klávesy.
    """
    workspace_manager.require_workspace(workspace_id)

    content = render_amber_mdin(req)
    filename = _mdin_filename(req)

    if not preview:
        thermostat, barostat = _resolve_coupling(req)
        logger.info(
            f"Generated AMBER mdin for workspace {workspace_id}: "
            f"{_ensemble_label(thermostat, barostat)} (thermostat={thermostat}, barostat={barostat}), "
            f"{req.duration_ns:g} ns, dt={req.dt_fs:g} fs, hmr={req.hmr}"
        )

    return PlainTextResponse(
        content=content,
        media_type="text/plain",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )
