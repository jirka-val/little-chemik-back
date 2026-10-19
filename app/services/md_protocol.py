"""
Celý MD běh pro Export (AmberMdinRequest.protocol): minimalizace, ohřev,
ekvilibrace, produkce rozdělená na segmenty a master skript run_md.sh, který
kroky spustí za sebou.

    mini0  minimalizace s restrainty na solut (relaxuje voda a ionty)
    mini1  minimalizace bez restraintů
    heat   ohřev 0 K -> temp0 v NVT, solut držený restrainty
    equil  ekvilibrace v NPT (hustota/box), bez restraintů
    prod_001 .. prod_N  produkční segmenty, každý navazuje na restart předchozího

Délky a síly restraintů níže jsou běžné výchozí hodnoty, ne převzaté z
konkrétního protokolu skupiny - až budou k dispozici referenční vstupy
(mini0.in, mini1.in ...), srovnat a upravit.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

from app.api.v1.endpoints.simulation import (
    AmberMdinRequest, _NTC_NTF, _fmt, _mdin_filename, ghbfix_wt_lines, render_amber_mdin,
)

RUN_SCRIPT = "run_md.sh"
TOPOLOGY_FILE = "structure.prmtop"
COORDINATES_FILE = "structure.crd"

MINI0_CYCLES = 5000
MINI1_CYCLES = 10000
MINI_RESTRAINT_WT = 10.0   # kcal/mol/A^2
HEAT_PS = 100.0
HEAT_DT_FS = 2.0           # ohřev nejvýš s 2 fs, i když produkce jede s HMR na 4 fs
HEAT_RESTRAINT_WT = 5.0
EQUIL_NS = 1.0

_SOLVENT_RESNAMES = {"WAT", "HOH", "SOL"}


def solute_residue_count(pdb_text: str) -> Optional[int]:
    """
    Počet reziduí solutu. Připravené PDB je seřazené solut -> ionty -> voda
    (PDBService.reorder_and_inject_eps), solut tedy končí prvním vodním nebo
    jednoatomovým (iont) reziduem. None, když PDB nemá žádný solut.
    """
    residues: List[Tuple[str, int]] = []  # (resname, počet atomů)
    last_key = None
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21], line[22:27])
        if key != last_key:
            residues.append((line[17:21].strip(), 0))
            last_key = key
        name, count = residues[-1]
        residues[-1] = (name, count + 1)
    solute = 0
    for name, count in residues:
        if name in _SOLVENT_RESNAMES or count == 1:
            break
        solute += 1
    return solute or None


def _restraint_mask(solute_residues: Optional[int]) -> str:
    # těžké atomy solutu; bez PDB aspoň všechno kromě vody
    if solute_residues:
        return f":1-{solute_residues} & !@H="
    return "!:WAT,HOH,SOL & !@H="


def _minimization(title: str, cycles: int, req: AmberMdinRequest, mask: Optional[str]) -> str:
    lines = [title, "&cntrl",
             "  imin=1,",
             f"  maxcyc={cycles}, ncyc={cycles // 2},",
             f"  ntb=1, cut={_fmt(req.cut)},",
             "  ntpr=100,"]
    if mask:
        lines.append(f"  ntr=1, restraint_wt={_fmt(MINI_RESTRAINT_WT)}, restraintmask='{mask}',")
    lines.append("/")
    return "\n".join(lines) + "\n"


def _heating(req: AmberMdinRequest, mask: str) -> str:
    dt_fs = min(req.dt_fs, HEAT_DT_FS)
    dt_ps = dt_fs / 1000.0
    nstlim = round(HEAT_PS / dt_ps)
    ntc = _NTC_NTF[req.constraints]
    lines = [f"Heating 0 -> {req.temp0:g} K, NVT, {HEAT_PS:g} ps, restrained solute", "&cntrl",
             "  imin=0, irest=0, ntx=1,",
             f"  nstlim={nstlim}, dt={_fmt(dt_ps)},",
             f"  ntb=1, ntp=0, cut={_fmt(req.cut)},",
             f"  ntc={ntc}, ntf={ntc}, tol={_fmt(req.tol)},",
             f"  ntt=3, tempi=0.0, temp0={_fmt(req.temp0)}, gamma_ln={_fmt(max(req.gamma_ln, 1.0))}, ig=-1,",
             f"  ntr=1, restraint_wt={_fmt(HEAT_RESTRAINT_WT)}, restraintmask='{mask}',",
             f"  ntpr={max(1, nstlim // 100)}, ntwx={max(1, nstlim // 100)}, ntwr={nstlim},",
             f"  ioutfm={req.ioutfm}, ntxo={req.ntxo},",
             "  nmropt=1,",
             "/"]
    # teplota roste lineárně přes prvních 80 % ohřevu, zbytek drží temp0
    ramp_end = round(nstlim * 0.8)
    ramp = ("&wt", f"  type='TEMP0', istep1=0, istep2={ramp_end}, value1=0.0, value2={_fmt(req.temp0)},", "/",
            "&wt", f"  type='TEMP0', istep1={ramp_end + 1}, istep2={nstlim}, value1={_fmt(req.temp0)}, "
                   f"value2={_fmt(req.temp0)},", "/")
    if req.ghbfix:
        lines += ghbfix_wt_lines(nstlim, ramp)
    else:
        lines += [*ramp, "&wt", "  type='END',", "/"]
    return "\n".join(lines) + "\n"


def _equilibration(req: AmberMdinRequest) -> str:
    # NPT s Monte Carlo barostatem, aby se ustálila hustota - i když produkce
    # běží v NVT. Termostat Langevin, pokud produkce žádný nemá.
    thermostat = req.thermostat if req.thermostat not in (None, "none") else "langevin"
    equil = req.model_copy(update={
        "duration_ns": EQUIL_NS, "segments": 1, "protocol": True,
        "ensemble": "NPT", "thermostat": thermostat,
        "barostat": req.barostat if req.barostat not in (None, "none") else "monte-carlo",
    })
    text = render_amber_mdin(equil)
    first, rest = text.split("\n", 1)
    return f"Equilibration NPT, {EQUIL_NS:g} ns\n{rest}"


def _run_script(req: AmberMdinRequest) -> str:
    segments = req.segments
    return f"""#!/usr/bin/env bash
# FORGE MD run: mini0 -> mini1 -> heat -> equil -> prod_001..prod_{segments:03d}
# Production: {req.duration_ns:g} ns = {segments} x {req.duration_ns / segments:g} ns segments.
#
# Usage:   ./{RUN_SCRIPT}                     (GPU, pmemd.cuda)
#          MD_ENGINE=pmemd.MPI ./{RUN_SCRIPT}
# The script can be re-run after an interruption: finished steps (marked by
# a .done file) are skipped and it continues with the first unfinished one.
set -euo pipefail
cd "$(dirname "$0")"

MD_ENGINE="${{MD_ENGINE:-pmemd.cuda}}"
TOP={TOPOLOGY_FILE}
CRD={COORDINATES_FILE}
SEGMENTS={segments}

# run_step NAME INPUT START_COORDS [REFERENCE_COORDS]
run_step() {{
    local name=$1 input=$2 start=$3 ref=${{4:-}}
    if [[ -f $name.done ]]; then
        echo "== $name: already done, skipping"
        return
    fi
    echo "== $name ($(date '+%Y-%m-%d %H:%M:%S'))"
    $MD_ENGINE -O -i "$input" -o "$name.out" -p "$TOP" -c "$start" \\
        -r "$name.rst7" -x "$name.nc" -inf "$name.mdinfo" ${{ref:+-ref "$ref"}}
    touch "$name.done"
}}

run_step mini0 mini0.in "$CRD" "$CRD"
run_step mini1 mini1.in mini0.rst7
run_step heat  heat.in  mini1.rst7 "$CRD"
run_step equil equil.in heat.rst7

prev=equil
for i in $(seq 1 "$SEGMENTS"); do
    name=$(printf "prod_%03d" "$i")
    run_step "$name" prod.in "$prev.rst7"
    prev=$name
done
echo "== all steps finished"
"""


def render_protocol(req: AmberMdinRequest, pdb_text: Optional[str]) -> List[Tuple[str, str]]:
    """(jméno souboru, obsah) všech vstupů běhu včetně run_md.sh."""
    prod = req.model_copy(update={"protocol": True})
    mask = _restraint_mask(solute_residue_count(pdb_text) if pdb_text else None)
    return [
        ("mini0.in", _minimization("Minimization 0: solvent and ions, solute restrained",
                                   MINI0_CYCLES, req, mask)),
        ("mini1.in", _minimization("Minimization 1: whole system", MINI1_CYCLES, req, None)),
        ("heat.in", _heating(req, mask)),
        ("equil.in", _equilibration(req)),
        (_mdin_filename(prod), render_amber_mdin(prod)),
        (RUN_SCRIPT, _run_script(req)),
    ]

