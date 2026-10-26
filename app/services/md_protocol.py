"""
Celý MD běh pro Export (AmberMdinRequest.protocol): přípravné kroky ze
seznamu `steps` (app/services/md_steps.py), produkce rozdělená na segmenty a
master skript run_md.sh, který kroky spustí za sebou. Výchozí kroky:

    mini0  minimalizace s restrainty na solut (relaxuje voda a ionty)
    mini1  minimalizace bez restraintů
    heat   ohřev 0 K -> temp0 v NVT, solut držený restrainty (konečná teplota
           jde nastavit, takže víc ohřevů za sebou dá libovolný průběh)
    equil  ekvilibrace v NPT (hustota/box), bez restraintů
    prod_001 .. prod_N  produkční segmenty, každý navazuje na restart předchozího

Každý krok startuje z restartu předchozího (první ze structure.crd), kroky s
restrainty mají jako referenci structure.crd.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

from app.api.v1.endpoints.simulation import (
    AmberMdinRequest, _NTC_NTF, _fmt, _mdin_filename, ghbfix_wt_lines, render_amber_mdin,
)
from app.core.exceptions import BadRequestError
from app.workspaces.manager import workspace_manager
from app.services.md_steps import (
    ProtocolStep, check_steps, default_steps, extra_items, is_md, step_temperature,
)

RUN_SCRIPT = "run_md.sh"
TOPOLOGY_FILE = "structure.prmtop"
COORDINATES_FILE = "structure.crd"
STRUCTURE_FILE = "structure.pdb"     # z něj maska solutu pro restrainty

HEAT_DT_FS = 2.0           # ohřev nejvýš s 2 fs, i když produkce jede s HMR na 4 fs
HEAVY_ATOMS_MASK = "!@H="

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


def _step_mask(step: ProtocolStep, solute_mask: str) -> Optional[str]:
    if step.restraint == "solute":
        return solute_mask
    if step.restraint == "heavy":
        return HEAVY_ATOMS_MASK
    if step.restraint == "custom":
        return step.restraint_mask.strip()
    return None


def _restraint_line(step: ProtocolStep, mask: Optional[str]) -> List[str]:
    if not mask:
        return []
    return [f"  ntr=1, restraint_wt={_fmt(step.restraint_wt)}, restraintmask='{mask}',"]


def _with_extra(text: str, step: ProtocolStep) -> str:
    """Doplní `extra` položky na konec &cntrl; klíč, který krok už nastavuje, je chyba."""
    items = extra_items(step.extra, step.name)
    if not items:
        return text
    lines = text.split("\n")
    end = lines.index("/")  # konec &cntrl (za ním můžou být &wt bloky)
    used = {part.split("=")[0].strip().lower()
            for line in lines[2:end] for part in line.split(",") if "=" in part}
    clash = sorted({key for key, _ in items} & used)
    if clash:
        raise BadRequestError(
            f"Step {step.name}: {', '.join(clash)} is already set by the step - change it in the step settings instead.",
            code="protocol_step_extra_clash",
        )
    lines[end:end] = ["  " + ", ".join(f"{k}={v}" for k, v in items) + ","]
    return "\n".join(lines)


def _minimization(step: ProtocolStep, req: AmberMdinRequest, mask: Optional[str]) -> str:
    ncyc = step.ncyc if step.ncyc is not None else step.maxcyc // 2
    title = f"Minimization {step.name}: {step.maxcyc} cycles" + (", restrained" if mask else "")
    lines = [title, "&cntrl",
             "  imin=1,",
             f"  maxcyc={step.maxcyc}, ncyc={ncyc},",
             f"  ntb=1, cut={_fmt(req.cut)},",
             "  ntpr=100,",
             *_restraint_line(step, mask),
             "/"]
    return "\n".join(lines) + "\n"


def _heating(step: ProtocolStep, req: AmberMdinRequest, mask: Optional[str], continues: bool) -> str:
    target = step_temperature(step, req.temp0)
    dt_fs = min(req.dt_fs, HEAT_DT_FS)
    dt_ps = dt_fs / 1000.0
    nstlim = max(1, round(step.length_ps / dt_ps))
    ntc = _NTC_NTF[req.constraints]
    title = (f"Heating {step.name}: {step.temp_start:g} -> {target:g} K, NVT, {step.length_ps:g} ps"
             + (", restrained" if mask else ""))
    lines = [title, "&cntrl",
             # po MD kroku navazuje na jeho rychlosti, jinak je přiřadí podle temp_start
             "  imin=0, irest=1, ntx=5," if continues else "  imin=0, irest=0, ntx=1,",
             f"  nstlim={nstlim}, dt={_fmt(dt_ps)},",
             f"  ntb=1, ntp=0, cut={_fmt(req.cut)},",
             f"  ntc={ntc}, ntf={ntc}, tol={_fmt(req.tol)},",
             ("  ntt=3, " if continues else f"  ntt=3, tempi={_fmt(step.temp_start)}, ")
             + f"temp0={_fmt(target)}, "
             f"gamma_ln={_fmt(max(req.gamma_ln, 1.0))}, ig=-1,",
             *_restraint_line(step, mask),
             f"  ntpr={max(1, nstlim // 100)}, ntwx={max(1, nstlim // 100)}, ntwr={nstlim},",
             f"  ioutfm={req.ioutfm}, ntxo={req.ntxo},",
             "  nmropt=1,",
             "/"]
    # teplota roste lineárně přes prvních 80 % ohřevu, zbytek drží cílovou
    ramp_end = max(1, round(nstlim * 0.8))
    ramp = ("&wt", f"  type='TEMP0', istep1=0, istep2={ramp_end}, value1={_fmt(step.temp_start)}, "
                   f"value2={_fmt(target)},", "/",
            "&wt", f"  type='TEMP0', istep1={ramp_end + 1}, istep2={nstlim}, value1={_fmt(target)}, "
                   f"value2={_fmt(target)},", "/")
    if req.ghbfix:
        lines += ghbfix_wt_lines(nstlim, ramp)
    else:
        lines += [*ramp, "&wt", "  type='END',", "/"]
    return "\n".join(lines) + "\n"


def _equilibration(step: ProtocolStep, req: AmberMdinRequest, mask: Optional[str], continues: bool) -> str:
    # Termostat Langevin, pokud produkce žádný nemá; v NPT Monte Carlo
    # barostat, pokud produkce žádný nemá - i když produkce běží v NVT.
    thermostat = req.thermostat if req.thermostat not in (None, "none") else "langevin"
    barostat = "none"
    if step.pressure:
        barostat = req.barostat if req.barostat not in (None, "none") else "monte-carlo"
    temp = step_temperature(step, req.temp0)
    equil = req.model_copy(update={
        "temp0": temp, "duration_ns": step.length_ps / 1000.0, "segments": 1,
        # navazuje na rychlosti předchozího MD kroku, po minimalizaci začíná znovu
        "protocol": continues, "run_type": "continue" if continues else "new",
        "ensemble": "NPT" if step.pressure else "NVT", "thermostat": thermostat, "barostat": barostat,
    })
    lines = render_amber_mdin(equil).split("\n")
    lines[0] = (f"Equilibration {step.name}: {'NPT' if step.pressure else 'NVT'}, {temp:g} K, {step.length_ps:g} ps"
                + (", restrained" if mask else ""))
    if not continues:
        lines = [line.replace("irest=0, ntx=1,", f"irest=0, ntx=1, tempi={_fmt(temp)},") for line in lines]
    lines[2:2] = _restraint_line(step, mask)  # hned za "&cntrl"
    return "\n".join(lines)


def _run_script(req: AmberMdinRequest, steps: List[ProtocolStep], masks: List[Optional[str]]) -> str:
    segments = req.segments
    calls = []
    start = '"$CRD"'
    for step, mask in zip(steps, masks):
        ref = ' "$CRD"' if mask else ""
        calls.append(f"run_step {step.name} {step.name}.in {start}{ref}")
        start = f"{step.name}.rst7"
    calls_text = "\n".join(calls)
    chain = " -> ".join(s.name for s in steps)
    return f"""#!/usr/bin/env bash
# FORGE MD run: {chain} -> prod_001..prod_{segments:03d}
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

{calls_text}

prev={steps[-1].name}
for i in $(seq 1 "$SEGMENTS"); do
    name=$(printf "prod_%03d" "$i")
    run_step "$name" prod.in "$prev.rst7"
    prev=$name
done
echo "== all steps finished"
"""


def render_protocol(req: AmberMdinRequest, pdb_text: Optional[str]) -> List[Tuple[str, str]]:
    """(jméno souboru, obsah) všech vstupů běhu včetně run_md.sh."""
    try:
        steps = check_steps(req.steps) if req.steps is not None else default_steps()
    except ValueError as e:
        raise BadRequestError(str(e), code="protocol_steps_invalid")
    prod = req.model_copy(update={"protocol": True})
    solute_mask = _restraint_mask(solute_residue_count(pdb_text) if pdb_text else None)
    files: List[Tuple[str, str]] = []
    masks: List[Optional[str]] = []
    previous_md = False
    for step in steps:
        mask = _step_mask(step, solute_mask)
        if step.kind == "minimization":
            text = _minimization(step, req, mask)
        elif step.kind == "heating":
            text = _heating(step, req, mask, continues=previous_md)
        else:
            text = _equilibration(step, req, mask, continues=previous_md)
        files.append((f"{step.name}.in", _with_extra(text, step)))
        masks.append(mask)
        previous_md = is_md(step)
    files.append((_mdin_filename(prod), render_amber_mdin(prod)))
    files.append((RUN_SCRIPT, _run_script(req, steps, masks)))
    return files


def render_workspace_protocol(req: AmberMdinRequest, workspace_id: str) -> List[Tuple[str, str]]:
    """render_protocol se strukturou workspace (Export i náhled v Simulation panelu)."""
    pdb_path = workspace_manager.get_file_path(workspace_id, STRUCTURE_FILE)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace") if pdb_path.exists() else None
    return render_protocol(req, pdb_text)
