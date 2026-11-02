"""
Přípravné kroky celého MD běhu (minimalizace, ohřev, ekvilibrace) jako seznam,
který jde v Simulation panelu upravit: přidat, ubrat, přeházet a nastavit.
Bez seznamu od klienta platí DEFAULT_STEPS - stejné kroky jako dřív.

Vstupy z kroků skládá app/services/md_protocol.py, tady je jen model,
výchozí seznam, kontrola seznamu a popis kroků pro report.

Až dorazí referenční vstupy skupiny (01_minimize_hydrogens.in ...), stačí
podle nich upravit DEFAULT_STEPS; parametry, které model nemá, jdou zatím
do `extra` (další řádky &cntrl).
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

StepKind = Literal["minimization", "heating", "equilibration"]
Restraint = Literal["none", "solute", "heavy", "custom"]

MAX_STEPS = 50
MAX_TEMPERATURE = 1000.0   # K - víc je u biomolekul jistě překlep
MIN_MD_LENGTH_PS = 1.0     # kratší MD krok by měl pár (nebo nula) kroků integrace
# Teplotní průběh ohřevu (temp_points): lomená čára, každý úsek = jeden &wt TEMP0.
MAX_TEMP_POINTS = 20
MIN_TEMP_SEGMENT_PS = 0.1  # kratší úsek by po zaokrouhlení na kroky integrace zmizel
_MD_KINDS = ("heating", "equilibration")
_EXTRA_ITEM = re.compile(r"^([a-z_][a-z0-9_]*)\s*=\s*([^=,/&]+)$", re.IGNORECASE)


class TempPoint(BaseModel):
    """Bod teplotního průběhu ohřevu: čas od začátku kroku (ps) a teplota (K)."""
    model_config = ConfigDict(allow_inf_nan=False)

    time_ps: float = Field(ge=0, le=1_000_000)
    temp: float = Field(ge=0, le=MAX_TEMPERATURE)


class ProtocolStep(BaseModel):
    # NaN a nekonečno by prošly mezemi (ge/le) až do mdin
    model_config = ConfigDict(allow_inf_nan=False)

    # jméno souborů kroku (NAME.in, NAME.rst7 ...) - jen bezpečné znaky
    name: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]{0,39}$")
    kind: StepKind
    # Restrainty: solute = těžké atomy solutu, heavy = všechny těžké atomy
    # (i voda a ionty), custom = vlastní AMBER maska.
    restraint: Restraint = "none"
    restraint_mask: str = Field("", max_length=200)
    restraint_wt: float = Field(10.0, gt=0, le=100_000)  # kcal/mol/A^2
    # minimalizace
    maxcyc: int = Field(5000, ge=1, le=10_000_000)
    ncyc: Optional[int] = Field(None, ge=0)          # None = polovina maxcyc
    # ohřev a ekvilibrace
    length_ps: float = Field(100.0, ge=MIN_MD_LENGTH_PS, le=1_000_000)
    temp_start: float = Field(0.0, ge=0, le=MAX_TEMPERATURE)   # ohřev: počáteční teplota
    # ohřev: konečná teplota, ekvilibrace: držená teplota; None = temp0 produkce
    temp_target: Optional[float] = Field(None, ge=0, le=MAX_TEMPERATURE)
    # ohřev: teplotní průběh jako lomená čára (první bod v 0 ps, poslední na
    # konci kroku). None = lineární rampa temp_start -> temp_target přes 80 %
    # kroku a pak držení. Když je zadaný, temp_start/temp_target se z něj odvodí.
    temp_points: Optional[List[TempPoint]] = Field(None, max_length=MAX_TEMP_POINTS)
    pressure: bool = True                            # ekvilibrace: NPT (jinak NVT)
    # další položky &cntrl ("drms=0.01, dx0=0.01") - parametry, které model nemá
    extra: str = Field("", max_length=500)

    @model_validator(mode="after")
    def _check(self):
        if self.ncyc is not None and self.ncyc > self.maxcyc:
            raise ValueError(f"step {self.name}: ncyc must not exceed maxcyc")
        if self.restraint == "custom" and not self.restraint_mask.strip():
            raise ValueError(f"step {self.name}: a custom restraint needs a mask")
        if "'" in self.restraint_mask or "\n" in self.restraint_mask:
            raise ValueError(f"step {self.name}: the restraint mask must not contain quotes or line breaks")
        extra_items(self.extra, self.name)
        if self.temp_points is not None:
            self._check_temp_points()
        return self

    def _check_temp_points(self):
        points = self.temp_points
        if self.kind != "heating":
            raise ValueError(f"step {self.name}: a temperature schedule is only for heating steps")
        if len(points) < 2:
            raise ValueError(f"step {self.name}: the temperature schedule needs at least 2 points")
        if points[0].time_ps != 0:
            raise ValueError(f"step {self.name}: the temperature schedule must start at 0 ps")
        for a, b in zip(points, points[1:]):
            if b.time_ps - a.time_ps < MIN_TEMP_SEGMENT_PS - 1e-9:
                raise ValueError(f"step {self.name}: schedule points must be at least "
                                 f"{MIN_TEMP_SEGMENT_PS:g} ps apart and in time order")
        if abs(points[-1].time_ps - self.length_ps) > 1e-6:
            raise ValueError(f"step {self.name}: the temperature schedule must end at the step length "
                             f"({self.length_ps:g} ps)")
        # začátek a konec ohřevu podle průběhu (tempi, navazování kroků, report)
        self.temp_start = points[0].temp
        self.temp_target = points[-1].temp


def extra_items(extra: str, step_name: str = "") -> List[tuple]:
    """'drms=0.01, dx0=0.01' -> [('drms', '0.01'), ('dx0', '0.01')]; ValueError u nesmyslu."""
    items = []
    for raw in re.split(r"[,\n]", extra):
        raw = raw.strip()
        if not raw:
            continue
        match = _EXTRA_ITEM.match(raw)
        if not match:
            raise ValueError(f"step {step_name}: extra settings must be 'key=value' items, got '{raw}'")
        items.append((match.group(1).lower(), match.group(2).strip()))
    return items


DEFAULT_STEPS: List[ProtocolStep] = [
    ProtocolStep(name="mini0", kind="minimization", restraint="solute", restraint_wt=10.0, maxcyc=5000),
    ProtocolStep(name="mini1", kind="minimization", maxcyc=10000),
    ProtocolStep(name="heat", kind="heating", restraint="solute", restraint_wt=5.0, length_ps=100.0),
    ProtocolStep(name="equil", kind="equilibration", length_ps=1000.0, pressure=True),
]


def default_steps() -> List[ProtocolStep]:
    return [step.model_copy() for step in DEFAULT_STEPS]


def check_steps(steps: List[ProtocolStep]) -> List[ProtocolStep]:
    """Seznam jako celek: počet, jedinečná jména, poslední krok s rychlostmi pro produkci."""
    if not steps:
        raise ValueError("the run needs at least one preparation step")
    if len(steps) > MAX_STEPS:
        raise ValueError(f"at most {MAX_STEPS} preparation steps are allowed")
    seen = set()
    for step in steps:
        key = step.name.lower()
        if key in seen:
            raise ValueError(f"step name '{step.name}' is used twice")
        if key == "prod" or key.startswith("prod_"):
            raise ValueError(f"step name '{step.name}' is reserved for production")
        seen.add(key)
    # Produkce navazuje na restart (irest=1) - potřebuje rychlosti z MD kroku.
    if steps[-1].kind not in _MD_KINDS:
        raise ValueError("the last preparation step must be heating or equilibration - production continues from it")
    return steps


def is_md(step: ProtocolStep) -> bool:
    return step.kind in _MD_KINDS


# ---------------------------------------------------------------------------
# Popis pro report
# ---------------------------------------------------------------------------

def _g(value: Any) -> str:
    return f"{value:g}" if isinstance(value, float) else str(value)


def _length(ps: float) -> str:
    return f"{_g(ps / 1000.0)} ns" if ps >= 1000 else f"{_g(ps)} ps"


def _restraint_text(step: ProtocolStep) -> str:
    target = {
        "solute": "solute heavy atoms",
        "heavy": "all heavy atoms",
        "custom": f"atoms '{step.restraint_mask}'",
    }.get(step.restraint)
    if not target:
        return "unrestrained"
    return f"{target} restrained with {_g(step.restraint_wt)} kcal/mol/A2"


def step_temperature(step: ProtocolStep, temp0: float) -> float:
    """Teplota na konci ohřevu / během ekvilibrace."""
    return step.temp_target if step.temp_target is not None else temp0


def describe_step(step: ProtocolStep, temp0: Optional[float]) -> str:
    if step.kind == "minimization":
        return f"energy minimization ({step.maxcyc} cycles, {_restraint_text(step)})"
    target = step.temp_target if step.temp_target is not None else temp0
    if step.kind == "heating":
        to = f" to {_g(target)} K" if target is not None else ""
        schedule = ""
        if step.temp_points and len(step.temp_points) > 2:
            listed = ", ".join(f"{_g(float(p.temp))} K at {_g(float(p.time_ps))} ps" for p in step.temp_points)
            schedule = f" along a piecewise-linear temperature schedule ({listed})"
        return (f"heating from {_g(step.temp_start)} K{to} over {_length(step.length_ps)}{schedule} "
                f"({_restraint_text(step)})")
    ensemble = "NPT" if step.pressure else "NVT"
    at = f" at {_g(target)} K" if target is not None else ""
    return f"{ensemble} equilibration{at} for {_length(step.length_ps)} ({_restraint_text(step)})"


def steps_from_settings(sim: Dict[str, Any]) -> List[ProtocolStep]:
    """Kroky z uloženého nastavení simulace (model_dump requestu); chybí-li, výchozí."""
    raw = sim.get("steps")
    if not raw:
        return default_steps()
    steps = []
    for item in raw:
        try:
            steps.append(ProtocolStep(**item))
        except ValueError:
            # uložené dřív, než platily dnešní meze - report je jen popíše
            steps.append(ProtocolStep.model_construct(**item))
    return steps
