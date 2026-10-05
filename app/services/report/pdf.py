"""
Souhrnný report přípravy systému jako PDF (Export, volba "Report").

Skládá se ze záznamů kroků (records.py), složení výsledné struktury
(structure.pdb + structure.forge_meta.json) a nastavení simulace z Exportu.
Písmo je vestavěná Helvetica (latin-1) - text se na ni převádí v _t().
"""

from __future__ import annotations

import json
import time
from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from fpdf import FPDF
from fpdf.fonts import FontFace

from app.core.config import settings
from app.workspaces.manager import workspace_manager

from .records import CONFORMATIONS, PREPARATION, SOURCE, TOPOLOGY, load_records

REPORT_FILENAME = "forge_report.pdf"

_GROUP_LABELS = {
    "P": "Protein", "R": "RNA", "D": "DNA", "W": "Water",
    "W3": "Water (3-point)", "W4": "Water (4-point)", "W5": "Water (5-point)",
    "I1": "Ions (Na+, K+, Cl-)", "I1+": "Ions (Li+, Rb+, Cs+, F-, Br-, I-)",
    "Im": "Ions (Mg2+)", "Im+": "Ions (multivalent metals)",
}
_MODE_LABELS = {"guided": "Guided", "standard": "Standard", "expert": "Expert"}
_CRYSTAL_WATER = {
    "remove_all": "Remove all (keep only solute)",
    "keep_water": "Keep crystal waters",
    "keep_all": "Keep crystal waters and ions",
}
_CONCENTRATION_MODE = {"water_ratio": "Ion-to-water ratio", "box_volume": "Periodic-box volume"}
_BOX_SHAPE = {"cube": "Cube", "octahedron": "Truncated octahedron", "truncated octahedron": "Truncated octahedron"}
_THERMOSTAT = {"langevin": "Langevin", "berendsen": "Berendsen", "andersen": "Andersen"}
_BAROSTAT = {"monte-carlo": "Monte Carlo", "berendsen": "Berendsen"}
_LIST_LIMIT = 40

_REPLACEMENTS = {
    "→": "->", "–": "-", "—": "-", "−": "-", "²": "2", "³": "3",
    "⁺": "+", "⁻": "-", "‘": "'", "’": "'", "“": '"', "”": '"',
    "…": "...", "·": "-",
}


def _t(value: Any) -> str:
    text = str(value)
    for src, dst in _REPLACEMENTS.items():
        text = text.replace(src, dst)
    return text.encode("latin-1", errors="replace").decode("latin-1")


def _num(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


# ---------------------------------------------------------------------------
# Složení výsledného systému
# ---------------------------------------------------------------------------

def system_composition(workspace_id: str) -> Dict[str, Any]:
    ws_dir = workspace_manager.get_workspace_dir(workspace_id)
    pdb_path = ws_dir / "structure.pdb"
    meta_path = ws_dir / "structure.forge_meta.json"
    out: Dict[str, Any] = {}
    if pdb_path.exists():
        atoms = 0
        for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                atoms += 1
            elif line.startswith("CRYST1"):
                try:
                    out["box"] = [float(line[6:15]), float(line[15:24]), float(line[24:33]),
                                  float(line[33:40]), float(line[40:47]), float(line[47:54])]
                except ValueError:
                    pass
        out["atoms"] = atoms
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            meta = {}
        groups = Counter(entry.get("group") for entry in meta.values())
        out["residues_by_group"] = dict(groups)
        out["ions"] = dict(Counter(
            entry.get("ff_resname") for entry in meta.values() if str(entry.get("group", "")).startswith("I")
        ))
    return out


# ---------------------------------------------------------------------------
# Sazba
# ---------------------------------------------------------------------------

class _Report(FPDF):
    def __init__(self, running_title: str = "FORGE - system preparation report"):
        super().__init__(format="A4")
        self.running_title = running_title
        self.set_auto_page_break(auto=True, margin=15)
        self.set_margins(15, 12, 15)
        self.alias_nb_pages()

    def header(self):
        self.set_font("Helvetica", "B", 9)
        self.set_text_color(110, 110, 110)
        self.cell(0, 6, _t(self.running_title), align="L")
        self.ln(8)
        self.set_text_color(0, 0, 0)

    def footer(self):
        self.set_y(-12)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(110, 110, 110)
        self.cell(0, 6, f"Page {self.page_no()}/{{nb}}", align="C")
        self.set_text_color(0, 0, 0)

    def h1(self, text: str):
        self.set_font("Helvetica", "B", 17)
        self.multi_cell(0, 9, _t(text), new_x="LMARGIN", new_y="NEXT")
        self.ln(1)

    def h2(self, text: str):
        self.ln(3)
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(232, 240, 238)
        self.cell(0, 7, _t(text), fill=True, new_x="LMARGIN", new_y="NEXT")
        self.ln(1.5)

    def para(self, text: str, size: float = 9.5, style: str = ""):
        self.set_font("Helvetica", style, size)
        self.multi_cell(0, 5, _t(text), new_x="LMARGIN", new_y="NEXT")

    def kv(self, rows: Iterable[Tuple[str, Any]]):
        rows = [(k, v) for k, v in rows if v not in (None, "", [], {})]
        if not rows:
            return
        self.set_font("Helvetica", "", 9.5)
        with self.table(col_widths=(55, 125), first_row_as_headings=False, borders_layout="HORIZONTAL_LINES",
                        line_height=5.5, padding=1) as table:
            for key, value in rows:
                row = table.row()
                row.cell(_t(key))
                row.cell(_t(value))
        self.ln(1)

    def grid(self, headings: Sequence[str], rows: List[Sequence[Any]], widths: Sequence[float]):
        if not rows:
            return
        self.set_font("Helvetica", "", 9)
        with self.table(col_widths=tuple(widths), borders_layout="HORIZONTAL_LINES", line_height=5, padding=1,
                        headings_style=FontFace(emphasis="BOLD", fill_color=(245, 245, 245))) as table:
            head = table.row()
            for h in headings:
                head.cell(_t(h))
            for values in rows:
                row = table.row()
                for v in values:
                    row.cell(_t(v))
        self.ln(1)

    def bullets(self, items: Sequence[str], limit: int = _LIST_LIMIT):
        self.set_font("Helvetica", "", 9.5)
        for item in items[:limit]:
            self.multi_cell(0, 5, _t(f"- {item}"), new_x="LMARGIN", new_y="NEXT")
        if len(items) > limit:
            self.multi_cell(0, 5, _t(f"... and {len(items) - limit} more"), new_x="LMARGIN", new_y="NEXT")


# ---------------------------------------------------------------------------
# Sekce
# ---------------------------------------------------------------------------

def _section_input(pdf: _Report, records: Dict[str, Dict[str, Any]]):
    source = records.get(SOURCE)
    conf = records.get(CONFORMATIONS)
    pdf.h2("1. Input structure")
    if source:
        origin = (f"PDB entry {source.get('pdb_code', '').upper()} (RCSB)" if source.get("kind") == "pdb"
                  else f"Uploaded file {source.get('filename')}")
        pdf.kv([("Source", origin)])
    else:
        pdf.para("Source not recorded.")
    if conf:
        selection = conf.get("selection") or {}
        pdf.kv([
            ("Model", conf.get("model")),
            ("Biological assembly", "Built from symmetry operators (REMARK 350)" if conf.get("apply_symmetry")
             else "Asymmetric unit"),
            ("Alternative locations", f"{len(selection)} residue(s) resolved" if selection else "None"),
            ("Removed chains", ", ".join(conf.get("remove_chains") or []) or "None"),
        ])
        if selection:
            pdf.para("Chosen alternative locations:", style="B")
            pdf.bullets([f"{res}: position {alt}" for res, alt in sorted(selection.items())])
    else:
        pdf.para("No model, assembly or alternative-location decisions were needed.")


def _ff_rows(ffs: Dict[str, Dict[str, Any]]) -> List[Sequence[Any]]:
    rows = []
    for group, ff in ffs.items():
        name = ff.get("display_name") or ff.get("ff_name") or "?"
        notes = []
        if ff.get("has_ghbfix"):
            notes.append("gHBfix restraints (ghbfix.f)")
        if ff.get("has_nbfix"):
            notes.append("NBfix pairs (not applied)")
        rows.append((_GROUP_LABELS.get(group, group), name, "; ".join(notes) or "-"))
    return rows


def _section_force_fields(pdf: _Report, records: Dict[str, Dict[str, Any]]):
    pdf.h2("2. Force fields")
    topo = records.get(TOPOLOGY)
    prep = records.get(PREPARATION)
    ffs = (topo or {}).get("force_fields") or (prep or {}).get("force_fields") or {}
    if not ffs:
        pdf.para("No force fields recorded.")
        return
    pdf.grid(("Molecule type", "Force field", "Notes"), _ff_rows(ffs), (45, 65, 70))
    if not topo:
        pdf.para("Topology was not generated yet - the list shows the force fields offered for preparation.",
                 size=8.5, style="I")


def _section_preparation(pdf: _Report, records: Dict[str, Dict[str, Any]]):
    pdf.h2("3. Structure preparation")
    prep = records.get(PREPARATION)
    if not prep:
        pdf.para("The structure was not prepared.")
        return
    s = prep.get("settings") or {}
    salts = [f"{s.get('positive_ion')}/{s.get('negative_ion')} {_num(s.get('ionic_strength'))} M"]
    salts += [f"{x.get('positive_ion')}/{x.get('negative_ion')} {_num(x.get('ionic_strength'))} M"
              for x in s.get("additional_salts") or []]
    rows: List[Tuple[str, Any]] = [
        ("pH", _num(s.get("ph"))),
        ("Crystal waters and ions", _CRYSTAL_WATER.get(s.get("crystal_water_mode"), s.get("crystal_water_mode"))),
        ("Clean up crystal ions", "Yes" if s.get("clean_crystal_ions") else "No"),
        ("Replace multivalent ions with Mg2+", "Yes" if s.get("replace_structural_multivalent_with_mg") else "No"),
        ("Solvation", "Water box and ions" if s.get("add_solvent") else "None (vacuum)"),
    ]
    if s.get("add_solvent"):
        rows += [
            ("Box shape", _BOX_SHAPE.get(s.get("box_shape"), s.get("box_shape"))),
            ("Box padding", f"{_num(s.get('box_padding_nm'))} nm"),
            ("Salts", ", ".join(salts)),
            ("Concentration basis", _CONCENTRATION_MODE.get(s.get("concentration_mode"), s.get("concentration_mode"))),
        ]
    decisions = s.get("structure_decisions") or {}
    flips = [d for d in decisions.get("amide_flips") or [] if d.get("apply")]
    rebuilds = [d for d in decisions.get("zero_occupancy") or [] if d.get("apply")]
    rows += [
        ("ASN/GLN amide flips", "Recommended flips applied automatically" if s.get("auto_amide_flips")
         else f"{len(flips)} chosen in Structure Check"),
        ("Rebuilt zero-occupancy side chains", len(rebuilds) or None),
        ("Protonation states set by hand", len(s.get("protonation_overrides") or []) or None),
        ("Side chains built interactively", "Yes" if prep.get("interactive_sidechains") else None),
    ]
    pdf.kv(rows)

    summary = prep.get("summary") or {}
    assignments = summary.get("protonation_assignments") or []
    if assignments:
        pdf.para("Non-default protonation states:", style="B")
        pdf.bullets([f"{a.get('residue')}: {a.get('default')} -> {a.get('assigned')}" for a in assignments])
    disulfides = summary.get("disulfide_bonds") or []
    if disulfides:
        pdf.para("Disulfide bonds:", style="B")
        pdf.bullets([f"{b.get('atom1')} - {b.get('atom2')} ({_num(b.get('distance_angstrom'))} A)" for b in disulfides])
    cleanup = summary.get("crystal_ion_cleanup")
    if cleanup and cleanup.get("input_ions"):
        pdf.para("Crystal ions:", style="B")
        pdf.kv([(k.replace("_", " ").capitalize(), _num(v) if not isinstance(v, (list, dict)) else
                 (", ".join(map(str, v)) if v else None)) for k, v in cleanup.items()])
    solv = summary.get("solvation")
    ions = summary.get("ion_addition")
    if solv or ions:
        pdf.para("Solvation and ions:", style="B")
        pdf.kv([
            ("Crystal waters kept", solv and f"{solv.get('retained_crystal_waters')} of {solv.get('input_crystal_waters')}"),
            ("Waters added", solv and solv.get("generated_waters")),
            # Ionty se pak umístí místo části vod - konečný počet je v sekci 4.
            ("Waters before ion placement", solv and solv.get("total_waters")),
            ("Neutralizing ions", ions and _ion_counts(ions.get("neutralization_ions"))),
            ("Ions added in total", ions and _ion_counts(ions.get("added_ions"))),
            ("Final system charge", ions and _num(ions.get("final_system_charge"))),
        ])
    warnings = summary.get("protonation_warnings") or []
    if warnings:
        pdf.para("Warnings:", style="B")
        pdf.bullets(warnings)


def _ion_counts(counts: Optional[Dict[str, Any]]) -> str:
    """{"I1:Na+": 8} -> "Na+ x8" (skupina z klíče builderu pryč)."""
    return ", ".join(f"{str(k).split(':')[-1]} x{v}" for k, v in (counts or {}).items())


def _section_system(pdf: _Report, composition: Dict[str, Any]):
    pdf.h2("4. Final system")
    if not composition:
        pdf.para("No structure available.")
        return
    groups = composition.get("residues_by_group") or {}
    rows: List[Tuple[str, Any]] = [("Atoms", composition.get("atoms"))]
    for group in ("P", "R", "D"):
        if groups.get(group):
            rows.append((f"{_GROUP_LABELS[group]} residues", groups[group]))
    waters = sum(v for k, v in groups.items() if str(k).startswith("W"))
    rows.append(("Water molecules", waters or None))
    ions = composition.get("ions") or {}
    rows.append(("Ions", ", ".join(f"{k} x{v}" for k, v in sorted(ions.items())) or None))
    box = composition.get("box")
    if box:
        rows.append(("Periodic box", f"a={box[0]:.2f} b={box[1]:.2f} c={box[2]:.2f} A, "
                                     f"angles {box[3]:.1f} / {box[4]:.1f} / {box[5]:.1f} deg"))
    pdf.kv(rows)


def _section_simulation(pdf: _Report, simulation: Optional[Dict[str, Any]], ensemble: Optional[str],
                        nstlim: Optional[int]):
    pdf.h2("5. Simulation input (AMBER)")
    if not simulation:
        pdf.para("Simulation settings were not exported.")
        return
    s = simulation
    thermostat = s.get("thermostat") or "none"
    barostat = s.get("barostat") or "none"
    rows: List[Tuple[str, Any]] = [
        ("Segment", "New run" if s.get("run_type") == "new" else "Continuation (restart)"),
        ("Ensemble", ensemble),
        ("Length", f"{_num(s.get('duration_ns'))} ns ({nstlim} steps)" if nstlim else f"{_num(s.get('duration_ns'))} ns"),
        ("Time step", f"{_num(s.get('dt_fs'))} fs"),
        ("Hydrogen mass repartitioning", "Yes" if s.get("hmr") else "No"),
        ("Thermostat", f"{_THERMOSTAT.get(thermostat, thermostat)}, {_num(s.get('temp0'))} K"
         if thermostat != "none" else "None"),
        ("Barostat", f"{_BAROSTAT.get(barostat, barostat)}, {_num(s.get('pres0'))} bar"
         if barostat != "none" else "None"),
        ("Constraints", s.get("constraints")),
        ("Nonbonded cutoff", f"{_num(s.get('cut'))} A"),
        ("Output", f"log every {_num(s.get('log_interval_ps'))} ps, trajectory every "
                   f"{_num(s.get('traj_interval_ps'))} ps, restart every {_num(s.get('restart_interval_ps'))} ps"),
        ("gHBfix restraints", "Yes (nmropt=1, DISANG=ghbfix.f)" if s.get("ghbfix") else "No"),
    ]
    pdf.kv(rows)


def _section_files(pdf: _Report, files: Sequence[str]):
    pdf.h2("6. Exported files")
    pdf.bullets(list(files) or ["(none)"])


def _section_references(pdf: _Report, records: Dict[str, Dict[str, Any]]):
    ffs = ((records.get(TOPOLOGY) or {}).get("force_fields")
           or (records.get(PREPARATION) or {}).get("force_fields") or {})
    entries = []
    seen = set()
    for ff in ffs.values():
        name = ff.get("display_name") or ff.get("ff_name")
        if not name or name in seen:
            continue
        seen.add(name)
        dois = [d.strip() for d in str(ff.get("reference_article_doi") or "").split(",") if d.strip()]
        parts = [name]
        if ff.get("creators"):
            authors = ff["creators"]
            parts.append("; ".join(authors[:6]) + (" et al." if len(authors) > 6 else ""))
        if dois:
            parts.append("Reference: " + ", ".join(f"doi:{d}" for d in dois))
        if ff.get("doi"):
            parts.append(f"Dataset: doi:{ff['doi']}")
        entries.append(". ".join(p.rstrip(".") for p in parts) + ".")
    if not entries:
        return
    pdf.h2("7. Force field references")
    pdf.bullets(entries, limit=100)


def render_report_pdf(
    workspace_id: str,
    simulation: Optional[Dict[str, Any]] = None,
    ensemble: Optional[str] = None,
    nstlim: Optional[int] = None,
    detail_level: Optional[str] = None,
    files: Sequence[str] = (),
) -> bytes:
    records = load_records(workspace_id)
    composition = system_composition(workspace_id)

    pdf = _Report()
    pdf.add_page()
    pdf.h1("System preparation report")
    pdf.kv([
        ("Generated", time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())),
        ("Application", f"FORGE {settings.VERSION}"),
        ("Detail level", _MODE_LABELS.get(detail_level or "", detail_level)),
    ])

    _section_input(pdf, records)
    _section_force_fields(pdf, records)
    _section_preparation(pdf, records)
    _section_system(pdf, composition)
    _section_simulation(pdf, simulation, ensemble, nstlim)
    _section_files(pdf, files)
    _section_references(pdf, records)
    return bytes(pdf.output())
