"""
Popis postupu jako souvislý text (druhé PDF v Exportu, vedle tabulkového
reportu v pdf.py): co uživatel načetl, co rozhodl a jak systém připravil,
ve stylu sekce Methods. Ze stejných dat jako report - záznamy kroků,
složení výsledné struktury a nastavení simulace.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Sequence

from app.core.config import settings

from .pdf import _CONCENTRATION_MODE, _Report, _num, system_composition
from .records import CONFORMATIONS, PREPARATION, SOURCE, TOPOLOGY, load_records

PROCEDURE_FILENAME = "forge_procedure.pdf"

_MOLECULE_NOUNS = {
    "P": "the protein", "R": "the RNA", "D": "the DNA", "W": "water",
    "I1": "Na+, K+ and Cl- ions", "I1+": "Li+, Rb+, Cs+, F-, Br- and I- ions",
    "Im": "Mg2+ ions", "Im+": "multivalent metal ions",
}
_BOX_ADJECTIVE = {"cube": "cubic", "octahedron": "truncated octahedral", "truncated octahedron": "truncated octahedral"}
_THERMOSTAT_TEXT = {
    "langevin": "a Langevin thermostat", "berendsen": "a Berendsen thermostat", "andersen": "an Andersen thermostat",
}
_BAROSTAT_TEXT = {"monte-carlo": "a Monte Carlo barostat", "berendsen": "a Berendsen barostat"}
_CONSTRAINTS_TEXT = {
    "h-bonds": "bonds involving hydrogen were constrained with SHAKE",
    "all-bonds": "all bonds were constrained with SHAKE",
    "none": "no bond constraints were used",
}


def _join(items: Sequence[str]) -> str:
    items = [i for i in items if i]
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _ions_text(counts: Optional[Dict[str, Any]]) -> str:
    return _join([f"{v} {str(k).split(':')[-1]}" for k, v in (counts or {}).items()])


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" + ("" if n == 1 else "s")


def _ff_name(ff: Dict[str, Any]) -> str:
    return ff.get("display_name") or ff.get("ff_name") or "?"


class _References:
    """Číslované citace FF v pořadí prvního výskytu v textu."""

    def __init__(self):
        self.entries: List[str] = []
        self._index: Dict[str, int] = {}

    def cite(self, ff: Dict[str, Any]) -> str:
        name = _ff_name(ff)
        dois = [d.strip() for d in str(ff.get("reference_article_doi") or "").split(",") if d.strip()]
        if not dois and not ff.get("creators"):
            return ""
        if name not in self._index:
            authors = ff.get("creators") or []
            text = name
            if authors:
                text += ". " + "; ".join(authors[:6]).rstrip(".") + (" et al" if len(authors) > 6 else "")
            if dois:
                text += ". " + ", ".join(f"doi:{d}" for d in dois)
            self.entries.append(text + ".")
            self._index[name] = len(self.entries)
        return f" [{self._index[name]}]"


# ---------------------------------------------------------------------------
# Odstavce
# ---------------------------------------------------------------------------

def _input_paragraph(records: Dict[str, Dict[str, Any]]) -> str:
    source = records.get(SOURCE) or {}
    if source.get("kind") == "pdb":
        text = f"The starting structure, PDB entry {source.get('pdb_code', '').upper()}, was downloaded from the RCSB Protein Data Bank."
    elif source.get("filename"):
        text = f"The starting structure was loaded from the uploaded file {source['filename']}."
    else:
        text = "The starting structure was loaded into FORGE."

    conf = records.get(CONFORMATIONS)
    if not conf:
        return text + " It contained a single model without alternative locations, so no conformational choices were needed."
    parts = [f"Model {conf.get('model', 1)} was used"]
    parts.append("the biological assembly was built from the symmetry operators in REMARK 350"
                 if conf.get("apply_symmetry") else "the asymmetric unit was kept as deposited")
    text += " " + _join(parts) + "."
    selection = conf.get("selection") or {}
    if selection:
        positions = sorted(set(selection.values()))
        if len(positions) == 1:
            text += (f" Alternative locations of {_plural(len(selection), 'residue')} were resolved by keeping "
                     f"position {positions[0]}.")
        else:
            chosen = _join([f"{res.replace('_', ' ')} ({alt})" for res, alt in sorted(selection.items())])
            text += f" Alternative locations were resolved for {_plural(len(selection), 'residue')}: {chosen}."
    removed = conf.get("remove_chains") or []
    if removed:
        text += f" Chain{'s' if len(removed) > 1 else ''} {_join(removed)} {'were' if len(removed) > 1 else 'was'} removed as identical copies of the molecule."
    return text


def _force_field_paragraph(ffs: Dict[str, Dict[str, Any]], refs: _References, topology_done: bool) -> str:
    if not ffs:
        return ""
    described = [f"{_MOLECULE_NOUNS.get(group, group)} by {_ff_name(ff)}{refs.cite(ff)}" for group, ff in ffs.items()]
    text = f"Force fields were assigned per molecule type: {_join(described)}."
    for ff in ffs.values():
        name = _ff_name(ff)
        if ff.get("has_ghbfix"):
            text += (f" The gHBfix hydrogen-bond correction of {name} was included as restraints read by AMBER from "
                     "ghbfix.f (nmropt=1, DISANG).")
        if ff.get("has_nbfix"):
            text += f" The NBfix pair corrections of {name} were not applied to the topology."
    if not topology_done:
        text += " (The topology had not been generated yet when this description was exported.)"
    return text


def _preparation_paragraphs(prep: Optional[Dict[str, Any]], water_name: Optional[str]) -> List[str]:
    if not prep:
        return ["The structure was not prepared further."]
    s = prep.get("settings") or {}
    summary = prep.get("summary") or {}
    out: List[str] = []

    crystal = {
        "remove_all": "All crystal waters and ions were removed",
        "keep_water": "Crystal waters were kept",
        "keep_all": "Crystal waters and ions were kept",
    }.get(s.get("crystal_water_mode"), "")
    text = (f"Missing atoms, including hydrogens, were built by the FORGE builder and protonation states were "
            f"assigned at pH {_num(s.get('ph'))}.")
    assignments = summary.get("protonation_assignments") or []
    if assignments:
        listed = _join([f"{a.get('residue')} as {a.get('assigned')}" for a in assignments[:10]])
        more = f" and {len(assignments) - 10} more" if len(assignments) > 10 else ""
        text += f" Non-default states were assigned to {listed}{more}."
    if s.get("protonation_overrides"):
        text += f" The protonation state of {_plural(len(s['protonation_overrides']), 'residue')} was set by hand."
    disulfides = summary.get("disulfide_bonds") or []
    if disulfides:
        pairs = _join([f"{b.get('atom1')}-{b.get('atom2')}" for b in disulfides])
        verb = "were" if len(disulfides) > 1 else "was"
        text += f" {_plural(len(disulfides), 'disulfide bond')} ({pairs}) {verb} formed."
    if s.get("auto_amide_flips"):
        text += " Recommended ASN/GLN amide flips were applied automatically."
    else:
        flips = [d for d in (s.get("structure_decisions") or {}).get("amide_flips") or [] if d.get("apply")]
        if flips:
            text += f" {_plural(len(flips), 'ASN/GLN amide')} {'were' if len(flips) > 1 else 'was'} flipped after review."
    if prep.get("interactive_sidechains"):
        text += " Missing side-chain conformations were set interactively before the structure was completed."
    if crystal:
        text += f" {crystal}."
    cleanup = summary.get("crystal_ion_cleanup") or {}
    if cleanup.get("input_ions"):
        removed = (cleanup.get("removed_monovalent") or 0) + (cleanup.get("removed_nonstructural_multivalent") or 0)
        text += (f" Of {cleanup['input_ions']} crystal ions, {removed} non-structural "
                 f"{'ion was' if removed == 1 else 'ions were'} removed")
        if cleanup.get("replaced_by_magnesium"):
            text += f" and {cleanup['replaced_by_magnesium']} replaced by Mg2+"
        text += "."
    out.append(text)

    if s.get("add_solvent"):
        solv = summary.get("solvation") or {}
        ions = summary.get("ion_addition") or {}
        shape = _BOX_ADJECTIVE.get(s.get("box_shape"), s.get("box_shape") or "")
        text = (f"The solute was placed in a {shape} periodic box with {_num(s.get('box_padding_nm'))} nm padding "
                f"and solvated with {solv.get('generated_waters', 'the required number of')} "
                f"{water_name or ''} water molecules.").replace("  ", " ")
        salts = [f"{s.get('positive_ion')}/{s.get('negative_ion')} at {_num(s.get('ionic_strength'))} M"]
        salts += [f"{x.get('positive_ion')}/{x.get('negative_ion')} at {_num(x.get('ionic_strength'))} M"
                  for x in s.get("additional_salts") or []]
        if ions.get("neutralization_ions"):
            text += f" The system was neutralized with {_ions_text(ions['neutralization_ions'])}"
            text += f" and salt was added ({_join(salts)}"
        else:
            text += f" Salt was added ({_join(salts)}"
        basis = _CONCENTRATION_MODE.get(s.get("concentration_mode"), "")
        text += f", concentration from the {basis.lower()})" if basis else ")"
        if ions.get("added_ions"):
            text += f", giving {_ions_text(ions['added_ions'])} in total"
        text += "."
        if ions.get("final_system_charge") is not None:
            text += f" The final net charge was {_num(ions['final_system_charge'])}."
        out.append(text)
    else:
        out.append("The system was not solvated (vacuum).")
    return out


def _system_paragraph(composition: Dict[str, Any]) -> str:
    if not composition:
        return ""
    groups = composition.get("residues_by_group") or {}
    parts = [f"{groups[g]} {_MOLECULE_NOUNS[g].replace('the ', '')} residues" for g in ("P", "R", "D") if groups.get(g)]
    waters = sum(v for k, v in groups.items() if str(k).startswith("W"))
    if waters:
        parts.append(f"{waters} water molecules")
    ions = composition.get("ions") or {}
    if ions:
        parts += [f"{v} {k}" for k, v in sorted(ions.items())]
    text = f"The final system contained {composition.get('atoms', 0)} atoms"
    if parts:
        text += f" ({_join(parts)})"
    box = composition.get("box")
    if box:
        text += (f" in a periodic box with a = {box[0]:.2f}, b = {box[1]:.2f}, c = {box[2]:.2f} A and angles "
                 f"{box[3]:.1f}, {box[4]:.1f}, {box[5]:.1f} deg")
    return text + "."


def _topology_paragraph(topo: Optional[Dict[str, Any]]) -> str:
    if not topo:
        return ""
    text = "An AMBER topology (structure.prmtop) and starting coordinates (structure.crd) were generated"
    if topo.get("hmr"):
        text += " with hydrogen mass repartitioning"
    return text + "."


def _simulation_paragraph(sim: Optional[Dict[str, Any]], ensemble: Optional[str], nstlim: Optional[int]) -> str:
    if not sim:
        return ""
    text = (f"An AMBER input file was prepared for {'a new' if sim.get('run_type') == 'new' else 'a continued'} "
            f"{_num(sim.get('duration_ns'))} ns {ensemble or ''} production segment"
            + (f" ({nstlim} steps)" if nstlim else "")
            + f" with a {_num(sim.get('dt_fs'))} fs time step").replace("  ", " ")
    if sim.get("hmr"):
        text += " enabled by hydrogen mass repartitioning"
    text += "."
    thermostat = sim.get("thermostat") or "none"
    barostat = sim.get("barostat") or "none"
    coupling = []
    if thermostat != "none":
        t = f"{_THERMOSTAT_TEXT.get(thermostat, thermostat)} at {_num(sim.get('temp0'))} K"
        if thermostat == "langevin":
            t += f" (collision frequency {_num(sim.get('gamma_ln'))} ps-1)"
        coupling.append(t)
    if barostat != "none":
        coupling.append(f"{_BAROSTAT_TEXT.get(barostat, barostat)} at {_num(sim.get('pres0'))} bar")
    if coupling:
        text += f" Temperature and pressure were controlled by {_join(coupling)}." if len(coupling) > 1 \
            else f" Temperature was controlled by {coupling[0]}."
    constraints = _CONSTRAINTS_TEXT.get(sim.get("constraints"), "")
    if constraints:
        text += f" During the run, {constraints}, and a {_num(sim.get('cut'))} A cutoff was used for nonbonded interactions."
    text += (f" Energies are written every {_num(sim.get('log_interval_ps'))} ps, coordinates every "
             f"{_num(sim.get('traj_interval_ps'))} ps and restart files every {_num(sim.get('restart_interval_ps'))} ps.")
    if sim.get("ghbfix"):
        text += " The gHBfix restraints are read from ghbfix.f (DISANG) and listed in ghbfix.out (LISTOUT)."
    return text


def render_procedure_pdf(
    workspace_id: str,
    simulation: Optional[Dict[str, Any]] = None,
    ensemble: Optional[str] = None,
    nstlim: Optional[int] = None,
) -> bytes:
    records = load_records(workspace_id)
    composition = system_composition(workspace_id)
    topo = records.get(TOPOLOGY)
    prep = records.get(PREPARATION)
    ffs = (topo or {}).get("force_fields") or (prep or {}).get("force_fields") or {}
    water = ffs.get("W")
    refs = _References()

    pdf = _Report(running_title="FORGE - procedure description")
    pdf.add_page()
    pdf.h1("Procedure description")
    pdf.para(f"Generated by FORGE {settings.VERSION} on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}.",
             size=8.5, style="I")

    sections = [
        ("Input structure", [_input_paragraph(records)]),
        ("Force fields", [_force_field_paragraph(ffs, refs, topo is not None)]),
        ("System preparation", _preparation_paragraphs(prep, _ff_name(water) if water else None)
         + [_system_paragraph(composition), _topology_paragraph(topo)]),
        ("Simulation input", [_simulation_paragraph(simulation, ensemble, nstlim)]),
    ]
    for title, paragraphs in sections:
        paragraphs = [p for p in paragraphs if p]
        if not paragraphs:
            continue
        pdf.h2(title)
        for paragraph in paragraphs:
            pdf.para(paragraph, size=10)
            pdf.ln(1.5)

    if refs.entries:
        pdf.h2("References")
        for i, entry in enumerate(refs.entries, 1):
            pdf.para(f"[{i}] {entry}", size=9)
    return bytes(pdf.output())
