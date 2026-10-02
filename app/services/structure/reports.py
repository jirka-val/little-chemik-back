"""Human-readable reports of a preparation run: titratable residues, the
histidine review (step 3.5) and the preparation summary shown after Prepare."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, Any, Dict, List, Optional






if TYPE_CHECKING:
    from app.services.structure.forge_service import ForgePreparationResult


def protonation_override_map(overrides: Optional[List[Dict[str, Any]]]) -> Dict[tuple, str]:
    """[{chain, resseq, icode, state}] z požadavku -> klíč rezidua builderu."""
    return {
        (o["chain"], int(o["resseq"]), o.get("icode") or ""): o["state"]
        for o in overrides or []
    }


def _format_residue_key(residue_key: tuple) -> str:
    chain, resseq, icode = residue_key
    return f"{chain}{resseq}{icode}".rstrip()


def _format_atom_key(atom_key: tuple) -> str:
    chain, resseq, icode, atom_name = atom_key
    return f"{_format_residue_key((chain, resseq, icode))}:{atom_name}"


# Builder píše do zpráv o konfliktech atomové klíče jako Python tuple
# ("('A', 57, '', 'ND1')") - pro uživatele je převedeme na "A57:ND1".
_ATOM_KEY_REPR = re.compile(r"\('([^']*)', (-?\d+), '([^']*)', '([^']+)'\)")
_RESIDUE_KEY_REPR = re.compile(r"\('([^']*)', (-?\d+), '([^']*)'\)")


def _readable_conflict(message: str) -> str:
    message = _ATOM_KEY_REPR.sub(lambda m: f"{m[1]}{m[2]}{m[3]}:{m[4]}", message)
    return _RESIDUE_KEY_REPR.sub(lambda m: f"{m[1]}{m[2]}{m[3]}", message)


def _family_base_names(family_names: tuple) -> List[str]:
    """
    Názvy stavů bez koncové předpony (NHIE -> HIE), když ji mají všechny
    členy rodiny - uživatel volí HID/HIE/HIP, builder si koncovou variantu
    dohledá sám (viz forced_states v assign_protonation_states).
    """
    for prefix in ("N", "C"):
        if family_names and all(n.startswith(prefix) and len(n) > 3 for n in family_names):
            return [n[1:] for n in family_names]
    return list(family_names)


def titratable_residue_report(prot: Any) -> List[Dict[str, Any]]:
    """
    Všechna titrovatelná rezidua (dnes HIS) s tím, PROČ builder zvolil daný
    stav - pro lidskou kontrolu v Structure kroku a ruční přepsání v Expert
    režimu. Evidence se skládá z reportu builderu podle atomových míst
    rezidua: pevné vodíkové vazby (partner má jednoznačnou roli), vazby mezi
    dvěma proměnnými místy a nejednoznačné kontakty.
    """
    families = [tuple(names) for names, _default in prot.family_defaults]
    out = []
    for a in prot.assignments:
        rk = tuple(a.residue_key)

        def own(atom_key) -> bool:
            return tuple(atom_key[:3]) == rk

        evidence = []
        for e in prot.fixed_evidence:
            if own(e.variable_site):
                evidence.append({
                    "kind": "fixed",
                    "site": e.variable_site[3],
                    "role": e.required_role,
                    "partner": _format_atom_key(e.partner_site),
                    "distance_angstrom": round(e.distance_angstrom, 2),
                })
        for c in prot.variable_contacts:
            if own(c.site1) or own(c.site2):
                mine, other = (c.site1, c.site2) if own(c.site1) else (c.site2, c.site1)
                evidence.append({
                    "kind": "variable",
                    "site": mine[3],
                    "partner": _format_atom_key(other),
                    "distance_angstrom": round(c.distance_angstrom, 2),
                })
        for c in prot.ambivalent_contacts:
            if own(c.variable_site):
                evidence.append({
                    "kind": "ambivalent",
                    "site": c.variable_site[3],
                    "partner": _format_atom_key(c.partner_site),
                    "distance_angstrom": round(c.distance_angstrom, 2),
                })

        conflicts = [_readable_conflict(c.message) for c in prot.conflicts if any(own(site) for site in c.sites)]
        unevaluable = [f"{key[3]}: {reason}" for key, reason in prot.unevaluable_sites if own(key)]

        if a.is_forced:
            source = "user"
        elif any(e["kind"] in ("fixed", "variable") for e in evidence):
            source = "hbond"
        else:
            source = "ph_default"

        # Proč stav stojí za lidskou kontrolu (krok 3.5 v Expert režimu).
        # Ruční volba uživatele se znovu nekontroluje - už o ní rozhodl.
        review_reasons = []
        if source != "user":
            if conflicts:
                review_reasons.append("conflict")
            if any(e["kind"] == "ambivalent" for e in evidence):
                review_reasons.append("ambivalent")
            if unevaluable:
                review_reasons.append("unevaluable")
            if source == "ph_default":
                review_reasons.append("no_hbond")

        family = next((f for f in families if a.new_resname in f), (a.new_resname,))
        base = dict(zip(family, _family_base_names(family)))
        chain, resseq, icode = rk
        out.append({
            "residue": _format_residue_key(rk),
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "original": a.old_resname,
            "assigned": base.get(a.new_resname, a.new_resname),
            "default": base.get(a.default_resname, a.default_resname),
            "options": _family_base_names(family),
            "source": source,
            "evidence": evidence,
            "conflicts": conflicts,
            "unevaluable": unevaluable,
            "review_reasons": review_reasons,
        })
    return out


def build_histidine_review(prot: Any) -> Dict[str, Any]:
    """
    HIS sekce kroku 3.5 (Expert): titrovatelná rezidua s nejistým stavem
    (prázdný seznam = builder o všech rozhodl jednoznačně). Obecné problémy,
    které nepatří ke konkrétnímu reziduu, jdou zvlášť do `general_issues`.
    """
    rows = titratable_residue_report(prot)
    per_residue = {msg for r in rows for msg in r["conflicts"]}
    general = [
        msg for msg in (_readable_conflict(c.message) for c in prot.conflicts)
        if msg not in per_residue
    ]
    return {
        "histidines": [r for r in rows if r["review_reasons"]],
        "total_histidines": len(rows),
        "general_issues": general + list(prot.warnings),
    }


# Sekce kroku 3.5, které vyžadují rozhodnutí - když jsou všechny prázdné,
# příprava se nezastaví (general_issues samy o sobě ne).
REVIEW_SECTIONS = ("histidines", "zero_occupancy", "amide_flips", "removed_heterogens")


def build_preparation_summary(result: ForgePreparationResult) -> Dict[str, Any]:
    """
    Serializuje report objekty (state_assignment/crystal_ion_cleanup/solvation/
    ion_addition), které builder počítá při KAŽDÉ přípravě, ale API je dřív
    zahazovalo - vraceli jsme jen {message, warnings, validation}, report se
    jen logoval na server (viz run_workflow logger.info volání níže) a nikdy
    se nedostal na frontend. To přesně odpovídá review krokům W5.1/W5.2 v
    FORGE_general_design_v5.xlsx ("Review ... state assignments" a "Review
    ... retained crystal species, salt conditions and net charge") - ty tam
    nejsou jako samostatná W2/W3 review obrazovka (to spec explicitně
    nechce), ale musí být VIDĚT někde na konci přípravy.

    U protonation_assignments vracíme jen ne-defaultní přiřazení (ne
    kompletní inventář všech reziduí) - u velké struktury by kompletní seznam
    byl tisíce nezajímavých řádků, zatímco přesně tohle přiřazení "odlišné od
    výchozího" je to, co si review krok podle Excel dokumentu žádá zvýraznit.
    """
    summary: Dict[str, Any] = {}

    sa = result.state_assignment
    if sa is not None:
        covalent = sa.covalent
        summary["disulfide_bonds"] = [
            {
                "atom1": _format_atom_key(b.atom1),
                "atom2": _format_atom_key(b.atom2),
                "distance_angstrom": round(b.distance_angstrom, 3),
            }
            for b in covalent.bonds
        ]
        summary["covalent_state_changes"] = [
            {"residue": _format_residue_key(res_key), "from": old, "to": new}
            for res_key, old, new in covalent.state_changes
        ]
        summary["covalent_missing_bond_atoms"] = [
            {"residue": _format_residue_key(res_key), "resname": resname, "expected_atom": atom}
            for res_key, resname, atom in covalent.missing_bond_atoms
        ]

        prot = sa.protonation
        summary["protonation_assignments"] = [
            {
                "residue": _format_residue_key(a.residue_key),
                "default": a.default_resname,
                "assigned": a.new_resname,
            }
            for a in prot.assignments
            if not a.is_default
        ]
        summary["protonation_conflicts"] = [
            {"kind": c.kind, "message": _readable_conflict(c.message)} for c in prot.conflicts
        ]
        summary["protonation_warnings"] = list(prot.warnings)

    c = result.crystal_ion_cleanup
    if c is not None:
        summary["crystal_ion_cleanup"] = {
            "input_ions": c.input_ions,
            "removed_monovalent": c.removed_monovalent,
            "removed_nonstructural_multivalent": c.removed_nonstructural_multivalent,
            "retained_structural_monovalent": c.retained_structural_monovalent,
            "retained_structural_multivalent": c.retained_structural_multivalent,
            "replaced_by_magnesium": c.replaced_by_magnesium,
        }

    s = result.solvation
    if s is not None:
        summary["solvation"] = {
            "box_shape": s.box_shape,
            "padding_angstrom": s.padding_angstrom,
            "input_crystal_waters": s.input_crystal_waters,
            "retained_crystal_waters": s.retained_crystal_waters,
            "generated_waters": s.generated_waters,
            "total_waters": s.total_waters,
        }

    i = result.ion_addition
    if i is not None:
        summary["ion_addition"] = {
            "neutralization_ions": dict(i.neutralization_ions),
            "added_ions": dict(i.added_ions),
            "final_system_charge": round(i.final_system_charge, 4),
        }

    return summary
