"""
Citace force fieldů pro report a popis postupu.

Katalog IDA má u FF jen seznam DOI a autory datasetu, a to ne vždy správně.
U RNA FF proto citujeme podle tabulky od autorů FF (ff.xlsx, 10/2026): každý
FF celou řadou článků, ze kterých vychází, v pořadí od nejstaršího. FF, které
tabulka nepokrývá, se citují z katalogu jako dřív.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

# klíč -> (autor a rok, název, DOI)
ARTICLES: Dict[str, Tuple[str, str, str]] = {
    "cornell1995": ("Cornell et al. (1995)",
                    "A Second Generation Force Field for the Simulation of Proteins, Nucleic Acids, and Organic "
                    "Molecules", "10.1021/ja00124a002"),
    "wang2000": ("Wang et al. (2000)",
                 "How Well Does a Restrained Electrostatic Potential (RESP) Model Perform in Calculating "
                 "Conformational Energies of Organic and Biological Molecules?",
                 "10.1002/1096-987X(200009)21:12<1049::AID-JCC3>3.0.CO;2-F"),
    "perez2007": ("Pérez et al. (2007)",
                  "Refinement of the AMBER Force Field for Nucleic Acids: Improving the Description of α/γ "
                  "Conformers", "10.1529/biophysj.106.097782"),
    "zgarbova2011": ("Zgarbová et al. (2011)",
                     "Refinement of the Cornell et al. Nucleic Acids Force Field Based on Reference Quantum Chemical "
                     "Calculations of Glycosidic Torsion Profiles", "10.1021/ct200162x"),
    "yildirim2010": ("Yildirim et al. (2010)",
                     "Reparameterization of RNA χ Torsion Parameters for the AMBER Force Field and Comparison to NMR "
                     "Spectra for Cytidine and Uridine", "10.1021/ct900604a"),
    "steinbrecher2012": ("Steinbrecher et al. (2012)", "Revised AMBER Parameters for Bioorganic Phosphates",
                         "10.1021/ct300613v"),
    "bergonzo2015": ("Bergonzo and Cheatham (2015)",
                     "Improved Force Field Parameters Lead to a Better Description of RNA Structure",
                     "10.1021/acs.jctc.5b00444"),
    "aytenfisu2017": ("Aytenfisu et al. (2017)",
                      "Revised RNA Dihedral Parameters for the Amber Force Field Improve RNA Molecular Dynamics",
                      "10.1021/acs.jctc.6b00870"),
    "tan2018": ("Tan et al. (2018)", "RNA Force Field with Accuracy Comparable to State-of-the-Art Protein Force Fields",
                "10.1073/pnas.1713027115"),
    "kuhrova2019": ("Kührová et al. (2019)",
                    "Improving the Performance of the Amber RNA Force Field by Tuning the Hydrogen-Bonding Interactions",
                    "10.1021/acs.jctc.8b00955"),
    "frohlking2022": ("Fröhlking et al. (2022)",
                      "Automatic Learning of Hydrogen-Bond Fixes in the AMBER RNA Force Field",
                      "10.1021/acs.jctc.2c00200"),
    "mlynsky2023": ("Mlýnský et al. (2023)",
                    "Simple Adjustment of Intranucleotide Base-Phosphate Interaction in the OL3 AMBER Force Field "
                    "Improves RNA Simulations", "10.1021/acs.jctc.3c00990"),
    "raguette2024": ("Raguette et al. (2024)",
                     "Adjusting the Energy Profile for CH–O Interactions Leads to Improved Stability of RNA Stem-Loop "
                     "Structures in MD Simulations", "10.1021/acs.jpcb.4c01910"),
    "mlynsky2026": ("Mlýnský et al. (2026)",
                    "From gHBfix to NBfix: Reweighting-Driven Refinement of Hydrogen-Bond Interactions in RNA Force "
                    "Fields", "10.1021/acs.jctc.6c00924"),
}

_OL3 = ("cornell1995", "wang2000", "perez2007", "zgarbova2011")
_OL3CP = _OL3 + ("steinbrecher2012",)

# ff_name (bez ohledu na velikost písmen) -> citované články
FF_REFERENCES: Dict[str, Tuple[str, ...]] = {
    "OL3": _OL3,
    "OL3CP": _OL3CP,
    "OL3-LJBB": _OL3CP + ("bergonzo2015",),
    "OL3-HR2.7": _OL3 + ("raguette2024",),
    "OL3CP_gHBfix19": _OL3CP + ("kuhrova2019",),
    "OL3CP_gHBfix21": _OL3CP + ("frohlking2022",),
    "OL3CP_NBfix19": _OL3CP + ("kuhrova2019", "mlynsky2026"),
    "OL3CP_NBfix-0BPh": _OL3CP + ("mlynsky2023",),
    "OL3CP_gHBfix19_NBfix-0BPh": _OL3CP + ("kuhrova2019", "mlynsky2023"),
    "OL3CP_gHBfix21_NBfix-0BPh": _OL3CP + ("frohlking2022", "mlynsky2023"),
    "FF99BSC0-YIL": ("cornell1995", "wang2000", "perez2007", "yildirim2010"),
    "ROC": ("aytenfisu2017",),
    "DESRES": ("tan2018",),
}
_BY_NAME = {name.upper(): keys for name, keys in FF_REFERENCES.items()}


def format_reference(key: str) -> str:
    label, title, doi = ARTICLES[key]
    return f"{label}. {title.rstrip('.')}. doi:{doi}."


def _catalog_references(ff: Dict[str, Any]) -> List[str]:
    """Jedna citace z katalogu IDA: autoři datasetu a všechna DOI."""
    dois = [d.strip() for d in str(ff.get("reference_article_doi") or "").split(",") if d.strip()]
    authors = ff.get("creators") or []
    if not dois and not authors:
        return []
    parts = []
    if authors:
        parts.append("; ".join(authors[:6]).rstrip(".") + (" et al" if len(authors) > 6 else ""))
    if dois:
        parts.append(", ".join(f"doi:{d}" for d in dois))
    return [". ".join(parts) + "."]


def ff_references(ff: Dict[str, Any]) -> List[str]:
    """Citace FF podle tabulky, jinak z katalogu; prázdné, když není co citovat."""
    for name in (ff.get("ff_name"), ff.get("display_name")):
        keys = _BY_NAME.get(str(name or "").upper())
        if keys:
            return [format_reference(k) for k in keys]
    return _catalog_references(ff)
