"""Which force-field groups (solute, water, ions) a structure needs, and
which ions can be selected."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

from app.utils.adams4sims_processing_library.utils.alias import resn_alias

from .residues import load_converting_dictionary
from .sequence import build_sequence_tokens

_POLYMER_GROUP_LABELS = {"P": "protein", "R": "RNA", "D": "DNA"}


def resolve_ion_mol_type(resname: str, conv: Optional[dict] = None) -> Optional[str]:
    """
    Jediné povolené místo pro převod iontového resname na FORGE mol_type
    (I1/I1+/Im/Im+) - viz varování v required_ff_groups o třech nezávislých
    hardcoded kopiích tohohle mapování, které dřív v repu existovaly a
    způsobily reálný pád (KeyError na 1JJ2, Mg2+ nerozpoznané jako Im).
    Vrací None, pokud converting_dictionary žádný takový iont nezná.
    """
    conv = conv if conv is not None else load_converting_dictionary()
    for mol_type in ("I1", "I1+", "Im", "Im+"):
        if resname in conv.get(mol_type, {}):
            return mol_type
    return None


def list_ion_options() -> Dict[str, List[str]]:
    """
    Všechny iontové resnames, které converting_dictionary.json zná po
    mol_type skupině - tohle je čistě chemická identita (jaký mol_type ion
    má), NE záruka, že pro něj v aktuální FF katalogové sadě existuje reálný
    force field. Kdo potřebuje jen skutečně stavitelné ionty (typicky
    dropdown na frontendu), musí výsledek protnout s
    FFCatalogService.get_buildable_ion_resnames() - viz GET
    /api/validation/ions. Bez toho průniku by šlo vybrat ion, u kterého
    builder vždy spadne na KeyError (viz docstring
    get_buildable_ion_resnames pro konkrétní příklady).
    """
    conv = load_converting_dictionary()
    return {
        mol_type: sorted(conv.get(mol_type, {}).keys())
        for mol_type in ("I1", "I1+", "Im", "Im+")
    }


def required_ff_groups(
    pdb_text: str,
    add_solvent_and_ions: bool = True,
    salts: Optional[List[Dict[str, Any]]] = None,
    require_mg: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """
    Zjistí, jaké FORGE mol_type skupiny (P/R/D/W/I1/I1+/Im/Im+) tahle
    konkrétní struktura reálně potřebuje, ať se to dá zkontrolovat proti
    ff_selections ještě PŘED spuštěním (drahého, u velkých struktur i
    několikaminutového) buildeu - viz forge_service.prepare_structure.

    Jediný zdroj pravdy je converting_dictionary.json (stejný, jaký uvnitř
    používá i builder), ne samostatný hardcoded seznam iontů - takových už
    v repu dřív byly nezávislé kopie (pdb_service.get_molecule_types a
    validation._MONOVALENT_ION_MOL_TYPE, ten druhý od teď nahrazený
    resolve_ion_mol_type() níže) a právě jejich vzájemná neshoda (žádná
    neznala "Im") byla přímou příčinou pádu "KeyError: Ion parameters
    missing for Im:Mg2+" na 1JJ2 -
    uživatel vybral I1+ místo Im a nikde nebylo vidět, že Im je potřeba.

    Ionty "Im"/"Im+" jsou v konverzním slovníku pojmenované matoucně -
    navzdory "m" v názvu jde o DVOJMOCNÉ KATIONTY (Im obsahuje Mg2+), ne o
    aniony. Nespoléhat na název, vždy číst přímo ze slovníku.
    """
    conv = load_converting_dictionary()
    result: Dict[str, Dict[str, Any]] = {}

    sequence_data = build_sequence_tokens(pdb_text, chain=None, fill_gaps=True)
    polymer_examples: Dict[str, List[str]] = {}
    for chain_data in sequence_data.get("chains", {}).values():
        for token in chain_data.get("tokens", []):
            group = token.get("group")
            if group in _POLYMER_GROUP_LABELS and not token.get("is_gap"):
                examples = polymer_examples.setdefault(group, [])
                name = token.get("pdb_resname")
                if name and name not in examples and len(examples) < 5:
                    examples.append(name)

    for group, examples in polymer_examples.items():
        result[group] = {
            "reason": f"{_POLYMER_GROUP_LABELS[group]} residues present (e.g. {', '.join(examples)})",
        }

    if add_solvent_and_ions:
        result["W"] = {"reason": "solvation requested (add_solvent_and_ions=True)"}

        ion_mol_type = {}
        for mol_type in ("I1", "I1+", "Im", "Im+"):
            for resname in conv.get(mol_type, {}):
                ion_mol_type[resname] = mol_type

        ion_counts: Dict[str, int] = {}
        for line in pdb_text.splitlines():
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            # Konverzní slovník vede ionty pod jejich kanonickým FF názvem
            # (Mg2+, Na+, ...), ne pod raw PDB zkratkou (MG, NA, ...) -
            # stejný alias krok jako u polymerních reziduí v _infer_group.
            canonical = resn_alias(resname)
            if canonical in ion_mol_type:
                ion_counts[canonical] = ion_counts.get(canonical, 0) + 1
            elif resname in ion_mol_type:
                ion_counts[resname] = ion_counts.get(resname, 0) + 1

        by_group: Dict[str, List[str]] = {}
        for resname, count in ion_counts.items():
            by_group.setdefault(ion_mol_type[resname], []).append(f"{resname} (x{count})")

        salt_mol_types: Set[str] = set()
        if salts:
            for spec in salts:
                for side in ("cation", "anion"):
                    mol_type = (spec.get(side) or {}).get("mol_type")
                    if mol_type:
                        salt_mol_types.add(mol_type)
        else:
            # Prázdný/chybějící seznam solí = jen výchozí neutralizace
            # K+/Cl- (viz INTEGRATION_CONTRACT.md "Salt input") - ta vždy
            # potřebuje I1, i když struktura sama žádné krystalové ionty
            # nemá.
            salt_mol_types.add("I1")

        for mol_type in salt_mol_types:
            by_group.setdefault(mol_type, []).append("default/requested neutralization or salt")

        if require_mg:
            by_group.setdefault("Im", []).append("structural multivalent ion -> Mg2+ replacement requested")

        for mol_type, reasons in by_group.items():
            result[mol_type] = {"reason": f"ions: {', '.join(reasons)}"}

    return result
