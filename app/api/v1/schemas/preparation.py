"""
Request model of a structure preparation (POST /api/sidechains/start) and its
translation into the arguments the FORGE builder takes.
"""

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field

from app.core.exceptions import BadRequestError
from app.services.analysis import resolve_ion_mol_type
from app.services.ff_catalog_service import catalog_service


def _resolve_buildable_ion_mol_type(resname: str) -> Optional[str]:
    """
    Jako resolve_ion_mol_type(), ale navíc ověří, že pro ten resname
    existuje v aktuálním FF katalogu reálný force field (viz
    FFCatalogService.get_buildable_ion_resnames) - jinak by volba prošla
    validací požadavku, ale build by vždy spadl na KeyError hluboko v
    builderu. GET /validation/ions už nabízí jen tuhle buildable množinu,
    tohle je obrana pro případ zastaralého frontendu nebo přímého API volání.
    """
    mol_type = resolve_ion_mol_type(resname)
    if mol_type is None:
        return None
    buildable = catalog_service.get_buildable_ion_resnames()
    return mol_type if resname in buildable.get(mol_type, set()) else None


def _build_salt_specs(positive_ion: str, negative_ion: str, ionic_strength: float) -> List[Dict[str, Any]]:
    if ionic_strength <= 0:
        return []
    # mol_type se čte z converting_dictionary.json (přes resolve_ion_mol_type),
    # ne z lokálního hardcoded mapování - repo už jednou mělo tři nezávislé
    # kopie tohodle mapování a jejich nesoulad (Mg2+ nikde jako "Im") způsobil
    # pád na 1JJ2, viz docstring analysis.ff_requirements.required_ff_groups.
    cation_mol_type = _resolve_buildable_ion_mol_type(positive_ion)
    anion_mol_type = _resolve_buildable_ion_mol_type(negative_ion)
    if cation_mol_type is None:
        raise BadRequestError(f"Positive ion '{positive_ion}' is unknown or has no force field in the current catalog.")
    if anion_mol_type is None:
        raise BadRequestError(f"Negative ion '{negative_ion}' is unknown or has no force field in the current catalog.")
    return [{
        "cation": {"mol_type": cation_mol_type, "resname": positive_ion},
        "anion": {"mol_type": anion_mol_type, "resname": negative_ion},
        "concentration": ionic_strength,
    }]


# --- Data Models ---


class SaltStep(BaseModel):
    """Další sůl přidaná po té hlavní (např. NaCl + KCl) - viz additional_salts."""
    positive_ion: str = Field(..., description="Cation resname, see GET /validation/ions.")
    negative_ion: str = Field(..., description="Anion resname, see GET /validation/ions.")
    ionic_strength: float = Field(..., ge=0, description="Concentration (M) of this salt.")


class ProtonationOverride(BaseModel):
    chain: str
    resseq: int
    icode: str = ""
    state: str = Field(..., description="Target state name within the residue's family, e.g. HID/HIE/HIP.")


class ResidueDecision(BaseModel):
    chain: str
    resseq: int
    icode: str = ""
    apply: bool = Field(..., description="True = flip the amide / rebuild zero-occupancy atoms; False = keep as is.")


class StructureDecisionsPayload(BaseModel):
    amide_flips: List[ResidueDecision] = Field(default_factory=list)
    zero_occupancy: List[ResidueDecision] = Field(default_factory=list)
    acknowledged_heterogens: List[str] = Field(
        default_factory=list, description="Keys 'chain:resseq:icode:resname' of removed ligands the user saw."
    )


class PreparationRequest(BaseModel):
    workspace_id: str = Field(...)
    ff_selections: Dict[str, Any] = Field(
        ...,
        description="mol_type -> FF metadata z IDA API, stejný tvar jako /api/topology/.../generate. "
                    "Builder potřebuje silové pole už pro stavbu chybějících atomů, ne až pro topologii.",
    )
    ph: float = Field(7.0)
    crystal_water_mode: Literal["remove_all", "keep_water", "keep_all"] = Field("remove_all")
    add_solvent: bool = Field(False)
    box_padding_nm: float = Field(1.0)
    ionic_strength: float = Field(0.15, description="Salt concentration (M). The system is automatically neutralized.")
    positive_ion: str = Field(
        "Na+", description="Cation resname to add for neutralization/salt. See GET /validation/ions for valid options."
    )
    negative_ion: str = Field(
        "Cl-", description="Anion resname to add for neutralization/salt. See GET /validation/ions for valid options."
    )
    box_shape: Literal["cube", "octahedron", "truncated octahedron"] = Field("cube")
    clean_crystal_ions: bool = Field(
        True, description="Remove non-structural crystallographic ions (e.g. cryo/buffer ions) before solvation."
    )
    replace_structural_multivalent_with_mg: bool = Field(
        False,
        description="Replace retained structural multivalent ions (e.g. Zn2+, Ca2+) with Mg2+. "
                    "Requires a force field selection for the 'Im' ion group.",
    )
    concentration_mode: Literal["water_ratio", "box_volume"] = Field(
        "water_ratio", description="How ionic_strength is interpreted when placing salt ions."
    )
    protonation_overrides: List[ProtonationOverride] = Field(
        default_factory=list,
        description="Protonation states forced by the user (Expert mode), e.g. HIS -> HIE. "
                    "The builder skips its own decision for these residues.",
    )
    structure_decisions: StructureDecisionsPayload = Field(
        default_factory=StructureDecisionsPayload,
        description="User decisions from the Structure Check step (amide flips, zero-occupancy "
                    "rebuilds, acknowledged removed ligands).",
    )
    review_structure: bool = Field(
        False,
        description="Expert mode: stop before building (status 'structure_review' from "
                    "/sidechains/start) when something needs a manual check.",
    )
    auto_amide_flips: bool = Field(
        False,
        description="Guided/Standard mode: apply the suggested ASN/GLN amide flips without asking.",
    )
    additional_salts: List[SaltStep] = Field(
        default_factory=list,
        description="Further salts added after the primary one (positive_ion/negative_ion/ionic_strength). "
                    "Neutralization always uses the primary salt; builder places each salt in order.",
    )


def build_request_salt_specs(request: "PreparationRequest") -> List[Dict[str, Any]]:
    """
    Všechny soli z požadavku v pořadí: hlavní (positive_ion/negative_ion/
    ionic_strength) + additional_salts. Builder (add_ions_to_solvated_molecule)
    je umisťuje postupně a neutralizuje podle první z nich. Kroky s nulovou
    koncentrací se vynechají (stejně jako dřív u hlavní soli).
    """
    specs = _build_salt_specs(request.positive_ion, request.negative_ion, request.ionic_strength)
    for step in request.additional_salts:
        specs += _build_salt_specs(step.positive_ion, step.negative_ion, step.ionic_strength)
    return specs


# Popisky ze staršího PDBFixer-orientovaného UI na SolvationSettings.box_shape,
# které builder skutečně zná (viz SUPPORTED_BOX_SHAPES ve forge_molecule_solvation.py).
BOX_SHAPE_MAP = {
    "cube": "cubic",
    "octahedron": "truncated_octahedron",
    "truncated octahedron": "truncated_octahedron",
}
