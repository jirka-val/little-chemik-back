import json
import logging
import traceback
import time
import aiofiles
from typing import Any, Dict, List, Literal, Optional

from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.core.exceptions import AppBaseException, BadRequestError, InternalError
from app.services.analysis_service import list_ion_options, resolve_ion_mol_type
from app.services.ff_catalog_service import catalog_service
from app.services.validation.service import ValidationService
from app.services.structure.forge_service import ForgeStructureService, build_preparation_summary
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
validation_service = ValidationService()
forge_service = ForgeStructureService()


def _resolve_buildable_ion_mol_type(resname: str) -> Optional[str]:
    """
    Jako resolve_ion_mol_type(), ale navíc ověří, že pro ten resname
    existuje v aktuálním FF katalogu reálný force field (viz
    FFCatalogService.get_buildable_ion_resnames) - jinak by volba prošla
    /prepare validací, ale build by vždy spadl na KeyError hluboko v
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
    # pád na 1JJ2, viz docstring analysis_service.required_ff_groups.
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

class ValidationRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    label: str = Field("molecule_from_front", description="Molecule state identifier")


class FixAltLocRequest(BaseModel):
    workspace_id: str = Field(..., description="Workspace ID containing the molecule")
    selections: Dict[str, str] = Field(
        ...,
        description="Variant selection map, e.g., {'A-42': 'B'}",
        example={"A-42": "B", "A-15": "A"}
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


# --- Endpoints ---

@router.get("/ions", summary="Vrátí vybíratelné ionty pro salt/neutralizaci po mol_type skupině")
async def get_ion_options():
    """
    Zdroj pravdy pro pos./neg. ion dropdowny na frontendu - converting_dictionary.json
    (chemická identita, ne z hardcoded seznamu - viz poznámka u
    required_ff_groups o třech nezávislých kopiích, co tenhle seznam dřív
    měl a jak to skončilo) PROTNUTÉ s FFCatalogService.get_buildable_ion_resnames()
    (jestli pro ten ion v aktuálním katalogu vůbec existuje force field).
    Bez tohohle průniku by šlo z dropdownu vybrat ion (typicky exotický
    lanthanid/aktinid z Im+, nebo nábojem-sufixovanou variantu tam, kde
    katalog zná jen starší holý symbol), který by vždy skončil neřešitelnou
    chybou - buď 409 "missing force field" (nikde není co vybrat), nebo
    KeyError hluboko v builderu. "I1"/"I1+" jsou monovalentní, "Im"/"Im+"
    dvojmocné+ (Im obsahuje Mg2+).
    """
    known = list_ion_options()
    buildable = catalog_service.get_buildable_ion_resnames()
    return {
        mol_type: sorted(set(resnames) & buildable.get(mol_type, set()))
        for mol_type, resnames in known.items()
    }


@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Asynchronně načte PDB soubor a provede úvodní analýzu struktury.
    Detekuje alternativní lokace, chybějící atomy a kompatibilitu s forcefieldem.
    """
    workspace_manager.require_workspace(request.workspace_id)
    logger.info(f"Validating structure for workspace {request.workspace_id} (label: {request.label})...")

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)

        # Non-blocking I/O pro čtení souboru
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Delegace CPU-bound validace do threadpoolu pro zamezení blokování event loopu
        result = await run_in_threadpool(
            validation_service.validate_pdb_content,
            pdb_content,
            request.label
        )

        summary = result.get("summary", {})
        analysis = result.get("analysis", {})
        logger.info(
            f"Validation finished for workspace {request.workspace_id}: "
            f"ready_for_hpc={summary.get('is_ready_for_hpc')}, "
            f"errors={len(analysis.get('errors', []))}, warnings={len(analysis.get('warnings', []))}"
        )
        return result

    except Exception as e:
        logger.exception(f"Validation error for workspace {request.workspace_id}: {str(e)}")
        raise InternalError(f"Validation error: {str(e)}")


@router.post("/preview-selection", summary="Náhled geometrie před aplikací")
async def preview_selection(request: FixAltLocRequest):
    """
    Provede náhled geometrických úprav na základě uživatelských selekcí
    bez trvalého zápisu do souboru.
    """
    workspace_manager.require_workspace(request.workspace_id)

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Izolovaný výpočet kontinuity řetězce
        issues = await run_in_threadpool(
            validation_service.conf_manager.validate_continuity,
            pdb_content,
            request.selections
        )

        is_safe = len(issues) == 0

        return {
            "is_ok": is_safe,
            "issues": issues,
            "message": "Selection is geometrically valid" if is_safe else "Critical chain gaps detected"
        }

    except Exception as e:
        logger.exception(f"Selection preview error for workspace {request.workspace_id}")
        raise InternalError(str(e))


@router.post("/apply-selections", summary="Aplikuje výběr konformací a ověří kontinuitu")
async def apply_selections(request: FixAltLocRequest):
    """
    Aplikuje vybrané konformace a asynchronně přepíše zdrojový PDB soubor.
    """
    workspace_manager.require_workspace(request.workspace_id)
    logger.info(
        f"Applying {len(request.selections)} AltLoc selection(s) for workspace "
        f"{request.workspace_id} (overwrites structure.pdb)..."
    )

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        result = await run_in_threadpool(
            validation_service.apply_alt_loc_selection,
            pdb_content,
            request.selections
        )

        # Zápis upravené struktury zpět na disk (non-blocking)
        if "pdb_content" in result:
            async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
                await f.write(result["pdb_content"])

        logger.info(f"AltLoc selections applied and written for workspace {request.workspace_id}.")
        return result

    except Exception as e:
        logger.exception(f"Error applying selections for workspace {request.workspace_id}")
        raise InternalError(f"Error applying selections: {str(e)}")


# Popisky ze staršího PDBFixer-orientovaného UI na SolvationSettings.box_shape,
# které builder skutečně zná (viz SUPPORTED_BOX_SHAPES ve forge_molecule_solvation.py).
_BOX_SHAPE_MAP = {
    "cube": "cubic",
    "octahedron": "truncated_octahedron",
    "truncated octahedron": "truncated_octahedron",
}


@router.post("/prepare", summary="Kompletní příprava: Protonace, Solvatace, Ionty")
async def prepare_molecule(request: PreparationRequest):
    """
    Spouští výpočetně náročný proces přípravy struktury přes FORGE builder
    (stavy/protonace, stavba chybějících atomů, solvatace, ionty).
    Endpoint využívá threadpool pro paralelizaci a zachování odezvy serveru.
    """
    start_time = time.time()
    logger.info(f"Started molecule preparation for workspace: {request.workspace_id} (pH: {request.ph})")

    workspace_manager.require_workspace(request.workspace_id)

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        logger.info(f"Starting ForgeStructureService (solvent: {request.add_solvent}, ions: {request.ionic_strength})")

        # Spuštění primární chemické transformace v dedikovaném vlákně
        result = await run_in_threadpool(
            forge_service.prepare_structure,
            pdb_text=pdb_content,
            ff_selections=request.ff_selections,
            ph=request.ph,
            add_solvent_and_ions=request.add_solvent,
            salts=_build_salt_specs(request.positive_ion, request.negative_ion, request.ionic_strength),
            box_shape=_BOX_SHAPE_MAP.get(request.box_shape),
            box_padding_angstrom=request.box_padding_nm * 10.0,
            keep_crystal_waters=request.crystal_water_mode != "remove_all",
            crystal_water_mode=request.crystal_water_mode,
            clean_crystal_ions=request.clean_crystal_ions,
            replace_structural_multivalent_with_mg=request.replace_structural_multivalent_with_mg,
            concentration_mode=request.concentration_mode,
        )

        if not result.pdb_text:
            logger.error("ForgeStructureService returned empty output.")
            raise InternalError("The resulting PDB content is empty.")

        logger.info(f"Structure prepared successfully in {time.time() - start_time:.2f}s. Starting write process.")

        async with aiofiles.open(pdb_path, "w", encoding="utf-8") as f:
            await f.write(result.pdb_text)

        # Sidecar s autoritativním ff_resname/group - viz TopologyService._load_forge_meta.
        meta_path = pdb_path.with_name(pdb_path.name.replace(".pdb", ".forge_meta.json"))
        async with aiofiles.open(meta_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(result.forge_meta))

        # Finální validace integrity po modifikaci
        validation_results = await run_in_threadpool(
            validation_service.validate_pdb_content,
            result.pdb_text,
            label="prepared_state"
        )

        return {
            "message": "Structure successfully prepared.",
            "warnings": result.warnings,
            "validation": validation_results,
            "preparation_summary": build_preparation_summary(result),
        }

    except AppBaseException:
        # ForgeMissingDOFError apod. - necháme propadnout ke globálnímu
        # app_exception_handler (viz app/main.py), který jim dá jejich vlastní
        # status_code (409 pro missing_dof) místo generické 500 níže.
        raise
    except Exception as e:
        logger.error(
            f"Critical failure during preparation of workspace {request.workspace_id}:\n{traceback.format_exc()}")
        raise InternalError(str(e))