import json
import logging
import traceback
import time
import aiofiles
from typing import Any, Dict, List, Literal

from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field

from app.core.exceptions import AppBaseException, InternalError
from app.services.validation.service import ValidationService
from app.services.structure.forge_service import ForgeStructureService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
validation_service = ValidationService()
forge_service = ForgeStructureService()

# Mapování dnešních jednoduchých pojmenovaných iontů na FORGE mol_type - všechny
# monovalentní ionty v Literal seznamu níže jsou v converting_dictionary vedené
# pod "I1" (Na+/K+/Cl-) nebo "I1+" (Li+/Rb+/Cs+/F-/Br-/I-).
_MONOVALENT_ION_MOL_TYPE = {
    "Na+": "I1", "K+": "I1", "Cl-": "I1",
    "Li+": "I1+", "Cs+": "I1+", "Rb+": "I1+", "F-": "I1+", "Br-": "I1+", "I-": "I1+",
}


def _build_salt_specs(positive_ion: str, negative_ion: str, ionic_strength: float) -> List[Dict[str, Any]]:
    if ionic_strength <= 0:
        return []
    return [{
        "cation": {"mol_type": _MONOVALENT_ION_MOL_TYPE.get(positive_ion, "I1"), "resname": positive_ion},
        "anion": {"mol_type": _MONOVALENT_ION_MOL_TYPE.get(negative_ion, "I1"), "resname": negative_ion},
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
    positive_ion: Literal["Na+", "K+", "Li+", "Cs+", "Rb+"] = Field("Na+")
    negative_ion: Literal["Cl-", "F-", "Br-", "I-"] = Field("Cl-")
    box_shape: Literal["cube", "octahedron", "truncated octahedron"] = Field("cube")


# --- Endpoints ---

@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Asynchronně načte PDB soubor a provede úvodní analýzu struktury.
    Detekuje alternativní lokace, chybějící atomy a kompatibilitu s forcefieldem.
    """
    workspace_manager.require_workspace(request.workspace_id)

    try:
        pdb_path = workspace_manager.get_file_path(request.workspace_id)

        # Non-blocking I/O pro čtení souboru
        async with aiofiles.open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Delegace CPU-bound validace do threadpoolu pro zamezení blokování event loopu
        return await run_in_threadpool(
            validation_service.validate_pdb_content,
            pdb_content,
            request.label
        )

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
            "validation": validation_results
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