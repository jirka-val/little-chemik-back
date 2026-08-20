import logging
from typing import Any, Dict, List

import aiofiles
from fastapi import APIRouter
from fastapi.concurrency import run_in_threadpool

from app.core.exceptions import InternalError
from app.services.analysis_service import required_ff_groups
from app.services.forcefield_service import ForceFieldService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
ff_validator = ForceFieldService()


@router.get("/{workspace_id}", summary="Získá dostupné forcefieldy pro danou molekulu")
async def get_my_forcefields(workspace_id: str):
    """
    Asynchronně načte PDB soubor a v dedikovaném vlákně zanalyzuje chemické
    složení (typ reziduí, přítomné ionty) pro detekci kompatibilních
    forcefieldů.

    Dřív se tu ionty vždy paušálně nabízely přes obecné "I" (viz
    ForceFieldService.get_matching_forcefields), které se rozbalilo na
    VŠECHNY čtyři iontové mol_type podskupiny (I1/I1+/Im/Im+) najednou -
    frontend tak dostal jeden nerozlišený seznam ~15 iontových FF bez
    ponětí, které konkrétní podskupiny struktura reálně potřebuje. To byl
    přímý zdroj pádu na 1JJ2 (uživatel vybral I1+ místo potřebných I1+Im,
    protože nikde nebylo vidět, že by měl vybírat dvě různé položky).
    required_ff_groups teď zjistí přesně, jaké mol_type skupiny (podle
    reálně přítomných polymerů/iontů) jsou potřeba, takže se dají FF
    seskupit podle skupiny a chybějící pokrytí je vidět dopředu, ne až po
    pádu v /prepare.
    """
    workspace_manager.require_workspace(workspace_id)

    try:
        # Načteme PDB ze souboru ASYNCHRONNĚ
        path = workspace_manager.get_file_path(workspace_id)
        async with aiofiles.open(path, "r", encoding="utf-8") as f:
            pdb_content = await f.read()

        # Přesné mol_type skupiny, které tahle konkrétní struktura
        # potřebuje - viz docstring výše. Solvataci/ionty nabízíme vždy
        # jako volitelnou (add_solvent_and_ions=True), i když ji uživatel
        # nakonec nepoužije - FF pro ně je potřeba vybrat předem.
        required = await run_in_threadpool(required_ff_groups, pdb_content, True, None)

        search_types = set(required.keys())
        search_types.add("W")  # obecné "W" -> get_matching_forcefields rozbalí na W3/W4/W5

        ffs = await run_in_threadpool(ff_validator.get_matching_forcefields, sorted(search_types))

        ffs_by_group: Dict[str, List[Any]] = {}
        for ff in ffs:
            for mol_type in ff.get("molecule_type") or []:
                ffs_by_group.setdefault(mol_type, []).append(ff)

        return {
            "detected_types": sorted(required.keys()),
            "required_groups": required,
            "forcefields_by_group": ffs_by_group,
            # Zachováno pro zpětnou kompatibilitu s dosavadním plochým seznamem.
            "forcefields": ffs,
        }

    except Exception as e:
        logger.exception(f"Error fetching forcefields for workspace {workspace_id}")
        raise InternalError(f"Internal server error while fetching forcefields: {str(e)}")