import logging
from typing import Any, Dict, List

import aiofiles
from fastapi import APIRouter, Header
from fastapi.concurrency import run_in_threadpool

from app.core.config import settings
from app.core.exceptions import AppBaseException, BadRequestError, ExternalServiceError, ForbiddenError, InternalError
from app.services.analysis_service import required_ff_groups
from app.services.ff_catalog_service import catalog_service
from app.services.ff_classification_service import classification_service
from app.services.forcefield_service import ForceFieldService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

router = APIRouter()
ff_validator = ForceFieldService()

SOLUTE_MOL_TYPES = ("D", "R", "P")
ION_MOL_TYPES = ("I1", "I1+", "Im", "Im+")
_TIER_SORT_ORDER = {"recommended": 0, "supported": 1, "obsolete": 2, "new_unclassified": 3}
_MODES = ("guided", "standard", "expert")


def _classify_and_sort(ffs: List[Dict[str, Any]], classify_fn, mode: str) -> List[Dict[str, Any]]:
    """
    Obohatí každý FF objekt o "tier"/"is_default" a podle módu vyfiltruje a
    seřadí:
      - guided:   jen FF označený jako default (panel jen potvrzuje výběr),
      - standard: recommended + supported, v tomto pořadí, obsolete/
                   new_unclassified skryté,
      - expert:   úplně vše (recommended/supported/obsolete/new_unclassified),
                   nic není preferováno - uživatel musí vybrat sám.
    """
    enriched = []
    for ff in ffs:
        name = ff_validator.ff_name(ff)
        info = classify_fn(name)
        enriched.append({**ff, "tier": info["tier"], "is_default": info["is_default"]})

    if mode == "guided":
        enriched = [f for f in enriched if f["is_default"]]
    elif mode == "standard":
        enriched = [f for f in enriched if f["tier"] in ("recommended", "supported")]

    enriched.sort(key=lambda f: _TIER_SORT_ORDER.get(f["tier"], 99))
    return enriched


def _build_water_profiles(ffs_by_group: Dict[str, List[Any]], mode: str) -> List[Dict[str, Any]]:
    """
    Nahrazuje dřívější plochý seznam W3/W4/W5 FF "profily" ze
    force_fields.json - svazkem konkrétní vody s jejími kompatibilními ionty
    pro I1/I1+/Im/Im+ (viz plán, sekce Profily vody a iontů). Profil se
    nabídne, jen pokud je jeho vodní FF vůbec v aktuálním katalogu; resolvnuté
    ionty jsou jen ty, které katalog reálně obsahuje (chybějící se v ion_ffs
    prostě vynechají, frontend to nemá bránit ve výběru vody samotné).
    """
    water_by_name: Dict[str, Any] = {}
    for wt in ("W3", "W4", "W5"):
        for ff in ffs_by_group.get(wt, []):
            water_by_name[ff_validator.ff_name(ff)] = ff

    ion_by_name: Dict[str, Dict[str, Any]] = {
        mt: {ff_validator.ff_name(ff): ff for ff in ffs_by_group.get(mt, [])} for mt in ION_MOL_TYPES
    }

    profiles = []
    for profile in classification_service.all_water_profiles():
        water_ff = water_by_name.get(profile["water"])
        if not water_ff:
            continue

        ion_ffs = {}
        for mol_type, ion_name in (profile.get("ions") or {}).items():
            ion_ff = ion_by_name.get(mol_type, {}).get(ion_name)
            if ion_ff:
                ion_ffs[mol_type] = ion_ff

        profiles.append({**profile, "water_ff": water_ff, "ion_ffs": ion_ffs})

    if mode == "guided":
        profiles = [p for p in profiles if p["is_default"]]
    elif mode == "standard":
        profiles = [p for p in profiles if p["tier"] in ("recommended", "supported")]

    profiles.sort(key=lambda p: _TIER_SORT_ORDER.get(p["tier"], 99))
    return profiles


@router.get("/{workspace_id}", summary="Získá dostupné forcefieldy pro danou molekulu")
async def get_my_forcefields(workspace_id: str, mode: str = "standard"):
    """
    Vrací FF dostupné pro tuhle strukturu, seskupené a obohacené o tier
    (recommended/supported/obsolete/new_unclassified) podle data/force_fields.json,
    filtrované podle `mode` (guided/standard/expert - viz _classify_and_sort).

    Na rozdíl od dřívější verze NEVOLÁ IDA přímo - čte z lokálního katalogového
    snapshotu (FFCatalogService), který se obnovuje buď ručně
    (POST /catalog/refresh), nebo automaticky jednou denně na pozadí. To byl
    přímý zdroj pomalosti tohohle endpointu předtím.

    Dřív se tu ionty vždy paušálně nabízely přes obecné "I", které se
    rozbalilo na VŠECHNY čtyři iontové mol_type podskupiny (I1/I1+/Im/Im+)
    najednou - frontend tak dostal jeden nerozlišený seznam ~15 iontových FF
    bez ponětí, které konkrétní podskupiny struktura reálně potřebuje. To byl
    přímý zdroj pádu na 1JJ2 (uživatel vybral I1+ místo potřebných I1+Im,
    protože nikde nebylo vidět, že by měl vybírat dvě různé položky).
    required_ff_groups teď zjistí přesně, jaké mol_type skupiny (podle
    reálně přítomných polymerů/iontů) jsou potřeba, takže se dají FF
    seskupit podle skupiny a chybějící pokrytí je vidět dopředu, ne až po
    pádu v /prepare.
    """
    if mode not in _MODES:
        raise BadRequestError(f"Unknown mode '{mode}', expected one of {_MODES}.")

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
        search_types.add("W")  # obecné "W" -> filter_forcefields rozbalí na W3/W4/W5

        all_ffs = catalog_service.get_forcefields()
        ffs = ForceFieldService.filter_forcefields(all_ffs, sorted(search_types))

        ffs_by_group: Dict[str, List[Any]] = {}
        for ff in ffs:
            for mol_type in ff.get("molecule_type") or []:
                ffs_by_group.setdefault(mol_type, []).append(ff)

        enriched_by_group: Dict[str, List[Any]] = {}
        for group_key in required.keys():
            if group_key in SOLUTE_MOL_TYPES:
                classify_fn = lambda name, g=group_key: classification_service.classify_solute(g, name)
            elif group_key in ION_MOL_TYPES:
                classify_fn = lambda name, g=group_key: classification_service.classify_ion(g, name)
            else:
                continue  # "W" se řeší samostatně níže jako water_profiles
            enriched_by_group[group_key] = _classify_and_sort(ffs_by_group.get(group_key, []), classify_fn, mode)

        water_profiles = _build_water_profiles(ffs_by_group, mode) if "W" in required else []

        return {
            "mode": mode,
            "catalog_fetched_at": catalog_service.fetched_at(),
            "detected_types": sorted(required.keys()),
            "required_groups": required,
            "forcefields_by_group": enriched_by_group,
            "water_profiles": water_profiles,
            # Zachováno pro zpětnou kompatibilitu s dosavadním plochým seznamem.
            "forcefields": ffs,
        }

    except AppBaseException:
        raise
    except Exception as e:
        logger.exception(f"Error fetching forcefields for workspace {workspace_id}")
        raise InternalError(f"Internal server error while fetching forcefields: {str(e)}")


@router.post("/catalog/refresh", summary="Ručně obnoví lokální katalog force fieldů z IDA")
async def refresh_catalog():
    """
    Spouští se z tlačítka "Refresh" ve FF panelu. Dělá přesně to, co jinak
    dělá noční background job (viz ff_catalog_refresher.py) - stáhne čerstvý
    seznam z IDA, uloží ho jako lokální snapshot a zařadí nově objevené FF do
    new_unclassified. Není admin-gated - jde jen o ruční spuštění téže
    neškodné, idempotentní operace, kterou beztak dělá noční job automaticky.
    """
    try:
        snapshot = await run_in_threadpool(catalog_service.refresh_catalog)
    except Exception as e:
        logger.exception("Manual FF catalog refresh failed")
        raise ExternalServiceError(f"Nepodařilo se obnovit katalog force fieldů z IDA: {e}")

    return {"fetched_at": snapshot["fetched_at"], "count": len(snapshot["forcefields"])}


@router.patch("/classification", summary="(admin) Přeřadí force field (D/R/P) do jiného tieru")
async def patch_classification(payload: Dict[str, str], x_admin_token: str = Header(default="")):
    """
    Ruční zařazení FF čekajícího v new_unclassified (nebo přeřazení
    existujícího) do recommended/supported/obsolete. Chráněno sdíleným
    tokenem (Settings.ADMIN_TOKEN) - dokud není v .env nastavený, endpoint je
    zamčený úplně.

    Zatím pokrývá jen solute skupiny (D/R/P). Voda/ionty se řadí ručně přímo
    v data/force_fields.json - přeřazení "profilu" je kurátorské rozhodnutí
    (jaká sada iontů k dané vodě patří), ne prosté přesunutí jména mezi
    seznamy, viz ForceFieldClassificationService.set_solute_tier docstring.
    """
    if not settings.ADMIN_TOKEN or x_admin_token != settings.ADMIN_TOKEN:
        raise ForbiddenError()

    group = payload.get("group")
    ff_name = payload.get("ff_name")
    tier = payload.get("tier")
    if not group or not ff_name or not tier:
        raise BadRequestError("Request body must contain 'group', 'ff_name' and 'tier'.")

    try:
        classification_service.set_solute_tier(group, ff_name, tier)
    except ValueError as e:
        raise BadRequestError(str(e))

    return {"status": "ok", "group": group, "ff_name": ff_name, "tier": tier}
