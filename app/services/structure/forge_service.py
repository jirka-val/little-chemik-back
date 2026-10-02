"""
Most bere hydrataci/stavbu chybějících atomů/solvataci/ionty výhradně přes
vendorovaný FORGE builder (app/builder), místo starého PDBFixer/OpenMM
HydrogenationService. Builder je striktní vůči svému vstupu (viz
app/builder/INTEGRATION_CONTRACT.md) - očekává už upstream vyčištěnou
strukturu (jeden model, vyřešené AltLocs, správně rozpoznané gap/terminální
varianty) a sám žádnou opravu identity reziduí ani přemostění chybějícího
úseku řetězce neprovádí.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

from app.core.logging import console_logger

from forge_workflow import (  # noqa: E402
    WorkflowResources,
    WorkflowSettings,
    WorkflowResult,
    run_forge_workflow,
)
from forge_molecule_ions import IonPlacementSettings, load_salt_specifications  # noqa: E402
from forge_molecule_parser import (  # noqa: E402
    Molecule,
    build_molecule_from_forge_json,
    load_json,
)
from forge_molecule_solvation import SolvationVdwParameters, SolvationSettings, load_solvation_template  # noqa: E402
from forge_molecule_state_assignment import assign_molecule_states  # noqa: E402

from app.core.config import settings
from app.core.exceptions import AppBaseException
from app.services.structure.groups import BUILDER_GROUPS, WATER_GROUPS
from app.services.structure.pdb_writer import build_forge_meta, molecule_to_pdb
from app.services.structure.reports import REVIEW_SECTIONS, build_histidine_review, protonation_override_map
from app.services.analysis import build_sequence_tokens, required_ff_groups
from app.services.forcefield_service import ForceFieldService
from app.services.structure.structure_review import apply_structure_edits, find_amide_flips, find_zero_occupancy

_DATA_DIR = settings.BASE_DIR / "data"

# Rezidua, která builder umí sám rozpoznat a klasifikovat (polymer/voda/iont) -
# viz _strip_unrecognized_heterogens níže. Ionty jsou tu záměrně, na rozdíl od
# dřívější PDBFixer cesty (HydrogenationService, odstraněna) - builder existující krystalové ionty
# umí sám vyhodnotit a případně nahradit (result.crystal_ion_cleanup), takže je
# netřeba (a nechceme je) stripovat spolu s ligandy.
_KNOWN_POLYMER_RESNAMES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
})
_WATER_RESNAMES = frozenset({"HOH", "WAT", "SOL"})
_KNOWN_ION_RESNAMES = frozenset({
    "NA", "CL", "K", "MG", "CA", "LI", "RB", "CS", "ZN", "F", "BR", "I",
})


def _resolve_mol_type(key: str, ff_data: Dict[str, Any]) -> str:
    """
    ff_selections je klíčované tím, co detekuje pdb_service.get_molecule_types
    ("W" pro vodu, "I" pro ionty obecně) - to je ale legacy vokabulář z doby
    před FORGE builderem. Builder potřebuje přesný podtyp (W3/W4/W5 podle
    počtu bodů modelu vody, I1/I1+/Im/Im+ podle mocenství iontu), protože na
    něm závisí i to, pod jakým klíčem se zaregistrují LJ parametry pro
    solvataci (viz SolvationVdwParameters.from_force_field_directories -
    mol_type se odvozuje z názvu adresáře "{ff_name}_{mol_type}").
    Použití obecného klíče přímo by LJ parametry zaregistrovalo pod "W"
    místo "W3" a solvatace by pak spadla na KeyError ("LJ sigma missing for
    W3:WAT:O") - přesně tenhle bug řeší tahle funkce.
    """
    if key in BUILDER_GROUPS:
        return key

    candidates = ff_data.get("molecule_type") or []
    precise = [c for c in candidates if c in BUILDER_GROUPS]

    if len(precise) == 1:
        return precise[0]
    if len(precise) > 1:
        raise ValueError(
            f"Ambiguous FORGE mol_type for ff_selections key {key!r}: "
            f"force field declares multiple candidate groups {precise!r} - "
            "cannot infer a single one automatically."
        )
    raise ValueError(
        f"Cannot resolve a FORGE mol_type for ff_selections key {key!r} "
        f"(force field's own molecule_type={candidates!r} contains no "
        "recognized builder group)."
    )


def _strip_unrecognized_heterogens(pdb_text: str, crystal_water_mode: str) -> str:
    """
    Zrcadlí chování dřívějšího PDBFixer.removeHeterogens() (HydrogenationService),
    které se při migraci na FORGE builder ztratilo. Builder neumí stavět
    libovolné krystalizační ligandy/aditiva (GOL, SO4, EDO, ...) - jen je
    tiše propustí do výstupu (Molecule.passthrough_atoms) beze změny, a
    TopologyService na nich pak spadne s KeyError, protože pro ně neexistuje
    žádná FF šablona (potvrzeno reálným pádem na 3DVZ: KeyError('GOL')).
    """
    if crystal_water_mode == "keep_all":
        return pdb_text

    keep_water = crystal_water_mode == "keep_water"
    out_lines = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in _WATER_RESNAMES:
                if not keep_water:
                    continue
            elif resname not in _KNOWN_ION_RESNAMES and resname not in _KNOWN_POLYMER_RESNAMES:
                continue
        out_lines.append(line)
    return "\n".join(out_lines)


def find_removed_heterogens(pdb_text: str, crystal_water_mode: str) -> List[Dict[str, Any]]:
    """
    Ligandy/aditiva, které _strip_unrecognized_heterogens zahodí (builder pro
    ně nemá šablonu) - pro krok 3.5, ať to uživatel ví. Voda sem nepatří,
    tu řídí crystal_water_mode.
    """
    if crystal_water_mode == "keep_all":
        return []
    found: Dict[tuple, Dict[str, Any]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in _WATER_RESNAMES or resname in _KNOWN_ION_RESNAMES or resname in _KNOWN_POLYMER_RESNAMES:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        chain, icode = line[21].strip() or "?", line[26].strip()
        entry = found.setdefault((chain, resseq, icode, resname), {
            "key": f"{chain}:{resseq}:{icode}:{resname}",
            "residue": f"{chain}{resseq}{icode}",
            "resname": resname,
            "atoms": 0,
        })
        entry["atoms"] += 1
    return list(found.values())



class ForgeMissingDOFError(AppBaseException):
    """
    Builder narazil na chybějící stupeň volnosti (typicky nekompletní postranní
    řetězec/báze bez jednoznačně určitelné geometrie) a zastavil se, než cokoliv
    dostavil natvrdo. Tohle je legitimní modelovaná hranice, ne chyba parsování -
    viz INTEGRATION_CONTRACT.md. Vyžaduje rozhodnutí uživatele, ne tichou opravu.
    """

    status_code = 409
    code = "missing_dof"

    def __init__(self, step: Any, molecule: Molecule):
        reason_atom = step.reason_atom
        residue = molecule.chains[reason_atom.chain_id].residues[reason_atom.residue_index]
        super().__init__(
            "Structure requires a manual decision the builder cannot make automatically.",
            payload={
                "chain": residue.chain_id,
                "resseq": residue.resseq,
                "icode": residue.icode,
                "ff_resname": residue.ff_resname,
                "original_resname": residue.original_resname,
                "atom_name": reason_atom.atom_name,
            },
        )


class ForgeMissingForceFieldError(AppBaseException):
    """
    ff_selections nepokrývá všechny FORGE mol_type skupiny, které tahle
    struktura reálně potřebuje. Dva zdroje:

    1. Statická před-kontrola (viz prepare_structure) přes
       analysis.ff_requirements.required_ff_groups - odchytí to DŘÍV, než se vůbec
       spustí drahý (u velkých struktur i několikaminutový) builder run.
    2. Bezpečnostní síť kolem run_forge_workflow() - pokud se přesto
       během buildu/solvatace/iontů narazí na chybějící MM/LJ/iontové
       parametry (KeyError z app/builder), překlopí se sem místo syrové 500.

    V obou případech jde o stejnou upstream chybu (chybí nebo je špatně
    vybrané FF pro konkrétní mol_type skupinu, typicky ionty - "Im" pro
    Mg2+ je matoucí název, snadno se zamění za "I1+"), ne o pád kódu -
    proto 409 (konzistentní s ForgeMissingDOFError), ne 500.
    """

    status_code = 409
    code = "missing_force_field"

    def __init__(self, missing: Dict[str, Any], detail: Optional[str] = None):
        if missing:
            groups = ", ".join(sorted(missing.keys()))
            message = f"ff_selections is missing coverage for required mol_type group(s): {groups}."
        else:
            message = "The builder could not find MM/LJ parameters for an atom or ion in the selected force fields."
        if detail:
            message += f" ({detail})"
        super().__init__(message, payload={"missing_groups": missing, "detail": detail})


@dataclass(frozen=True)
class _StaticResources:
    converting_dictionary: Any
    building_template: Any
    state_definitions: Any
    water_template: Any


@lru_cache(maxsize=1)
def _static_resources() -> _StaticResources:
    return _StaticResources(
        converting_dictionary=load_json(_DATA_DIR / "converting_dictionary.json"),
        building_template=load_json(_DATA_DIR / "building_template_v1.json"),
        state_definitions=load_json(_DATA_DIR / "protonation_states_v1.json"),
        water_template=load_solvation_template(_DATA_DIR / "solvation_template_v1.json"),
    )


@lru_cache(maxsize=32)
def _cached_ff_parameters(directories: tuple) -> SolvationVdwParameters:
    return SolvationVdwParameters.from_force_field_directories([Path(d) for d in directories])



@dataclass
class ForgeWorkflowRun:
    """
    Výsledek run_forge_workflow() spolu se zdroji/nastavením, kterými byl
    spuštěn - sidechain_service.py je potřebuje znovu (mm_parameters pro MM
    optimalizaci side-chainů, settings/salts pro navazující solvataci/ionty
    po přijetí GUI voleb).
    """

    result: Optional[WorkflowResult]
    resources: WorkflowResources
    settings: WorkflowSettings
    salts: List[Any]
    # Vyplněné (a result=None), když se run_workflow zastavil před stavbou
    # kvůli kroku 3.5 "Structure Check" v Expert režimu - viz _structure_review.
    structure_review: Optional[Dict[str, Any]] = None


@dataclass
class StructureDecisions:
    """
    Rozhodnutí uživatele z kroku 3.5 (mimo HIS, ty jdou přes
    protonation_overrides). *_decided = o čem už rozhodl (i "ponechat"),
    aby se na totéž při další přípravě znovu neptal.
    """

    flips: List[tuple]
    amide_flips_decided: List[tuple]
    rebuilds: List[tuple]
    zero_occupancy_decided: List[tuple]
    acknowledged_heterogens: set

    @classmethod
    def from_payload(cls, payload: Optional[Dict[str, Any]]) -> "StructureDecisions":
        payload = payload or {}

        def key(item: Dict[str, Any]) -> tuple:
            return (item["chain"], int(item["resseq"]), item.get("icode") or "")

        flips = payload.get("amide_flips") or []
        zero = payload.get("zero_occupancy") or []
        return cls(
            flips=[key(i) for i in flips if i.get("apply")],
            amide_flips_decided=[key(i) for i in flips],
            rebuilds=[key(i) for i in zero if i.get("apply")],
            zero_occupancy_decided=[key(i) for i in zero],
            acknowledged_heterogens=set(payload.get("acknowledged_heterogens") or []),
        )


@dataclass
class ForgePreparationResult:
    pdb_text: str
    forge_meta: Dict[str, Dict[str, Any]]
    warnings: List[str]
    state_assignment: Any
    crystal_ion_cleanup: Any
    solvation: Any
    ion_addition: Any



class ForgeStructureService:
    """Bridges upstream-cleaned PDB structures (app/services/analysis) to app/builder."""

    def __init__(self):
        self.ff_service = ForceFieldService()

    def _check_ff_coverage(
        self,
        pdb_text: str,
        ff_selections: Dict[str, Any],
        add_solvent_and_ions: bool,
        salts: Optional[List[Dict[str, Any]]],
        require_mg: bool = False,
    ) -> None:
        """
        Ověří PŘED spuštěním buildu, že ff_selections pokrývá všechny
        mol_type skupiny, které tahle konkrétní struktura potřebuje (viz
        analysis.ff_requirements.required_ff_groups) - ať se chybějící/špatně
        vybrané FF (typicky ionty, "Im" pro Mg2+ vs "I1+") odhalí hned,
        ne až po několikaminutovém běhu buildu/solvatace pádem s KeyError.
        """
        covered = set()
        for key, ff_data in ff_selections.items():
            try:
                covered.add(_resolve_mol_type(key, ff_data))
            except ValueError:
                continue

        required = required_ff_groups(
            pdb_text, add_solvent_and_ions=add_solvent_and_ions, salts=salts, require_mg=require_mg
        )
        missing = {}
        for mol_type, info in required.items():
            # "W" je záměrně obecný požadavek (nezáleží, jaký konkrétní
            # vodní model - viz required_ff_groups), zatímco ff_selections
            # se přes _resolve_mol_type vždy rozřeší na přesný podtyp
            # (W3/W4/W5) - porovnávat je proto nutné přes celou skupinu.
            if mol_type == "W":
                if not covered & WATER_GROUPS:
                    missing[mol_type] = info
            elif mol_type not in covered:
                missing[mol_type] = info
        if missing:
            raise ForgeMissingForceFieldError(missing)

    def _resolve_force_field_parameters(self, ff_selections: Dict[str, Any]) -> SolvationVdwParameters:
        if not ff_selections:
            raise ValueError("ff_selections must not be empty - the builder needs a force field per mol_type.")

        directories = []
        for key, ff_data in ff_selections.items():
            resolved_mol_type = _resolve_mol_type(key, ff_data)
            target_dir = self.ff_service.prepare_forge_force_field_directory(ff_data, resolved_mol_type)
            directories.append(str(target_dir))

        # Explicitní seznam adresářů (from_force_field_directories), NE
        # from_force_field_root nad celou sdílenou ff_cache_forge cache - ta je
        # společná napříč workspace i FF volbami a from_force_field_root by při
        # dvou různých FF se stejným mol_type spadl na ValueError (duplicitní
        # mol_type v jednom rootu).
        return _cached_ff_parameters(tuple(sorted(directories)))

    def _build_resources(self, ff_selections: Dict[str, Any]) -> WorkflowResources:
        static = _static_resources()
        return WorkflowResources(
            converting_dictionary=static.converting_dictionary,
            building_template=static.building_template,
            state_definitions=static.state_definitions,
            water_template=static.water_template,
            force_field_parameters=self._resolve_force_field_parameters(ff_selections),
        )

    @staticmethod
    def _structure_review(
        original_pdb: str,
        crystal_water_mode: str,
        decisions: "StructureDecisions",
        structure_data: Dict[str, Any],
        resources: WorkflowResources,
        settings: WorkflowSettings,
    ) -> Optional[Dict[str, Any]]:
        """
        Krok 3.5 "Structure Check": vrátí sekce s nevyřízenými položkami,
        nebo None, když není o čem rozhodovat. Nulová obsazenost, amidy a
        heterogeny se hledají v původním PDB (bez položek, o kterých už
        uživatel rozhodl); HIS na struktuře s už aplikovanými rozhodnutími.
        """
        review = {
            **ForgeStructureService._histidine_review(structure_data, resources, settings),
            "zero_occupancy": find_zero_occupancy(original_pdb, decided=decisions.zero_occupancy_decided),
            "amide_flips": find_amide_flips(original_pdb, decided=decisions.amide_flips_decided),
            "removed_heterogens": [
                h for h in find_removed_heterogens(original_pdb, crystal_water_mode)
                if h["key"] not in decisions.acknowledged_heterogens
            ],
        }
        if not any(review[section] for section in REVIEW_SECTIONS):
            return None
        return review

    @staticmethod
    def _histidine_review(
        structure_data: Dict[str, Any], resources: WorkflowResources, settings: WorkflowSettings
    ) -> Dict[str, Any]:
        """
        Samostatné přiřazení stavů (stejné volání jako na začátku
        run_forge_workflow) jen kvůli reportu - bez stavby atomů, takže je
        levné. Builder tím zůstává beze změny: zastavení je věc služby.
        """
        molecule = build_molecule_from_forge_json(structure_data, resources.converting_dictionary)
        _molecule, report = assign_molecule_states(
            molecule,
            resources.converting_dictionary,
            resources.building_template,
            resources.state_definitions,
            pH=settings.pH,
            covalent_cutoff_angstrom=settings.covalent_cutoff_angstrom,
            hydrogen_bond_settings=settings.hydrogen_bond,
            forced_protonation_states=settings.protonation_overrides,
            modify_myself=True,
        )
        return build_histidine_review(report.protonation)

    def run_workflow(
        self,
        pdb_text: str,
        ff_selections: Dict[str, Any],
        ph: float = 7.0,
        add_solvent_and_ions: bool = True,
        salts: Optional[List[Dict[str, Any]]] = None,
        box_shape: Optional[str] = None,
        box_padding_angstrom: Optional[float] = None,
        keep_crystal_waters: Optional[bool] = None,
        crystal_water_mode: str = "remove_all",
        clean_crystal_ions: Optional[bool] = None,
        replace_structural_multivalent_with_mg: Optional[bool] = None,
        concentration_mode: Optional[str] = None,
        protonation_overrides: Optional[List[Dict[str, Any]]] = None,
        structure_decisions: Optional[Dict[str, Any]] = None,
        review_structure: bool = False,
        auto_amide_flips: bool = False,
    ) -> "ForgeWorkflowRun":
        """
        Sdílené jádro mezi neinteraktivním `prepare_structure()` a interaktivním
        side-chain flow (viz sidechain_service.py) - ff-coverage kontrola, sestavení
        WorkflowResources/WorkflowSettings a samotné spuštění run_forge_workflow().
        Rozhodnutí, co dělat s `result.stopped_at_missing_dof` (409 vs. otevření
        interaktivní session), zůstává na volajícím.

        review_structure=True (Expert, krok 3.5 "Structure Check"): když je o
        čem rozhodovat (nejistý HIS, nulová obsazenost, otočený amid,
        odstraněný ligand), workflow se vůbec nespustí a vrátí se
        `structure_review` s result=None. Frontend po kontrole pošle přípravu
        znovu s rozhodnutími (protonation_overrides, structure_decisions) a
        review_structure=False.

        auto_amide_flips=True (Guided/Standard): doporučená otočení amidů
        ASN/GLN (find_amide_flips) se použijí bez ptaní.
        """
        logger.info(
            f"FORGE: Preparing structure (pH={ph}, add_solvent_and_ions={add_solvent_and_ions}, "
            f"crystal_water_mode={crystal_water_mode}, ff_selections={list(ff_selections.keys())})"
        )
        console_logger.info("Preparing structure...")

        decisions = StructureDecisions.from_payload(structure_decisions)
        original_pdb = pdb_text
        flips = list(decisions.flips)
        if auto_amide_flips:
            auto = find_amide_flips(original_pdb, decided=decisions.amide_flips_decided)
            flips += [(f["chain"], f["resseq"], f["icode"]) for f in auto]
            if auto:
                names = ", ".join(f"{f['resname']} {f['residue']}" for f in auto)
                logger.info(f"FORGE: Auto-flipped amides: {names}")
                console_logger.info(f"Flipped {len(auto)} ASN/GLN amide(s) to fit hydrogen bonds: {names}.")
        pdb_text = apply_structure_edits(
            pdb_text, flip_residues=flips, rebuild_residues=decisions.rebuilds
        )
        pdb_text = _strip_unrecognized_heterogens(pdb_text, crystal_water_mode)
        sequence_data = build_sequence_tokens(pdb_text, chain=None, fill_gaps=True)
        structure_data = {"pdb_text": pdb_text, "missing_atoms": sequence_data}

        self._check_ff_coverage(
            pdb_text, ff_selections, add_solvent_and_ions, salts,
            require_mg=bool(replace_structural_multivalent_with_mg),
        )

        resources = self._build_resources(ff_selections)
        salt_specs = load_salt_specifications({"salts": salts or []})

        solvation_kwargs = {}
        if box_shape is not None:
            solvation_kwargs["box_shape"] = box_shape
        if box_padding_angstrom is not None:
            solvation_kwargs["padding_angstrom"] = box_padding_angstrom
        if keep_crystal_waters is not None:
            solvation_kwargs["keep_crystal_waters"] = keep_crystal_waters

        ion_kwargs = {}
        if clean_crystal_ions is not None:
            ion_kwargs["clean_crystal_ions"] = clean_crystal_ions
        if replace_structural_multivalent_with_mg is not None:
            ion_kwargs["replace_structural_multivalent_with_mg"] = replace_structural_multivalent_with_mg
        if concentration_mode is not None:
            ion_kwargs["concentration_mode"] = concentration_mode

        settings = WorkflowSettings(
            pH=ph,
            add_solvent_and_ions=add_solvent_and_ions,
            solvation=SolvationSettings(**solvation_kwargs),
            ions=IonPlacementSettings(**ion_kwargs),
            protonation_overrides=protonation_override_map(protonation_overrides),
        )

        if review_structure:
            review = self._structure_review(
                original_pdb, crystal_water_mode, decisions, structure_data, resources, settings
            )
            if review is not None:
                counts = ", ".join(f"{s}={len(review[s])}" for s in REVIEW_SECTIONS if review[s])
                logger.info(f"FORGE: Stopped for structure check - {counts}.")
                console_logger.warning("Preparation paused - the structure needs a manual check.")
                return ForgeWorkflowRun(
                    result=None, resources=resources, settings=settings, salts=salt_specs,
                    structure_review=review,
                )

        try:
            result: WorkflowResult = run_forge_workflow(
                structure_data,
                resources,
                salts=salt_specs,
                settings=settings,
            )
        except KeyError as exc:
            # Bezpečnostní síť pro chybějící MM/LJ/iontové parametry, které
            # _check_ff_coverage výše z nějakého důvodu neodchytila (např.
            # konkrétní rezidum/atom chybí ve vybraném FF, i když formálně
            # mol_type skupina pokrytá je). app/builder tyhle KeyError hlásí
            # ve třech rozlišitelných formátech - viz forge_molecule_ions.py
            # _ion_params a forge_molecule_solvation.py/forge_molecule_builder.py
            # atom_params. Cokoliv jiného (skutečný programátorský bug)
            # necháváme propadnout jako dřív.
            detail = exc.args[0] if exc.args else str(exc)
            if isinstance(detail, str) and (
                "parameters missing" in detail or "LJ sigma missing" in detail
            ):
                console_logger.error("Structure preparation failed - missing force field parameters.")
                raise ForgeMissingForceFieldError({}, detail=detail) from exc
            console_logger.error("Structure preparation failed.")
            raise

        if result.stopped_at_missing_dof:
            logger.info(
                f"FORGE: Build stopped at a missing degree of freedom - "
                f"{len(result.remaining_plan.steps)} step(s) remaining, awaiting user decision."
            )
            console_logger.warning("Build paused - a residue needs a manual decision.")
        else:
            logger.info("FORGE: Build complete (all residues resolved).")
            if result.crystal_ion_cleanup:
                c = result.crystal_ion_cleanup
                logger.info(
                    f"FORGE: Crystal ion cleanup - input={c.input_ions}, "
                    f"removed_monovalent={c.removed_monovalent}, "
                    f"removed_nonstructural_multivalent={c.removed_nonstructural_multivalent}, "
                    f"retained={c.retained_structural_monovalent + c.retained_structural_multivalent}"
                )
            if result.solvation:
                s = result.solvation
                logger.info(
                    f"FORGE: Solvation complete - box={s.box_shape} (padding={s.padding_angstrom}A), "
                    f"crystal_waters_retained={s.retained_crystal_waters}, "
                    f"generated_waters={s.generated_waters}, total_waters={s.total_waters}"
                )
            if result.ion_addition:
                i = result.ion_addition
                logger.info(
                    f"FORGE: Ion placement complete - neutralization={i.neutralization_ions}, "
                    f"added={i.added_ions}, final_system_charge={i.final_system_charge}"
                )

        return ForgeWorkflowRun(result=result, resources=resources, settings=settings, salts=salt_specs)

    def prepare_structure(
        self,
        pdb_text: str,
        ff_selections: Dict[str, Any],
        ph: float = 7.0,
        add_solvent_and_ions: bool = True,
        salts: Optional[List[Dict[str, Any]]] = None,
        box_shape: Optional[str] = None,
        box_padding_angstrom: Optional[float] = None,
        keep_crystal_waters: Optional[bool] = None,
        crystal_water_mode: str = "remove_all",
        clean_crystal_ions: Optional[bool] = None,
        replace_structural_multivalent_with_mg: Optional[bool] = None,
        concentration_mode: Optional[str] = None,
        protonation_overrides: Optional[List[Dict[str, Any]]] = None,
        structure_decisions: Optional[Dict[str, Any]] = None,
        auto_amide_flips: bool = False,
    ) -> ForgePreparationResult:
        """
        Spustí kompletní FORGE zpracování (stavy/protonace, stavba chybějících
        atomů, solvatace, ionty) na už upstream vyčištěné struktuře (jeden model,
        vyřešené AltLocs, aplikovaná symetrie - viz analysis.structure_prep.process_structure).

        Vyhodí ForgeMissingDOFError, pokud builder narazí na chybějící stupeň
        volnosti, který nejde bezpečně dostavět - to volající musí propustit
        uživateli, ne potichu obejít. Pro interaktivní dostavění bezpečných
        side-chain větví viz sidechain_service.SidechainSessionService, který
        používá run_workflow() přímo místo tohohle wrapperu.
        """
        run = self.run_workflow(
            pdb_text,
            ff_selections,
            ph=ph,
            add_solvent_and_ions=add_solvent_and_ions,
            salts=salts,
            box_shape=box_shape,
            box_padding_angstrom=box_padding_angstrom,
            keep_crystal_waters=keep_crystal_waters,
            crystal_water_mode=crystal_water_mode,
            clean_crystal_ions=clean_crystal_ions,
            replace_structural_multivalent_with_mg=replace_structural_multivalent_with_mg,
            concentration_mode=concentration_mode,
            protonation_overrides=protonation_overrides,
            structure_decisions=structure_decisions,
            auto_amide_flips=auto_amide_flips,
        )
        result = run.result

        if result.stopped_at_missing_dof:
            raise ForgeMissingDOFError(result.remaining_plan.steps[0], result.molecule)

        for warning in result.molecule.warnings:
            logger.warning(f"FORGE: {warning}")
        logger.info("FORGE: Structure preparation finished successfully.")

        if result.molecule.warnings:
            console_logger.info(f"Structure prepared successfully ({len(result.molecule.warnings)} warning(s)).")
        else:
            console_logger.info("Structure prepared successfully.")

        return ForgePreparationResult(
            pdb_text=molecule_to_pdb(result.molecule),
            forge_meta=build_forge_meta(result.molecule),
            warnings=list(result.molecule.warnings),
            state_assignment=result.state_assignment,
            crystal_ion_cleanup=result.crystal_ion_cleanup,
            solvation=result.solvation,
            ion_addition=result.ion_addition,
        )
