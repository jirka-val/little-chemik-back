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
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

_BUILDER_DIR = Path(__file__).resolve().parents[2] / "builder"
if str(_BUILDER_DIR) not in sys.path:
    sys.path.insert(0, str(_BUILDER_DIR))

from forge_workflow import (  # noqa: E402
    WorkflowResources,
    WorkflowSettings,
    WorkflowResult,
    run_forge_workflow,
)
from forge_molecule_ions import load_salt_specifications  # noqa: E402
from forge_molecule_parser import (  # noqa: E402
    Molecule,
    Residue,
    format_pdb_atom_line,
    infer_element,
    load_json,
)
from forge_molecule_solvation import SolvationVdwParameters, SolvationSettings, load_solvation_template  # noqa: E402

from app.core.exceptions import AppBaseException
from app.services.analysis_service import build_sequence_tokens, required_ff_groups
from app.services.forcefield_service import ForceFieldService

_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_ION_GROUPS = frozenset({"I1", "I1+", "Im", "Im+"})
_WATER_GROUPS = frozenset({"W3", "W4", "W5"})

# Rezidua, která builder umí sám rozpoznat a klasifikovat (polymer/voda/iont) -
# viz _strip_unrecognized_heterogens níže. Ionty jsou tu záměrně, na rozdíl od
# staré PDBFixer cesty (hydrogenation.py) - builder existující krystalové ionty
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

# Přesné skupiny, jak je builder sám rozlišuje (viz converting_dictionary.json
# top-level klíče). ff_selections od frontendu je dnes klíčované obecným
# detekovaným typem ("W", "I"), ne tímhle - viz _resolve_mol_type níže.
_ALL_BUILDER_GROUPS = frozenset({"R", "D", "P"}) | _WATER_GROUPS | _ION_GROUPS


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
    if key in _ALL_BUILDER_GROUPS:
        return key

    candidates = ff_data.get("molecule_type") or []
    precise = [c for c in candidates if c in _ALL_BUILDER_GROUPS]

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
    Zrcadlí chování staré PDBFixer.removeHeterogens() (viz hydrogenation.py),
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


class ForgeWriterError(RuntimeError):
    """Vyhozeno, když by výsledná molekula nešla bezpečně zapsat do fixed-column PDB."""


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
       analysis_service.required_ff_groups - odchytí to DŘÍV, než se vůbec
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


def _pdb_safe_resname(residue: Residue) -> str:
    """
    Vrátí resname zapsatelné do 3sloupcového PDB pole. RNA/DNA varianty
    (RU3/RA5/...) do 3 znaků vždy vejdou. Proteinové terminální varianty
    (CGLU/NPHE/...) jsou 4 znaky - format_pdb_atom_line je NEOŘEZÁVÁ, jen by
    tiše posunul všechny další sloupce na řádku, takže se vždy vrací
    original_resname (builder ho u téhle mutace nepřepisuje).
    """
    if len(residue.ff_resname) <= 3:
        return residue.ff_resname
    original = residue.original_resname
    if original and len(original) <= 3:
        return original
    raise ForgeWriterError(
        f"Residue {residue.chain_id}:{residue.resseq}{residue.icode} "
        f"({residue.ff_resname!r}) has no PDB-safe (<=3 char) representation."
    )


def _format_ter_line(serial: int, resname: str, chain_id: str, resseq: int, icode: str) -> str:
    """
    Standardní PDB TER záznam - signalizuje downstream nástrojům (Mol* mimo
    jiné), že tady polymerní řetězec končí. Bez něj hrozí, že se poslední
    reziduum jednoho chainu a první reziduum dalšího vyhodnotí jako přerušený
    (gap) polymer téhož řetězce, což se ve vieweru projeví tečkovanou
    "vazbou" mezi dvěma chainy i po nastavení jejich terminality - viz
    molecule_to_pdb().
    """
    return f"TER   {serial:5d}      {resname:>3s} {chain_id[:1]:1s}{resseq:4d}{icode[:1]:1s}"


def _cryst1_line(molecule: Molecule) -> Optional[str]:
    import math

    if molecule.periodic_box is None:
        return None
    vectors = molecule.periodic_box.vectors
    lengths = [math.sqrt(sum(x * x for x in vector)) for vector in vectors]

    def angle(left: int, right: int) -> float:
        dot = sum(vectors[left][i] * vectors[right][i] for i in range(3))
        cosine = max(-1.0, min(1.0, dot / (lengths[left] * lengths[right])))
        return math.degrees(math.acos(cosine))

    alpha, beta, gamma = angle(1, 2), angle(0, 2), angle(0, 1)
    return (
        f"CRYST1{lengths[0]:9.3f}{lengths[1]:9.3f}{lengths[2]:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
    )


def _chain_sort_key(molecule: Molecule, chain_id: str):
    residues = molecule.chains[chain_id].residues
    is_water = bool(residues) and all(r.group in _WATER_GROUPS for r in residues)
    is_ion = bool(residues) and all(r.group in _ION_GROUPS for r in residues)
    return (2 if is_water else 1 if is_ion else 0, chain_id)


def _pdb_serial(serial: int) -> int:
    """
    Sériové číslo atomu ve fixed-column PDB smí mít nejvýš 5 číslic
    (sloupce 7-11). format_pdb_atom_line() to samo nehlídá - `f"{serial:5d}"`
    u čísla >= 100000 tiše přeteče na 6 znaků a posune všechny další sloupce
    na řádku o jeden doprava, takže resname/chain/souřadnice skončí na
    špatné pozici. Potvrzeno pádem na solvatovaném 1JJ2 (925 179 atomů):
    OpenMM (přes PDBFixer ve StructureChecker) na takhle posunutém řádku
    spadne na "Misaligned residue name". Sériové číslo je čistě kosmetický
    popisek (nic downstream ho nepoužívá jako identitu - všude se pracuje
    přes chain/resseq/atom name), takže cyklické zabalení zpátky do rozsahu
    1-99999 je bezpečné a zachová platný fixed-column formát i nad hranicí
    legacy PDB limitu.
    """
    return ((serial - 1) % 99999) + 1


def molecule_to_pdb(molecule: Molecule) -> str:
    """
    Zapíše výsledek FORGE builderu do PDB textu pro Molstar/downstream nástroje.
    Reimplementace vzoru z (nevendorovaného) forge_builder_v0/examples/forge_workflow_cli.py,
    doplněná o bezpečnou volbu resname (viz _pdb_safe_resname) - ta v příkladovém
    CLI writeru chybí a u proteinových terminálních variant/nahrazených iontů by
    tiše poškodila fixed-column formát.
    """
    lines: List[str] = []
    cryst1 = _cryst1_line(molecule)
    if cryst1:
        lines.append(cryst1)

    serial = 1
    for chain_id in sorted(molecule.chains, key=lambda cid: _chain_sort_key(molecule, cid)):
        # _chain_sort_key()[0] je 0 jen pro "normální" (ne čistě voda/ionty)
        # chainy - TER dává smysl jen pro tyhle polymerní řetězce, HETATM
        # voda/ionty žádnou spojitost/gap logiku ve vieweru nespouští.
        chain_is_polymer = _chain_sort_key(molecule, chain_id)[0] == 0
        last_written: Optional[tuple] = None
        for residue in molecule.chains[chain_id].residues:
            is_ion = residue.group in _ION_GROUPS
            hetero = is_ion or residue.group in _WATER_GROUPS
            for atom in residue.atoms.values():
                if atom.coord is None:
                    continue
                if is_ion:
                    # Monatomární ionty: FF resname (Mg2+, Cs+, ...) do 3 sloupců
                    # nevejde. Builder sám tenhle případ řeší přes element symbol
                    # (viz forge_molecule_ions._ion_element) - držíme se stejné
                    # konvence, ať je psaní iontů konzistentní s tím, co builder
                    # sám interně považuje za jejich identitu.
                    element = atom.element or infer_element(atom.name)
                    if not element:
                        raise ForgeWriterError(
                            f"Ion {residue.chain_id}:{residue.resseq} ({residue.ff_resname}) "
                            "has no resolvable element symbol."
                        )
                    atom_name_out = element.upper()
                    resname_out = element.upper()
                else:
                    atom_name_out = atom.name
                    resname_out = _pdb_safe_resname(residue)

                lines.append(
                    format_pdb_atom_line(
                        serial=_pdb_serial(serial),
                        record_name="HETATM" if hetero else "ATOM",
                        atom_name=atom_name_out,
                        resname=resname_out,
                        chain_id=residue.chain_id,
                        resseq=residue.resseq,
                        icode=residue.icode,
                        coord=atom.coord,
                        occupancy=atom.occupancy if atom.occupancy is not None else 1.0,
                        bfactor=atom.bfactor if atom.bfactor is not None else 0.0,
                        element=atom.element,
                        altloc="",
                    )
                )
                serial += 1
                if chain_is_polymer:
                    last_written = (resname_out, residue.resseq, residue.icode)

        if chain_is_polymer and last_written is not None:
            resname_out, resseq, icode = last_written
            lines.append(_format_ter_line(_pdb_serial(serial), resname_out, chain_id, resseq, icode))
            serial += 1

    for record in molecule.passthrough_atoms:
        is_ion = record.group in _ION_GROUPS
        if is_ion:
            element = record.element or infer_element(record.atom_name)
            atom_name_out = (element or record.atom_name).upper()
            resname_out = (element or record.resname[:3]).upper()
        else:
            atom_name_out = record.atom_name
            resname_out = record.resname[-3:] if len(record.resname) > 3 else record.resname

        lines.append(
            format_pdb_atom_line(
                serial=_pdb_serial(serial),
                record_name="HETATM",
                atom_name=atom_name_out,
                resname=resname_out,
                chain_id=record.chain_id,
                resseq=record.resseq,
                icode=record.icode,
                coord=record.coord,
                occupancy=record.occupancy if record.occupancy is not None else 1.0,
                bfactor=record.bfactor if record.bfactor is not None else 0.0,
                element=record.element,
                altloc="",
            )
        )
        serial += 1

    lines.extend(("END", ""))
    return "\n".join(lines)


def build_forge_meta(molecule: Molecule) -> Dict[str, Dict[str, Any]]:
    """
    Sidecar metadata (chain:resseq:icode -> autoritativní ff_resname/group) pro
    TopologyService. Builder do PDB textu zapisuje jen 3znakovou reprezentaci
    (viz _pdb_safe_resname), takže proteinové terminální varianty jako CGLU/NPHE
    by se po zpětném parsování PDB ztratily. Tohle je jediné místo, kde je
    plná (i 4znaková) identita reziduí po doběhnutí state-assignmentu dostupná.
    """
    meta: Dict[str, Dict[str, Any]] = {}
    for chain in molecule.chains.values():
        for residue in chain.residues:
            key = f"{residue.chain_id}:{residue.resseq}:{residue.icode}"
            meta[key] = {"ff_resname": residue.ff_resname, "group": residue.group}
    return meta


@dataclass
class ForgeWorkflowRun:
    """
    Výsledek run_forge_workflow() spolu se zdroji/nastavením, kterými byl
    spuštěn - sidechain_service.py je potřebuje znovu (mm_parameters pro MM
    optimalizaci side-chainů, settings/salts pro navazující solvataci/ionty
    po přijetí GUI voleb).
    """

    result: WorkflowResult
    resources: WorkflowResources
    settings: WorkflowSettings
    salts: List[Any]


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
    """Bridges upstream-cleaned PDB structures (analysis_service) to app/builder."""

    def __init__(self):
        self.ff_service = ForceFieldService()

    def _check_ff_coverage(
        self,
        pdb_text: str,
        ff_selections: Dict[str, Any],
        add_solvent_and_ions: bool,
        salts: Optional[List[Dict[str, Any]]],
    ) -> None:
        """
        Ověří PŘED spuštěním buildu, že ff_selections pokrývá všechny
        mol_type skupiny, které tahle konkrétní struktura potřebuje (viz
        analysis_service.required_ff_groups) - ať se chybějící/špatně
        vybrané FF (typicky ionty, "Im" pro Mg2+ vs "I1+") odhalí hned,
        ne až po několikaminutovém běhu buildu/solvatace pádem s KeyError.
        """
        covered = set()
        for key, ff_data in ff_selections.items():
            try:
                covered.add(_resolve_mol_type(key, ff_data))
            except ValueError:
                continue

        required = required_ff_groups(pdb_text, add_solvent_and_ions=add_solvent_and_ions, salts=salts)
        missing = {}
        for mol_type, info in required.items():
            # "W" je záměrně obecný požadavek (nezáleží, jaký konkrétní
            # vodní model - viz required_ff_groups), zatímco ff_selections
            # se přes _resolve_mol_type vždy rozřeší na přesný podtyp
            # (W3/W4/W5) - porovnávat je proto nutné přes celou skupinu.
            if mol_type == "W":
                if not covered & _WATER_GROUPS:
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
    ) -> "ForgeWorkflowRun":
        """
        Sdílené jádro mezi neinteraktivním `prepare_structure()` a interaktivním
        side-chain flow (viz sidechain_service.py) - ff-coverage kontrola, sestavení
        WorkflowResources/WorkflowSettings a samotné spuštění run_forge_workflow().
        Rozhodnutí, co dělat s `result.stopped_at_missing_dof` (409 vs. otevření
        interaktivní session), zůstává na volajícím.
        """
        pdb_text = _strip_unrecognized_heterogens(pdb_text, crystal_water_mode)
        sequence_data = build_sequence_tokens(pdb_text, chain=None, fill_gaps=True)
        structure_data = {"pdb_text": pdb_text, "missing_atoms": sequence_data}

        self._check_ff_coverage(pdb_text, ff_selections, add_solvent_and_ions, salts)

        resources = self._build_resources(ff_selections)
        salt_specs = load_salt_specifications({"salts": salts or []})

        solvation_kwargs = {}
        if box_shape is not None:
            solvation_kwargs["box_shape"] = box_shape
        if box_padding_angstrom is not None:
            solvation_kwargs["padding_angstrom"] = box_padding_angstrom
        if keep_crystal_waters is not None:
            solvation_kwargs["keep_crystal_waters"] = keep_crystal_waters

        settings = WorkflowSettings(
            pH=ph,
            add_solvent_and_ions=add_solvent_and_ions,
            solvation=SolvationSettings(**solvation_kwargs),
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
                raise ForgeMissingForceFieldError({}, detail=detail) from exc
            raise

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
    ) -> ForgePreparationResult:
        """
        Spustí kompletní FORGE zpracování (stavy/protonace, stavba chybějících
        atomů, solvatace, ionty) na už upstream vyčištěné struktuře (jeden model,
        vyřešené AltLocs, aplikovaná symetrie - viz analysis_service.process_structure).

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
        )
        result = run.result

        if result.stopped_at_missing_dof:
            raise ForgeMissingDOFError(result.remaining_plan.steps[0], result.molecule)

        return ForgePreparationResult(
            pdb_text=molecule_to_pdb(result.molecule),
            forge_meta=build_forge_meta(result.molecule),
            warnings=list(result.molecule.warnings),
            state_assignment=result.state_assignment,
            crystal_ion_cleanup=result.crystal_ion_cleanup,
            solvation=result.solvation,
            ion_addition=result.ion_addition,
        )
