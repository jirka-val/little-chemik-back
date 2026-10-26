"""
Protokol celé přípravy (soubor *.forge-protocol.json): nastavení všech kroků
stránky od konformací po export, nezávislé na konkrétní struktuře. Uloží se
z aktuálního stavu stránky a použije na jinou strukturu - načtením do panelů,
nebo automatickým během (frontend, src/core/run-protocol.ts).

Tady je formát, kontrola a převod starších verzí formátu. Protonace
konkrétních reziduí jsou volitelné - builder ty, které na strukturu nesedí,
ignoruje a nahlásí (forced_state_not_applicable).
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, List, Literal, Optional, Tuple, Type

from pydantic import BaseModel, Field, ValidationError

from app.api.v1.endpoints.simulation import AmberMdinRequest
from app.api.v1.schemas.preparation import PreparationSettings, ProtonationOverride
from app.core.exceptions import BadRequestError

FORMAT = "forge-protocol"
FORMAT_VERSION = 1
FILE_SUFFIX = ".forge-protocol.json"


class Conformations(BaseModel):
    # Modely, symetrie, AltLocs a identické kopie podle doporučení backendu
    # (jako Guided). Konkrétní výběr patří ke konkrétní struktuře.
    strategy: Literal["recommended"] = "recommended"


class ForceFieldRef(BaseModel):
    # id z katalogu IDA; podle jména se dohledá, když se id změní
    id: Optional[int] = None
    name: str = Field(..., min_length=1, max_length=200)


class Preparation(PreparationSettings):
    # Nahraný už solvatovaný systém: ponechat jeho box, nebo solvatovat znovu.
    existing_solvent: Literal["keep", "rebuild"] = "keep"


class ProtonationChoice(ProtonationOverride):
    # jméno rezidua při uložení - jen pro přehled, builder ho nepotřebuje
    resname: str = Field("", max_length=10)


class ExportChoice(BaseModel):
    pdb: bool = True
    top: bool = True
    crd: bool = True
    mdin: bool = True
    report: bool = True
    procedure: bool = True
    protocol: bool = True    # přibalit protokol k exportu


class RunProtocol(BaseModel):
    format: Literal["forge-protocol"] = FORMAT
    format_version: int = FORMAT_VERSION
    name: str = Field(..., min_length=1, max_length=100)
    notes: str = Field("", max_length=2000)
    created: Optional[str] = Field(None, max_length=40)
    mode: Literal["guided", "standard", "expert"] = "standard"
    conformations: Conformations = Field(default_factory=Conformations)
    # mol_type (R, P, W, I1 ...) -> FF
    forcefields: Dict[str, ForceFieldRef] = Field(default_factory=dict)
    preparation: Preparation = Field(default_factory=Preparation)
    protonation: List[ProtonationChoice] = Field(default_factory=list, max_length=5000)
    simulation: AmberMdinRequest = Field(default_factory=AmberMdinRequest)
    export: ExportChoice = Field(default_factory=ExportChoice)


# Převody starších verzí formátu: {verze: funkce(dict) -> dict o verzi výš}.
_UPGRADES: Dict[int, Any] = {}


def _unknown_keys(raw: Any, model: Type[BaseModel], path: str = "") -> List[str]:
    """Klíče, které formát nezná (překlep, novější verze) - jen jako varování."""
    if not isinstance(raw, dict):
        return []
    out = []
    for key, value in raw.items():
        field = model.model_fields.get(key)
        where = f"{path}{key}"
        if field is None:
            out.append(where)
            continue
        sub = field.annotation
        if isinstance(sub, type) and issubclass(sub, BaseModel):
            out += _unknown_keys(value, sub, where + ".")
    return out


def _error_text(e: ValidationError) -> str:
    first = e.errors()[0]
    where = ".".join(str(part) for part in first["loc"]) or "protocol"
    return f"{where}: {first['msg']}"


def parse_protocol(raw: Any) -> Tuple[RunProtocol, List[str]]:
    """Načte protokol (i starší verzi formátu); vrací ho a varování."""
    if not isinstance(raw, dict) or raw.get("format") != FORMAT:
        raise BadRequestError("This is not a FORGE protocol file.", code="protocol_invalid")
    version = raw.get("format_version")
    if not isinstance(version, int) or version < 1:
        raise BadRequestError("The protocol has no valid format version.", code="protocol_invalid")
    if version > FORMAT_VERSION:
        raise BadRequestError(
            f"The protocol was saved by a newer version of the application (format {version}, "
            f"this one reads up to {FORMAT_VERSION}).",
            code="protocol_too_new",
        )
    data = dict(raw)
    while data["format_version"] < FORMAT_VERSION:
        data = _UPGRADES[data["format_version"]](data)
    warnings = [f"Unknown setting '{key}' was ignored." for key in _unknown_keys(data, RunProtocol)]
    try:
        protocol = RunProtocol.model_validate(data)
    except ValidationError as e:
        raise BadRequestError(f"The protocol is not valid - {_error_text(e)}", code="protocol_invalid")
    return protocol, warnings


def protocol_fingerprint(protocol: RunProtocol) -> str:
    """Krátký otisk obsahu (bez jména a poznámek) - pro report, ať jde dohledat, čím data vznikla."""
    content = protocol.model_dump(mode="json", exclude={"name", "notes", "created"})
    digest = hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest()
    return digest[:12]


def protocol_filename(protocol: RunProtocol) -> str:
    # jen ASCII - jméno jde do ZIPu i do hlavičky Content-Disposition
    safe = "".join(c if (c.isascii() and c.isalnum()) or c in "-_" else "_" for c in protocol.name.strip())[:60].strip("_")
    return f"{safe or 'protocol'}{FILE_SUFFIX}"
