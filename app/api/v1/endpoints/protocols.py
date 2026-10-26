"""
Protokol celé přípravy (app/services/run_protocol.py): kontrola souboru, který
uživatel načítá, a výchozí protokol pro nové uložení.
"""

from typing import Any, Dict, List

from fastapi import APIRouter, Body
from pydantic import BaseModel

from app.services.run_protocol import RunProtocol, parse_protocol, protocol_filename, protocol_fingerprint

router = APIRouter()


class CheckedProtocol(BaseModel):
    protocol: RunProtocol
    warnings: List[str]
    fingerprint: str
    filename: str


@router.post("/check", response_model=CheckedProtocol, summary="Zkontroluje a doplní protokol přípravy")
async def check_protocol(raw: Dict[str, Any] = Body(...)):
    """
    Načtený soubor (i starší verze formátu) -> úplný protokol s výchozími
    hodnotami, varování (neznámé položky) a otisk pro report. Chyba = 400
    s kódem protocol_invalid / protocol_too_new.
    """
    protocol, warnings = parse_protocol(raw)
    return CheckedProtocol(protocol=protocol, warnings=warnings,
                           fingerprint=protocol_fingerprint(protocol), filename=protocol_filename(protocol))
