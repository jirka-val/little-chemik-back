from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, Any, Optional
from app.services.validation.service import ValidationService

router = APIRouter()
validation_service = ValidationService()

# --- Pydantic Modely pro validaci vstupů ---

class ValidationRequest(BaseModel):
    pdb_content: str = Field(..., description="Obsah PDB souboru v textové podobě")
    label: str = Field("molecule_from_front", description="Identifikátor stavu molekuly")

class FixAltLocRequest(BaseModel):
    pdb_content: str = Field(..., description="Původní obsah PDB souboru")
    selections: Dict[str, str] = Field(
        ...,
        description="Mapa výběru variant, např: {'A-42': 'B'}",
        example={"A-42": "B", "A-15": "A"}
    )

# --- Endpointy ---

@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Provede úvodní analýzu molekuly. Detekuje alternativní lokace,
    chybějící atomy a kompatibilitu s forcefieldem.
    """
    try:
        return validation_service.validate_pdb_content(request.pdb_content, request.label)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Chyba při validaci: {str(e)}")

@router.post("/apply-selections", summary="Aplikuje výběr konformací a ověří kontinuitu")
async def apply_selections(request: FixAltLocRequest):
    """
    Tento endpoint je klíčový pro profesionální přípravu:
    1. Vyfiltruje PDB dle výběru (vymaže AltLoc identifikátory).
    2. Okamžitě provede geometrickou kontrolu (C-N vzdálenosti).
    3. Pokud výběr rozbije molekulu, vrátí kritickou chybu v sekci validation.
    """
    try:
        # Voláme novou metodu, která vrací {pdb_content, validation}
        result = validation_service.apply_alt_loc_selection(
            request.pdb_content,
            request.selections
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Chyba při aplikaci selekcí: {str(e)}")


# app/api/v1/endpoints/validation.py

@router.post("/preview-selection", summary="Náhled geometrie před aplikací")
async def preview_selection(request: FixAltLocRequest):
    try:
        # 1. Získáme seznam problémů s kontinuitou (vzdálenosti vazeb mezi rezidui)
        issues = validation_service.conf_manager.validate_continuity(
            request.pdb_content,
            request.selections
        )

        # 2. Vyhodnotíme, zda je výběr bezpečný
        # Za kritické považujeme cokoli, co validate_continuity vrátí
        is_safe = len(issues) == 0

        return {
            "is_ok": is_safe,
            "issues": issues,
            "message": "Výběr je geometricky v pořádku" if is_safe else "Detekovány kritické mezery v řetězci"
        }
    except Exception as e:
        print(f"DEBUG: Preview error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))