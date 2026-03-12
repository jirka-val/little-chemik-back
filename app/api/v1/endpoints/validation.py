from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, Any, Optional
from app.services.validation.service import ValidationService

# Importujeme manažera pracovních prostorů
from app.workspaces.manager import workspace_manager

router = APIRouter()
validation_service = ValidationService()


# --- Pydantic Modely pro validaci vstupů (UPRAVENO PRO NOVÝ FRONTEND) ---

class ValidationRequest(BaseModel):
    workspace_id: str = Field(..., description="ID pracovního prostoru s molekulou")
    label: str = Field("molecule_from_front", description="Identifikátor stavu molekuly")


class FixAltLocRequest(BaseModel):
    workspace_id: str = Field(..., description="ID pracovního prostoru s molekulou")
    selections: Dict[str, str] = Field(
        ...,
        description="Mapa výběru variant, např: {'A-42': 'B'}",
        example={"A-42": "B", "A-15": "A"}
    )


# --- Endpointy ---

@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    """
    Přečte soubor z disku podle workspace_id a provede úvodní analýzu molekuly.
    Detekuje alternativní lokace, chybějící atomy a kompatibilitu s forcefieldem.
    """
    try:
        # 1. Zkontrolujeme, zda soubor existuje pomocí existující metody
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Soubor molekuly nebyl na serveru nalezen.")

        # 2. Získáme cestu k souboru a přečteme obsah
        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        with open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = f.read()

        # 3. Zvalidujeme obsah
        return validation_service.validate_pdb_content(pdb_content, request.label)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Chyba při validaci: {str(e)}")


@router.post("/preview-selection", summary="Náhled geometrie před aplikací")
async def preview_selection(request: FixAltLocRequest):
    """
    Přečte soubor z disku a provede náhled geometrických úprav.
    """
    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Soubor molekuly nebyl na serveru nalezen.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        with open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = f.read()

        # Získáme seznam problémů s kontinuitou
        issues = validation_service.conf_manager.validate_continuity(
            pdb_content,
            request.selections
        )

        is_safe = len(issues) == 0

        return {
            "is_ok": is_safe,
            "issues": issues,
            "message": "Výběr je geometricky v pořádku" if is_safe else "Detekovány kritické mezery v řetězci"
        }
    except HTTPException:
        raise
    except Exception as e:
        print(f"DEBUG: Preview error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/apply-selections", summary="Aplikuje výběr konformací a ověří kontinuitu")
async def apply_selections(request: FixAltLocRequest):
    """
    Aplikuje výběr uživatele, ověří geometrii a ULOŽÍ NOVÝ SOUBOR na serveru.
    """
    try:
        if not workspace_manager.workspace_exists(request.workspace_id):
            raise HTTPException(status_code=404, detail="Soubor molekuly nebyl na serveru nalezen.")

        pdb_path = workspace_manager.get_file_path(request.workspace_id)
        with open(pdb_path, "r", encoding="utf-8") as f:
            pdb_content = f.read()

        # Aplikujeme výběr (vrátí dict s novým 'pdb_content' a 'validation')
        result = validation_service.apply_alt_loc_selection(
            pdb_content,
            request.selections
        )

        # Přepíšeme starý soubor na disku nově upraveným obsahem
        if "pdb_content" in result:
            with open(pdb_path, "w", encoding="utf-8") as f:
                f.write(result["pdb_content"])

        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Chyba při aplikaci selekcí: {str(e)}")