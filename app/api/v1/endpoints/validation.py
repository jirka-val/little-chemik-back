from fastapi import APIRouter
from pydantic import BaseModel
from typing import Dict
from app.services.validation_service import ValidationService

router = APIRouter()
validation_service = ValidationService()

class ValidationRequest(BaseModel):
    pdb_content: str
    label: str = "molecule_from_front"

class FixAltLocRequest(BaseModel):
    pdb_content: str
    selections: Dict[str, str]  # Např: {"A-42-ARG": "B"}

@router.post("/check", summary="Zvaliduje stav molekuly a detekuje AltLocs")
async def check_molecule(request: ValidationRequest):
    return validation_service.validate_pdb_content(request.pdb_content, request.label)

@router.post("/apply-selections", summary="Vygeneruje čisté PDB podle vybraných konformací")
async def apply_selections(request: FixAltLocRequest):
    cleaned_pdb = validation_service.apply_alt_loc_selection(request.pdb_content, request.selections)
    # Po vyčištění doporučujeme znovu provést validaci
    new_validation = validation_service.validate_pdb_content(cleaned_pdb)
    return {
        "pdb_content": cleaned_pdb,
        "validation": new_validation
    }