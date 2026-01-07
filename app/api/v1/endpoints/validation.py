from fastapi import APIRouter
from pydantic import BaseModel
from app.services.validation_service import ValidationService

router = APIRouter()
validation_service = ValidationService()

class ValidationRequest(BaseModel):
    pdb_content: str
    label: str = "molecule_from_front"

@router.post("/check", summary="Zvaliduje aktuální stav molekuly z frontendu")
async def check_molecule(request: ValidationRequest):
    return validation_service.validate_pdb_content(request.pdb_content, request.label)