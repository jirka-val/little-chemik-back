from fastapi import APIRouter, HTTPException
from app.services.validation_service import ValidationService

router = APIRouter()
validation_service = ValidationService()


@router.get("/{pdb_code}", summary="Ověří, zda je molekula připravena pro HPC simulaci")
async def validate_molecule(pdb_code: str):
    report = await validation_service.validate_molecule(pdb_code)

    if "error" in report:
        raise HTTPException(status_code=404, detail=report["error"])

    return report