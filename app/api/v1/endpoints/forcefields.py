from fastapi import APIRouter, HTTPException
from app.services.pdb_service import PDBService
from app.services.validation.forcefield import ForceFieldValidator
from app.workspaces.manager import workspace_manager

router = APIRouter()
pdb_service = PDBService()
ff_validator = ForceFieldValidator()

@router.get("/{workspace_id}")
async def get_my_forcefields(workspace_id: str):
    # 1. Existuje ten workspace?
    if not workspace_manager.workspace_exists(workspace_id):
        raise HTTPException(status_code=404, detail="Soubor nenalezen")

    # 2. Načteme PDB ze souboru
    path = workspace_manager.get_file_path(workspace_id)
    with open(path, "r", encoding="utf-8") as f:
        pdb_content = f.read()

    # 3. Zjistíme, co v tom je za chemii (D, R, P...)
    types = pdb_service.get_molecule_types(pdb_content)

    # 4. Stáhneme ty správné FF
    ffs = ff_validator.get_matching_forcefields(types)

    return {
        "detected_types": types,
        "forcefields": ffs
    }