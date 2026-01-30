from fastapi import APIRouter, Body, HTTPException
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any

from app.services.analysis_service import build_sequence_tokens

router = APIRouter()

class SequenceRequest(BaseModel):
    pdb_text: str = Field(..., description="Raw PDB text (ATOM/HETATM records are enough).")
    chain: Optional[str] = Field(None, description="Chain identifier, e.g. 'A'")
    fill_gaps: bool = Field(True, description="Insert gap tokens for missing residue numbers.")

@router.post("/sequence", summary="Vrátí sekvenci residuí pro vykreslení řetězce")
async def sequence(payload: SequenceRequest = Body(...)) -> Dict[str, Any]:
    if not payload.pdb_text.strip():
        raise HTTPException(status_code=400, detail="pdb_text is empty")

    return build_sequence_tokens(
        pdb_text=payload.pdb_text,
        chain=payload.chain,
        fill_gaps=payload.fill_gaps,
    )
