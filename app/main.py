from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from app.services.pdb_service import PDBService
from app.core.logging import setup_logging
import logging

setup_logging()
logger = logging.getLogger("api")

app = FastAPI(title="Nucleic Acid Analysis API", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

pdb_service = PDBService()

@app.get("/")
async def root():
    return {"status": "online", "service": "Nucleic Acid Analysis"}

@app.get("/api/process/{pdb_code}")
async def process_molecule(pdb_code: str):
    try:
        content = await pdb_service.fetch_pdb_content(pdb_code.lower())
        return content
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception("Neočekávaná chyba při zpracování molekuly")
        raise HTTPException(status_code=500, detail="Interní chyba serveru.")