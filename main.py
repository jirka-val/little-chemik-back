from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
import os
import urllib.request

app = FastAPI()

# 1. Povolení CORS - nutné pro komunikaci s Vite (port 5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Backend běží! Použij /api/process/{pdb_code} pro načtení molekuly."}

# 2. Změna z /api/load/ na /api/process/, aby to odpovídalo frontendu
@app.get("/api/process/{pdb_code}", response_class=PlainTextResponse)
async def load_molecule(pdb_code: str):
    """
    Stáhne PDB soubor z RCSB (pokud ho ještě nemáme lokálně) a vrátí jeho obsah.
    """
    pdb_filename = f"{pdb_code}.pdb"
    print(f"--- Požadavek na načtení: {pdb_code} ---")

    try:
        # Pokud soubor nemáme, stáhneme ho
        if not os.path.exists(pdb_filename):
            print(f"Stahuji {pdb_code} z RCSB...")
            # Opravená URL (bez markdown závorek)
            url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
            urllib.request.urlretrieve(url, pdb_filename)

        # Načtení souboru
        with open(pdb_filename, "r") as f:
            pdb_content = f.read()

        return pdb_content

    except Exception as e:
        print(f"CHYBA: {e}")
        # Pokud PDB neexistuje na serveru RCSB, vrátíme chybu
        raise HTTPException(status_code=404, detail=f"Molekula {pdb_code} nebyla nalezena.")

if __name__ == "__main__":
    import uvicorn
    # Spuštění na portu 8000
    uvicorn.run(app, host="0.0.0.0", port=8000)