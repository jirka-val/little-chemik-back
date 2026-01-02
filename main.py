from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel
import os
import urllib.request

app = FastAPI()

# Povolení komunikace s frontendem (CORS)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Model pro data, která přijdou z frontendu
# Pydantic zajistí, že request.pdb_text bude string a request.residue_number integer
class DeleteRequest(BaseModel):
    pdb_text: str
    residue_number: int

def delete_specific_residue(pdb_text: str, residue_to_delete: int) -> str:
    """
    Funkce projde PDB text a odstraní řádky atomů, které patří
    specifickému číslu rezidua (residue_to_delete).
    """
    new_lines = []
    cut_count = 0

    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                res_seq = int(line[22:26].strip())

                if res_seq == residue_to_delete:
                    cut_count += 1
                    continue

            except ValueError:
                pass

        new_lines.append(line)

    print(f"DEBUG: Smazáno {cut_count} atomů rezidua číslo {residue_to_delete}")
    return "\n".join(new_lines)


@app.get("/")
def read_root():
    return {"message": "Backend běží! Použij /api/load/{pdb_code} pro načtení a POST /api/delete pro úpravu."}


@app.get("/api/load/{pdb_code}", response_class=PlainTextResponse)
async def load_molecule(pdb_code: str):
    """
    1. Prvotní načtení: Stáhne PDB soubor (pokud ho nemáme) a vrátí ho.
    Tuto funkci zavolá frontend jen jednou na začátku.
    """
    pdb_filename = f"{pdb_code}.pdb"
    print(f"--- Načítání originálu: {pdb_code} ---")

    try:
        if not os.path.exists(pdb_filename):
            print(f"Stahuji {pdb_code} z internetu...")
            url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
            urllib.request.urlretrieve(url, pdb_filename)

        with open(pdb_filename, "r") as f:
            pdb_content = f.read()

        return pdb_content

    except Exception as e:
        print(f"CHYBA při stahování: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/delete", response_class=PlainTextResponse)
async def delete_residue_endpoint(request: DeleteRequest):
    print(f"--- Požadavek na smazání rezidua: {request.residue_number} ---")

    try:
        modified_data = delete_specific_residue(request.pdb_text, request.residue_number)
        return modified_data

    except Exception as e:
        print(f"CHYBA při mazání: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)