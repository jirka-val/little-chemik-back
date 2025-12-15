from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
import os
import urllib.request

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


def simple_pdb_cut(pdb_text: str, res_start: int, res_end: int) -> str:
    """
    Jednoduchá funkce pro ořezání PDB souboru (mazání řádků).
    Nahrazuje složitou logiku, dokud nemáme ChemPy.
    """
    new_lines = []
    cut_count = 0

    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                res_seq = int(line[22:26].strip())

                if res_start <= res_seq <= res_end:
                    cut_count += 1
                    continue 

            except ValueError:
                pass

        new_lines.append(line)

    print(f"DEBUG: Ořezáno {cut_count} atomů v rozsahu {res_start}-{res_end}")
    return "\n".join(new_lines)


@app.get("/")
def read_root():
    return {"message": "Backend běží! Použij endpoint /api/process/{pdb_code}"}


@app.get("/api/process/{pdb_code}", response_class=PlainTextResponse)
async def process_molecule(pdb_code: str, res1: int = 0, res2: int = 0):
    """
    Hlavní funkce:
    1. Stáhne PDB soubor (pokud ho nemáme).
    2. Načte ho.
    3. Pokud jsou zadány parametry res1 a res2, provede ořezání.
    4. Vrátí text PDB souboru.
    """
    pdb_filename = f"{pdb_code}.pdb"
    print(f"--- Požadavek: {pdb_code} | Rozsah mazání: {res1} - {res2} ---")

    try:
        if not os.path.exists(pdb_filename):
            print(f"Stahuji {pdb_code} z internetu...")
            url = f"[https://files.rcsb.org/download/](https://files.rcsb.org/download/){pdb_code}.pdb"
            urllib.request.urlretrieve(url, pdb_filename)

        with open(pdb_filename, "r") as f:
            pdb_content = f.read()

        if res1 > 0 and res2 >= res1:
            modified_data = simple_pdb_cut(pdb_content, res1, res2)
            return modified_data
        else:
            return pdb_content

    except Exception as e:
        print(f"CHYBA: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
