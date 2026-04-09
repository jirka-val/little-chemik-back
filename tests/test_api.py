import pytest
import asyncio
import time
from httpx import AsyncClient, ASGITransport
from fastapi.testclient import TestClient
from app.main import app
from app.services.structure.hydrogenation import HydrogenationService

client = TestClient(app)
shared_workspace_id = None

def test_1_upload_valid_pdb():
    """Uploads a valid .pdb file and checks if workspace_id is returned."""
    dummy_pdb_content = b"ATOM      1  N   ALA A   1      -1.011   1.455  -0.082  1.00  0.00           N\n"
    files = {"file": ("test_molecule.pdb", dummy_pdb_content, "chemical/x-pdb")}

    response = client.post("/api/molecules/upload", files=files)
    assert response.status_code == 200

    data = response.json()
    assert "workspace_id" in data
    assert data["filename"] == "test_molecule.pdb"

    global shared_workspace_id
    shared_workspace_id = data["workspace_id"]

def test_2_upload_invalid_file():
    """Tries to upload an invalid format - server should reject with 400."""
    dummy_text = b"This is not a pdb, this is just text"
    files = {"file": ("invalid.txt", dummy_text, "text/plain")}

    response = client.post("/api/molecules/upload", files=files)
    assert response.status_code == 400
    # Adjusted to match the likely English error message
    assert "only .pdb files" in response.json()["detail"].lower() or "pdb" in response.json()["detail"].lower()

def test_3_download_workspace():
    """Tries to download data via the fast download endpoint for Molstar."""
    assert shared_workspace_id is not None
    response = client.get(f"/api/download/{shared_workspace_id}")

    assert response.status_code == 200
    assert response.headers["content-type"] == "chemical/x-pdb"
    assert "ATOM" in response.text

def test_4_download_nonexistent_workspace():
    """Verifies that the server returns 404 for a non-existent workspace."""
    response = client.get("/api/download/non-existent-uuid-12345")
    assert response.status_code == 404

def test_5_full_integration_fetch_and_hydrogenate():
    """Fetches 1CRN from PDB and sends it for hydrogenation."""
    response_fetch = client.get("/api/molecules/fetch-pdb/1crn")
    assert response_fetch.status_code == 200

    ws_id = response_fetch.json()["workspace_id"]
    response_hydrogens = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={"ph": 7.4})

    assert response_hydrogens.status_code == 200
    assert "Hydrogens successfully added" in response_hydrogens.json()["message"]

def test_6_fetch_invalid_pdb_code():
    """Verifies that an invalid PDB code returns 404."""
    response = client.get("/api/molecules/fetch-pdb/NON_EXISTENT_9999")
    assert response.status_code == 404
    assert "not exist" in response.json()["detail"].lower() or "error" in response.json()["detail"].lower()

def test_7_sequence_analysis_integration():
    """Uploads valid PDB and verifies that the sequence analysis is correct."""
    dummy_pdb_content = b"ATOM      1  N   ALA A   1      -1.011   1.455  -0.082  1.00  0.00           N\n"
    files = {"file": ("alanine.pdb", dummy_pdb_content, "chemical/x-pdb")}
    response_upload = client.post("/api/molecules/upload", files=files)
    assert response_upload.status_code == 200

    ws_id = response_upload.json()["workspace_id"]
    response_analysis = client.get(f"/api/analysis/sequence/{ws_id}")

    assert response_analysis.status_code == 200
    data = response_analysis.json()
    assert "sequence" in data
    assert "A" in data["sequence"]["chains"]
    assert data["sequence"]["chains"]["A"]["tokens"][0]["pdb_resname"] == "ALA"

def test_8_add_hydrogens_to_gibberish():
    """Tries to add hydrogens to gibberish file and verifies 500 error."""
    gibberish_content = b"This certainly is not a valid chemical structure! @#$%"
    files = {"file": ("gibberish.pdb", gibberish_content, "chemical/x-pdb")}
    upload_res = client.post("/api/molecules/upload", files=files)
    ws_id = upload_res.json()["workspace_id"]

    response = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={"ph": 7.0})
    assert response.status_code == 500
    assert "Structure modification failed" in response.json()["detail"]

def test_9_3dvz_sequence_missing_atoms_full_debug():
    """Diagnostic test for 3DVZ to verify token states."""
    response_fetch = client.get("/api/molecules/fetch-pdb/3dvz")
    assert response_fetch.status_code == 200
    workspace_id = response_fetch.json()["workspace_id"]

    response_analysis = client.get(f"/api/analysis/sequence/{workspace_id}")
    assert response_analysis.status_code == 200
    data = response_analysis.json()

    chain_a_tokens = data["sequence"]["chains"]["A"]["tokens"]
    u2647 = next((t for t in chain_a_tokens if str(t["resseq"]) == "2647"), None)

    assert u2647 is not None
    assert len(u2647["missing_atoms"]) > 0

def test_10_3dvz_rna_variants_and_atoms_detection():
    """Tests if the backend correctly identifies RNA variants in 3DVZ."""
    response_fetch = client.get("/api/molecules/fetch-pdb/3dvz")
    assert response_fetch.status_code == 200
    workspace_id = response_fetch.json()["workspace_id"]

    response_analysis = client.get(f"/api/analysis/sequence/{workspace_id}")
    assert response_analysis.status_code == 200
    data = response_analysis.json()

    tokens = data["sequence"]["chains"]["A"]["tokens"]

    # Check first residue (U 2647 -> RU5)
    u2647 = next((t for t in tokens if str(t["resseq"]) == "2647"), None)
    assert u2647 is not None
    assert u2647["ff_resname"] == "RU5"
    assert u2647["known"] is True

    # Check internal residue (G 2648 -> RG)
    g2648 = next((t for t in tokens if str(t["resseq"]) == "2648"), None)
    assert g2648 is not None
    assert g2648["ff_resname"] == "RG"
    assert g2648["known"] is True


@pytest.mark.asyncio
async def test_11_concurrency_no_blocking():
    """
    Verifies that a heavy CPU-bound request (sequence analysis) does not block
    the FastAPI event loop, allowing a light request (download) to be processed instantly.
    """
    # Získáme validní molekulu pomocí synchronního klienta, ať máme co analyzovat.
    response_fetch = client.get("/api/molecules/fetch-pdb/1crn")
    assert response_fetch.status_code == 200
    ws_id = response_fetch.json()["workspace_id"]

    # Použijeme asynchronní klient pro simulaci 2 uživatelů naráz
    async with AsyncClient(transport=ASGITransport(app=app), base_url="http://test") as ac:

        async def fetch_heavy():
            """Uživatel 1 - Požádá o těžký výpočet analýzy sekvence."""
            start = time.time()
            res = await ac.get(f"/api/analysis/sequence/{ws_id}")
            duration = time.time() - start
            return duration, res.status_code

        async def fetch_light():
            """Uživatel 2 - Zkusí si jen bleskově stáhnout soubor, zatímco se počítá analýza."""
            # Malé zpoždění, aby těžký úkol zaručeně odstartoval první
            await asyncio.sleep(0.05)
            start = time.time()
            res = await ac.get(f"/api/download/{ws_id}")
            duration = time.time() - start
            return duration, res.status_code

        # Spustíme oba dotazy současně
        results = await asyncio.gather(fetch_heavy(), fetch_light())

        heavy_duration, heavy_status = results[0]
        light_duration, light_status = results[1]

        assert heavy_status == 200
        assert light_status == 200

        # Můžeme si to i vypsat pro kontrolu při běhu pytestu s příznakem -s
        print(f"\nHEAVY TASK (CPU) duration: {heavy_duration:.4f}s")
        print(f"LIGHT TASK (I/O) duration: {light_duration:.4f}s")

        # Lehký úkol se musí dokončit buď dříve než těžký úkol, nebo v čase menším než 0.1s.
        # Kdyby backend stále blokoval event loop, lehký úkol by musel čekat na dokončení těžkého.
        assert light_duration < heavy_duration or light_duration < 0.1