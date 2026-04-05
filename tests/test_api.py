import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)
shared_workspace_id = None

from app.services.structure.hydrogenation import HydrogenationService



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


def test_11_inspect_hydrogenation_service():
    service = HydrogenationService()
    # Získáme všechny metody, které nejsou interní (nezačínají _)
    methods = [method for method in dir(service) if callable(getattr(service, method)) and not method.startswith("_")]

    print(f"\n--- DEBUG: Dostupné metody v HydrogenationService ---")
    for m in methods:
        print(f" -> {m}")

    # Tento assert nám v testu 'test_5' padá. Tady zjistíme proč.
    assert "add_hydrogen_atoms" in methods, f"Chyba: HydrogenationService postrádá metodu add_hydrogen_atoms! Dostupné: {methods}"