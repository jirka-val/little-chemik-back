import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)
shared_workspace_id = None


def test_1_upload_valid_pdb():
    """Zkusí nahrát validní .pdb soubor a zkontroluje, zda vrátí workspace_id."""
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
    """Zkusí nahrát špatný formát - server by ho měl odmítnout s chybou 400."""
    dummy_text = b"Tohle neni pdbcko, tohle je text"
    files = {"file": ("zavirak.txt", dummy_text, "text/plain")}

    response = client.post("/api/molecules/upload", files=files)
    assert response.status_code == 400
    assert "pouze .pdb soubory" in response.json()["detail"]


def test_3_download_workspace():
    """Pokusí se stáhnout data přes rychlý download endpoint pro Molstar."""
    assert shared_workspace_id is not None
    response = client.get(f"/api/download/{shared_workspace_id}")

    assert response.status_code == 200
    assert response.headers["content-type"] == "chemical/x-pdb"
    assert "ATOM" in response.text


def test_4_download_nonexistent_workspace():
    """Ověří, že server vrátí 404 při požadavku na stažení neexistujícího workspace."""
    response = client.get("/api/download/neexistujici-uuid-12345")
    assert response.status_code == 404


def test_5_full_integration_fetch_and_hydrogenate():
    """Zkusí stáhnout reálnou molekulu 1CRN z Protein Data Bank a pošle ji k přidání vodíků."""
    response_fetch = client.get("/api/molecules/fetch-pdb/1crn")
    assert response_fetch.status_code == 200

    ws_id = response_fetch.json()["workspace_id"]
    response_hydrogens = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={"ph": 7.4})

    assert response_hydrogens.status_code == 200
    assert "Vodíky byly úspěšně doplněny" in response_hydrogens.json()["message"]


def test_6_fetch_invalid_pdb_code():
    """Ověří, že zadání neexistujícího PDB kódu nesestřelí server a vrátí 404."""
    response = client.get("/api/molecules/fetch-pdb/NEEXISTUJICI_KOD_9999")
    assert response.status_code == 404
    assert "neexistuje" in response.json()["detail"].lower()


def test_7_sequence_analysis_integration():
    """Nahraje validní kousek PDB a ověří, že API správně provede analýzu sekvence."""
    dummy_pdb_content = b"ATOM      1  N   ALA A   1      -1.011   1.455  -0.082  1.00  0.00           N\n"
    files = {"file": ("alanin.pdb", dummy_pdb_content, "chemical/x-pdb")}
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
    """Pokusí se přidat vodíky do nesmyslného souboru a ověří, že server vrátí čistou chybu 500."""
    gibberish_content = b"Tohle urcite neni spravna chemicka struktura! @#$%"
    files = {"file": ("nesmysl.pdb", gibberish_content, "chemical/x-pdb")}
    upload_res = client.post("/api/molecules/upload", files=files)
    ws_id = upload_res.json()["workspace_id"]

    response = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={"ph": 7.0})
    assert response.status_code == 500
    assert "Chyba při úpravě struktury" in response.json()["detail"]


def test_9_3dvz_sequence_missing_atoms_full_debug():
    """
    Diagnostický test pro 3DVZ k ověření vnitřních stavů tokenu.
    """
    # 1. Stažení 3DVZ
    response_fetch = client.get("/api/molecules/fetch-pdb/3dvz")
    assert response_fetch.status_code == 200
    workspace_id = response_fetch.json()["workspace_id"]

    # 2. Volání analýzy sekvence
    response_analysis = client.get(f"/api/analysis/sequence/{workspace_id}")
    assert response_analysis.status_code == 200
    data = response_analysis.json()

    # 3. Získání tokenu pro U 2647
    chain_a_tokens = data["sequence"]["chains"]["A"]["tokens"]
    u2647 = next((t for t in chain_a_tokens if str(t["resseq"]) == "2647"), None)

    assert u2647 is not None, "Reziduum A 2647 nebylo v sekvenci nalezeno"

    # --- KOMPLETNÍ DIAGNOSTIKA ---
    print(f"\n" + "=" * 50)
    print(f"DEBUG DATA PRO REZIDUUM: {u2647['pdb_resname']} {u2647['resseq']}")
    print(f"=" * 50)
    print(f"1. Detekovaná FF varianta (ff_resname): {u2647.get('ff_resname')}")
    print(f"2. Skupina (group): {u2647.get('group')}")
    print(f"3. Je známo v JSONu (known): {u2647.get('known')}")
    print(f"4. Pozice v řetězci (position): {u2647.get('position')}")
    print(f"5. Atomy v PDB ({len(u2647['atoms'])}): {u2647['atoms']}")
    print(f"6. Chybějící atomy: {u2647['missing_atoms']}")

    # Kontrola, zda jsou přítomny aspoň základní atomy pro identifikaci
    has_o5_prime = "O5'" in u2647['atoms']
    print(f"7. Obsahuje O5': {has_o5_prime}")

    # Kontrola varování pro tento workspace (může obsahovat důvod selhání vyhledávání v JSONu)
    warnings = data["sequence"]["chains"]["A"].get("warnings", [])
    print(f"8. Varování pro řetězec A: {warnings}")
    print("=" * 50)

    # 4. KONTROLA
    assert len(u2647["missing_atoms"]) > 0, \
        f"Chyba: Reziduum {u2647['ff_resname']} by mělo hlásit chybějící vodíky a HO5'!"


def test_10_3dvz_rna_variants_and_atoms_detection():
    """
    Testuje, zda backend správně identifikuje RNA varianty (RU5 vs RU)
    ve struktuře 3DVZ a zda vidí chybějící atomy.
    """
    # 1. Stažení 3DVZ z PDB
    response_fetch = client.get("/api/molecules/fetch-pdb/3dvz")
    assert response_fetch.status_code == 200
    workspace_id = response_fetch.json()["workspace_id"]

    # 2. Analýza sekvence
    response_analysis = client.get(f"/api/analysis/sequence/{workspace_id}")
    assert response_analysis.status_code == 200
    data = response_analysis.json()

    # Získáme tokeny pro řetězec A
    tokens = data["sequence"]["chains"]["A"]["tokens"]

    # --- KONTROLA PRVNÍHO REZIDUA (U 2647 -> RU5) ---
    u2647 = next((t for t in tokens if str(t["resseq"]) == "2647"), None)
    assert u2647 is not None

    print(f"\nDiagnostika {u2647['pdb_resname']} {u2647['resseq']}:")
    print(f"  - FF varianta: {u2647['ff_resname']}")
    print(f"  - Je známo v JSONu: {u2647['known']}")
    print(f"  - Chybějící atomy: {u2647['missing_atoms']}")

    # Ověření, že se správně vybrala varianta RU5 (protože je první)
    assert u2647["ff_resname"] == "RU5", f"Očekáváno RU5, ale nalezeno {u2647['ff_resname']}"
    assert u2647["known"] is True, "Reziduum RU5 by mělo být v JSONu nalezeno (known: True)"
    # RU5 v 3DVZ nemá vodíky ani terminální HO5', seznam by neměl být prázdný
    assert len(u2647["missing_atoms"]) > 0

    # --- KONTROLA VNITŘNÍHO REZIDUA (G 2648 -> RG) ---
    g2648 = next((t for t in tokens if str(t["resseq"]) == "2648"), None)
    assert g2648 is not None

    print(f"\nDiagnostika {g2648['pdb_resname']} {g2648['resseq']}:")
    print(f"  - FF varianta: {g2648['ff_resname']}")
    print(f"  - Je známo v JSONu: {g2648['known']}")

    # Vnitřní reziduum by mělo být RG (ne RG5 ani RG3)
    assert g2648["ff_resname"] == "RG"
    assert g2648["known"] is True