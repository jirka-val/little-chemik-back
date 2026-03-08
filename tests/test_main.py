from fastapi.testclient import TestClient
# Importujeme vaši aplikaci z vašeho kódu
from app.main import app

# Vytvoříme falešného klienta (robota), který bude na aplikaci sahat
client = TestClient(app)


def test_health_check_root():
    """
    Testuje, zda kořenový endpoint '/' vrací status 200 a správný JSON.
    """
    # 1. AKCE: Falešný klient pošle GET požadavek na kořenovou URL
    response = client.get("/")

    # 2. OVĚŘENÍ (Asserts): Zkontrolujeme, co server vrátil

    # Tvrdím, že HTTP status musí být 200 (OK)
    assert response.status_code == 200

    # Získám si JSON odpověď jako Python slovník
    data = response.json()

    # Tvrdím, že v odpovědi musí být klíč 'status' s hodnotou 'online'
    assert data["status"] == "online"

    # Tvrdím, že tam musí být i informace o verzi (ověříme, že klíč existuje)
    assert "version" in data