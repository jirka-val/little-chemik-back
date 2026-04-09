from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


def test_health_check_root():
    """
    Testuje, zda kořenový endpoint '/' vrací status 200 a správný JSON.
    """
    response = client.get("/")
    assert response.status_code == 200

    data = response.json()

    assert data["status"] == "online"
    assert "version" in data