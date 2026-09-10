"""Health-check endpoint. Nahrazuje tests/test_main.py."""

import pytest

pytestmark = pytest.mark.integration


def test_root_health_check(client):
    response = client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "online"
    assert "version" in data
