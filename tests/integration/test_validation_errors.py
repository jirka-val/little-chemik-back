"""
FastAPI 422 validace vstupu (app/core/exceptions.py::validation_exception_handler).
"""

import pytest

pytestmark = pytest.mark.integration

_URL = "/api/simulation/nope/amber-mdin?preview=true"


@pytest.mark.parametrize("raw", ['{"temp0": Infinity}', '{"temp0": NaN}', '{"duration_ns": -Infinity}'])
def test_non_finite_numbers_are_422_not_500(client, raw):
    # Python json NaN/Infinity přijme; ozvěna vstupu v 422 je dřív shodila na 500
    response = client.post(_URL, content=raw, headers={"Content-Type": "application/json"})
    assert response.status_code == 422
    assert response.json()["detail"][0]["loc"][-1] in ("temp0", "duration_ns")


def test_invalid_step_points_at_the_step(client):
    response = client.post(_URL, json={"protocol": True, "steps": [
        {"name": "h1", "kind": "heating"}, {"name": "h2", "kind": "heating", "temp_target": 5000},
    ]})
    assert response.status_code == 422
    assert response.json()["detail"][0]["loc"] == ["body", "steps", 1, "temp_target"]
