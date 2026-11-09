"""
Velké JSON odpovědi mimo event loop.

FastAPI u vráceného dictu nejdřív projde každý objekt přes jsonable_encoder
(čistý Python) a pak teprve volá json.dumps - obojí v event loopu. U sekvence
připravené struktury (3,7 MB) to trvalo 240 ms, po které server nikomu
neodpovídal; json.dumps v C stejný obsah zakóduje za 20 ms (tests/load).

large_json_response() data zakóduje stejně jako Starlette JSONResponse, ale ve
threadpoolu. Data musí být čisté JSON typy (dict, list, str, int, float,
bool, None) - jinak json.dumps vyhodí TypeError.
"""

from __future__ import annotations

import json
from typing import Any

from fastapi.concurrency import run_in_threadpool
from fastapi.responses import Response


def _encode(content: Any) -> bytes:
    # Stejné volby jako starlette.responses.JSONResponse.render
    return json.dumps(content, ensure_ascii=False, allow_nan=False, indent=None, separators=(",", ":")).encode("utf-8")


async def large_json_response(content: Any, status_code: int = 200) -> Response:
    body = await run_in_threadpool(_encode, content)
    return Response(content=body, status_code=status_code, media_type="application/json")
