"""
Middleware pro deník akcí (viz app/services/incidents/history.py).

Ke každému API požadavku, který se týká existujícího workspace (UUID v cestě
nebo `workspace_id` v JSON těle), zapíše do jeho deníku čas, cestu, JSON
nastavení, status a dobu trvání. Před požadavky, které mohou strukturu
změnit (cokoli kromě GET), uloží snímek structure.pdb. Serverovou chybu
zapíše do errors.jsonl i s tracebackem, který sem předá app_exception_handler
přes request.state.

Při INCIDENT_REPORTS_ENABLED=false nedělá nic.
"""

import json
import time
import traceback

from fastapi import Request
from fastapi.concurrency import run_in_threadpool
from starlette.middleware.base import BaseHTTPMiddleware

from app.services.incidents import history

_SKIP_PREFIXES = ("/api/incidents",)


class IncidentJournalMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        path = request.url.path
        # OPTIONS = CORS preflight prohlížeče, nic neříká o krocích uživatele.
        if (
            not history.enabled()
            or request.method == "OPTIONS"
            or not path.startswith("/api/")
            or path.startswith(_SKIP_PREFIXES)
        ):
            return await call_next(request)

        body = None
        body_size = 0
        if request.method != "GET" and request.headers.get("content-type", "").startswith("application/json"):
            raw = await request.body()
            body_size = len(raw)
            try:
                body = json.loads(raw) if raw else None
            except ValueError:
                body = "[invalid JSON]"
        elif request.method != "GET" and "multipart/form-data" in request.headers.get("content-type", ""):
            body = "[multipart upload]"

        workspace_id = history.workspace_id_from(path, body)
        if workspace_id and request.method != "GET":
            # "/api/sidechains/start/<uuid>" -> "sidechains-start"
            label = path.removeprefix("/api/").replace(workspace_id, "").strip("/").replace("/", "-")
            await run_in_threadpool(history.snapshot_structure, workspace_id, label)

        started = time.perf_counter()
        entry = {
            "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "method": request.method,
            "path": path,
            "query": str(request.url.query) or None,
            "body": body,
            "body_size": body_size,
        }
        try:
            response = await call_next(request)
        except Exception as exc:
            if workspace_id:
                entry.update(status=500, duration_s=round(time.perf_counter() - started, 3))
                await run_in_threadpool(history.record_request, workspace_id, entry)
                await run_in_threadpool(history.record_error, workspace_id, {
                    "path": path,
                    "status": 500,
                    "message": str(exc),
                    "traceback": "".join(traceback.format_exception(exc)),
                })
            raise

        if workspace_id:
            entry.update(status=response.status_code, duration_s=round(time.perf_counter() - started, 3))
            await run_in_threadpool(history.record_request, workspace_id, entry)
            if response.status_code >= 500:
                await run_in_threadpool(history.record_error, workspace_id, {
                    "path": path,
                    "status": response.status_code,
                    "message": getattr(request.state, "incident_error_message", None),
                    "traceback": getattr(request.state, "incident_traceback", None),
                })
        return response
