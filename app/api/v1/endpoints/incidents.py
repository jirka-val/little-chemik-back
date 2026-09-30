"""
Hlášení chyb od uživatelů (viz app/services/incidents/store.py).

- GET  /api/incidents/config          - je funkce zapnutá? (frontend podle toho
                                        ukáže tlačítko "!" a nabídku reportu)
- POST /api/incidents                 - odeslání reportu; jen se souhlasem
- GET  /api/incidents                 - (admin) seznam reportů
- GET  /api/incidents/{id}/download   - (admin) celý report jako ZIP

Admin endpointy chrání Settings.ADMIN_TOKEN (hlavička X-Admin-Token), stejně
jako přeřazování force fieldů.
"""

import logging
from typing import Any, Dict, List, Literal, Optional

from fastapi import APIRouter, Header
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import Response
from pydantic import BaseModel, Field

from app.core.config import settings
from app.core.exceptions import AppBaseException, BadRequestError, ForbiddenError, NotFoundError
from app.services.incidents import history, store

logger = logging.getLogger(__name__)
router = APIRouter()


class IncidentReportsDisabledError(AppBaseException):
    status_code = 404
    code = "incident_reports_disabled"

    def __init__(self):
        super().__init__("Error reporting is turned off on this server.")


class IncidentStorageFullError(AppBaseException):
    status_code = 507
    code = "incident_storage_full"

    def __init__(self):
        super().__init__("The server cannot store more error reports right now. Please contact the developers directly.")


class IncidentReport(BaseModel):
    workspace_id: Optional[str] = None
    kind: Literal["error", "manual"]
    description: str = Field("", max_length=20_000)
    consent: bool
    error: Optional[Dict[str, Any]] = None
    client_log: List[Any] = Field(default_factory=list, max_length=1000)
    page: Dict[str, Any] = Field(default_factory=dict)
    # PNG snímek vieweru jako data URL (~ do 10 MB).
    screenshot: Optional[str] = Field(None, max_length=14_000_000)


def _require_admin(token: str) -> None:
    if not settings.ADMIN_TOKEN or token != settings.ADMIN_TOKEN:
        raise ForbiddenError()


@router.get("/config", summary="Je hlášení chyb zapnuté?")
async def incidents_config():
    return {"enabled": history.enabled()}


@router.post("", summary="Odešle report chyby (se souhlasem uživatele)")
async def submit_incident(report: IncidentReport):
    if not history.enabled():
        raise IncidentReportsDisabledError()
    if not report.consent:
        raise BadRequestError("The report can only be stored with your consent.")
    if report.kind == "manual" and not report.description.strip():
        raise BadRequestError("Please describe what is wrong.")

    try:
        incident_id = await run_in_threadpool(
            store.create_incident,
            workspace_id=report.workspace_id,
            kind=report.kind,
            description=report.description,
            error=report.error,
            client_log=report.client_log,
            page=report.page,
            screenshot=report.screenshot,
        )
    except store.IncidentStorageFull:
        raise IncidentStorageFullError()
    return {"incident_id": incident_id}


@router.get("", summary="(admin) Seznam nahlášených chyb")
async def list_incidents(x_admin_token: str = Header(default="")):
    _require_admin(x_admin_token)
    incidents = await run_in_threadpool(store.list_incidents)
    used = await run_in_threadpool(store.storage_used_bytes)
    return {
        "incidents": incidents,
        "storage_used_mb": round(used / (1024 * 1024), 1),
        "storage_limit_mb": settings.INCIDENTS_MAX_TOTAL_MB,
    }


@router.get("/{incident_id}/download", summary="(admin) Stáhne report jako ZIP")
async def download_incident(incident_id: str, x_admin_token: str = Header(default="")):
    _require_admin(x_admin_token)
    data = await run_in_threadpool(store.incident_zip, incident_id)
    if data is None:
        raise NotFoundError(f"Incident {incident_id} not found.")
    return Response(
        content=data,
        media_type="application/zip",
        headers={"Content-Disposition": f"attachment; filename=incident_{incident_id}.zip"},
    )
