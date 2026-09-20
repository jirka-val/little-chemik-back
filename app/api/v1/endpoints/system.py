# app/api/v1/endpoints/system.py
from fastapi import APIRouter

from app.core.logging import console_log_buffer

router = APIRouter()


@router.get("/logs", summary="Poslední stavové zprávy pro Console panel (polling)")
async def get_logs(since_id: int = 0, limit: int = 200):
    """
    Vrací jen krátké, obecné stavové zprávy (viz app.core.logging.console_logger) -
    NE kompletní aplikační log. Frontend pollováním (viz Console panel v sidebaru)
    postupně dohání aktuální stav bez opakovaného stahování celého bufferu.
    `last_id` v odpovědi se pošle jako `since_id` v příštím requestu.
    """
    entries = console_log_buffer.get_since(since_id=since_id, limit=limit)
    last_id = entries[-1]["id"] if entries else since_id
    return {"entries": entries, "last_id": last_id}
