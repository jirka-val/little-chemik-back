from fastapi import Request
from fastapi.responses import JSONResponse
import logging

logger = logging.getLogger("api")

class AppBaseException(Exception):
    def __init__(self, message: str, status_code: int = 500):
        self.message = message
        self.status_code = status_code

class MoleculeNotFoundError(AppBaseException):
    def __init__(self, pdb_code: str):
        super().__init__(f"Molekula s kódem {pdb_code} nebyla nalezena na serveru RCSB.", 404)

class ProcessingError(AppBaseException):
    def __init__(self, detail: str):
        super().__init__(f"Chyba při analýze molekuly: {detail}", 422)

async def app_exception_handler(request: Request, exc: AppBaseException):
    logger.error(f"Chyba aplikace na {request.url.path}: {exc.message}")
    return JSONResponse(
        status_code=exc.status_code,
        content={"error": exc.__class__.__name__, "message": exc.message}
    )