"""
Jednotná chybová obálka pro celé API.

Cíl: FRONTEND MUSÍ MÍT NA COKOLIV, CO SPADNE, JEDEN SPOLEHLIVÝ TVAR ODPOVĚDI:

    {"error": "<PythonTřídaVýjimky>", "code": "<stabilní_strojový_identifikátor>",
     "message": "<lidsky čitelný text>", ...volitelná strukturovaná data}

`code` je ten identifikátor, na který se má frontend dívat (ne `error` -
název Python třídy se může kdykoliv přejmenovat při refaktoru, `code` je
API kontrakt a musí zůstat stabilní). Endpointy by NEMĚLY vyhazovat
`fastapi.HTTPException` přímo (s výjimkou FastAPI/pydantic vlastní 422
validace vstupu, kterou nekontrolujeme) - místo toho vždy jednu z tříd
níže, nebo novou podtřídu `AppBaseException`, když se objeví další
opakující se případ.
"""

from fastapi import Request
from fastapi.responses import JSONResponse
import logging

logger = logging.getLogger("api")


class AppBaseException(Exception):
    status_code: int = 500
    code: str = "internal_error"

    def __init__(self, message: str, status_code: int | None = None, code: str | None = None,
                 payload: dict | None = None):
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        if code is not None:
            self.code = code
        self.payload = payload or {}


class WorkspaceNotFoundError(AppBaseException):
    """Workspace s daným ID neexistuje - typicky uplynulá/smazaná session nebo překlep v ID."""
    status_code = 404
    code = "workspace_not_found"

    def __init__(self, message: str = "Workspace not found. Have you uploaded a file first?"):
        super().__init__(message)


class NotFoundError(AppBaseException):
    """Obecné 404 pro konkrétní zdroj uvnitř existujícího workspace (soubor, reziduum...)."""
    status_code = 404
    code = "not_found"


class RemoteMoleculeNotFoundError(AppBaseException):
    """Zadaný PDB kód v RCSB databázi neexistuje."""
    status_code = 404
    code = "remote_molecule_not_found"

    def __init__(self, pdb_code: str):
        super().__init__(f"Molecule '{pdb_code}' does not exist in the RCSB database.",
                          payload={"pdb_code": pdb_code})


class ExternalServiceError(AppBaseException):
    """Externí služba (RCSB, IDA API) neodpověděla nebo vrátila neočekávanou chybu."""
    status_code = 502
    code = "external_service_error"


class BadRequestError(AppBaseException):
    """Vstup od uživatele je strukturálně špatně (špatný formát souboru, neplatný parametr...)."""
    status_code = 400
    code = "bad_request"


class InternalError(AppBaseException):
    """
    Generický 500 catch-all pro neočekávané selhání uvnitř zpracování (chyba
    v cizí knihovně, poškozený vstup, který prošel validací, atd.). Vždy nese
    původní zprávu výjimky v `message`, ať je vidět v logu i v odpovědi.
    """
    status_code = 500
    code = "internal_error"


async def app_exception_handler(request: Request, exc: AppBaseException):
    logger.error(f"Chyba aplikace na {request.url.path}: {exc.message}")
    return JSONResponse(
        status_code=exc.status_code,
        content={"error": exc.__class__.__name__, "code": exc.code, "message": exc.message, **exc.payload}
    )
