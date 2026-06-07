from fastapi import APIRouter
from app.api.v1.endpoints import molecules, validation, analysis, download, forcefields, topology, editor

api_router = APIRouter()
api_router.include_router(molecules.router, prefix="/molecules", tags=["molecules"])
api_router.include_router(validation.router, prefix="/validation", tags=["validation"])
api_router.include_router(analysis.router, prefix="/analysis", tags=["analysis"])
api_router.include_router(topology.router, prefix="/topology", tags=["topology"]) # PŘIDÁNO
api_router.include_router(forcefields.router, prefix="/forcefields", tags=["Force Fields"])
api_router.include_router(download.router, prefix="/download", tags=["download"])
api_router.include_router(editor.router, prefix="/editor", tags=["editor"])