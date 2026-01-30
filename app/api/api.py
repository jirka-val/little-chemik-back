from fastapi import APIRouter
from app.api.v1.endpoints import molecules, validation, analysis

api_router = APIRouter()

api_router.include_router(molecules.router, prefix="/v1/molecules", tags=["molecules"])
api_router.include_router(validation.router, prefix="/v1/validate", tags=["validation"])
api_router.include_router(analysis.router, prefix="/v1/analysis", tags=["analysis"])