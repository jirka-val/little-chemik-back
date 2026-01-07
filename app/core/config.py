from pydantic_settings import BaseSettings
from pathlib import Path
from typing import List


class Settings(BaseSettings):
    PROJECT_NAME: str = "Little Chemik API"
    VERSION: str = "1.2.0"
    API_V1_STR: str = "/api/v1"

    BACKEND_CORS_ORIGINS: List[str] = ["*"]

    BASE_DIR: Path = Path(__file__).resolve().parent.parent.parent
    PDB_DATA_DIR: Path = BASE_DIR / "data" / "pdb_files"

    RCSB_PDB_URL: str = "https://files.rcsb.org/download"

    class Config:
        env_file = ".env"
        case_sensitive = True


settings = Settings()