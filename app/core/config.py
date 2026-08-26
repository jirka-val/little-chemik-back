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

    FORCE_FIELDS_CLASSIFICATION_FILE: Path = BASE_DIR / "data" / "force_fields.json"
    FF_CATALOG_SNAPSHOT_FILE: Path = BASE_DIR / "data" / "ff_catalog.json"
    # Jak často (v sekundách) se má katalog FF automaticky obnovovat z IDA na
    # pozadí (viz app/workspaces/tasks/ff_catalog_refresher.py). Výchozí 24h
    # odpovídá Pavlovu "aktualizace jednou denně v noci".
    FF_CATALOG_REFRESH_INTERVAL_SECONDS: int = 24 * 60 * 60
    # Sdílený token pro admin operace (přeřazování FF mezi tiery). Prázdné =
    # endpoint je zamčený úplně, dokud si ho nasazení nenastaví v .env -
    # bezpečnější výchozí stav než "otevřeno pro každého".
    ADMIN_TOKEN: str = ""

    class Config:
        env_file = ".env"
        case_sensitive = True


settings = Settings()