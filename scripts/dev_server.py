"""
Lokální vývojový server s automatickým restartem po změně .py souborů v app/.

    .venv/Scripts/python scripts/dev_server.py [--port 8000]

Proč ne `uvicorn --reload`: na Windows uvicorn při restartu posílá starému
procesu CTRL_C_EVENT a pak na něj čeká bez omezení. Ten signál dorazí jen
procesům se sdílenou konzolí - server spuštěný na pozadí (z IDE, ze skriptu,
z agenta) ho nedostane a restart navždy visí na "Reloading...".
watchfiles starý proces ukončí natvrdo (os.kill -> TerminateProcess), takže
restart funguje vždy. Jen pro vývoj; produkce (Dockerfile) běží bez reloadu.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from watchfiles import PythonFilter, run_process

BACKEND_DIR = Path(__file__).resolve().parents[1]


def _serve(port: int) -> None:
    import os
    import sys

    os.chdir(BACKEND_DIR)
    sys.path.insert(0, str(BACKEND_DIR))
    import uvicorn

    uvicorn.run("app.main:app", host="127.0.0.1", port=port)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    run_process(
        BACKEND_DIR / "app",
        target=_serve,
        args=(args.port,),
        watch_filter=PythonFilter(),
        sigint_timeout=3,
    )


if __name__ == "__main__":
    main()
