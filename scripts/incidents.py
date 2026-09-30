"""
Výpis a stažení nahlášených chyb z běžícího FORGE serveru (admin API
/api/incidents, viz app/api/v1/endpoints/incidents.py).

Heslo se bere z proměnné prostředí FORGE_ADMIN_TOKEN, jinak z ADMIN_TOKEN v
little-chemik-back/.env (stejná hodnota jako ADMIN_TOKEN na serveru); adresa
z FORGE_URL (výchozí http://147.251.115.223).

    python scripts/incidents.py list
    python scripts/incidents.py download 20261001-101500-abc123
    python scripts/incidents.py download --all --out problem-pdb/incidents

Stažený report se rovnou rozbalí do složky <out>/<id>/. Jen standardní knihovna.
"""

import argparse
import io
import json
import os
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

DEFAULT_URL = "http://147.251.115.223"


def _request(base_url: str, token: str, path: str) -> bytes:
    request = urllib.request.Request(f"{base_url.rstrip('/')}{path}", headers={"X-Admin-Token": token})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return response.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 403:
            sys.exit("403 - wrong or missing token (FORGE_ADMIN_TOKEN must equal ADMIN_TOKEN on the server).")
        if exc.code == 404:
            sys.exit(f"404 - not found: {path}")
        sys.exit(f"HTTP {exc.code}: {exc.read()[:300]!r}")
    except urllib.error.URLError as exc:
        sys.exit(f"Cannot reach {base_url}: {exc.reason}")


def list_incidents(base_url: str, token: str) -> list:
    data = json.loads(_request(base_url, token, "/api/incidents"))
    incidents = data["incidents"]
    print(f"{len(incidents)} report(s), {data['storage_used_mb']} MB of {data['storage_limit_mb']} MB used\n")
    for item in incidents:
        text = (item.get("description") or item.get("error") or "").replace("\n", " ")
        print(f"{item['id']}  {item['kind']:<6}  {item['size_mb']:>7.2f} MB  {text[:90]}")
    return incidents


def download(base_url: str, token: str, incident_id: str, out_dir: Path) -> Path:
    data = _request(base_url, token, f"/api/incidents/{incident_id}/download")
    out_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        archive.extractall(out_dir)
    target = out_dir / incident_id
    print(f"{incident_id} -> {target}")
    return target


def _token_from_env_file() -> str:
    env_file = Path(__file__).resolve().parent.parent / ".env"
    try:
        lines = env_file.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ""
    for line in lines:
        key, _, value = line.partition("=")
        if key.strip() == "ADMIN_TOKEN":
            return value.strip().strip('"').strip("'")
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description="List / download FORGE error reports.")
    parser.add_argument("--url", default=os.environ.get("FORGE_URL", DEFAULT_URL))
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list all reports")
    dl = sub.add_parser("download", help="download and unpack reports")
    dl.add_argument("ids", nargs="*", help="report IDs")
    dl.add_argument("--all", action="store_true", help="download every report")
    dl.add_argument("--out", default="incidents", help="target directory (default ./incidents)")
    args = parser.parse_args()

    token = os.environ.get("FORGE_ADMIN_TOKEN", "") or _token_from_env_file()
    if not token:
        sys.exit("Set FORGE_ADMIN_TOKEN (or ADMIN_TOKEN in little-chemik-back/.env) to the server's ADMIN_TOKEN.")

    if args.command == "list":
        list_incidents(args.url, token)
        return

    ids = list(args.ids)
    if args.all:
        ids = [item["id"] for item in json.loads(_request(args.url, token, "/api/incidents"))["incidents"]]
    if not ids:
        sys.exit("Give report IDs or --all.")
    for incident_id in ids:
        download(args.url, token, incident_id, Path(args.out))


if __name__ == "__main__":
    main()
