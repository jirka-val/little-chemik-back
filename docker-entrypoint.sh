#!/bin/sh
# Start kontejneru: datové složky (i volumes z dřívějška, vytvořené rootem)
# předá uživateli "forge" a server spustí pod ním, ne jako root.
set -e
if [ "$(id -u)" = "0" ]; then
    mkdir -p /app/temp_workspaces /app/data/incidents /app/data/ff_cache /app/data/ff_cache_forge
    chown -R forge:forge /app/temp_workspaces /app/data
    exec setpriv --reuid=forge --regid=forge --init-groups "$@"
fi
exec "$@"
