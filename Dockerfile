# Python a všechny balíčky ve stejných verzích, se kterými prochází testy
# (requirements.txt + constraints.txt). Dřív miniforge3:latest + mamba bez
# verzí: každé sestavení mohlo tiše natáhnout jinou verzi OpenMM/numpy.
FROM python:3.12.10-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

COPY requirements.txt constraints.txt ./
RUN pip install -r requirements.txt -c constraints.txt

COPY . .

# Server běží jako neprivilegovaný uživatel "forge"; docker-entrypoint.sh
# při startu (jako root) předá datové složky a volumes tomuto uživateli.
RUN useradd --system --uid 10001 --user-group --no-create-home forge \
    && chmod +x /app/docker-entrypoint.sh

EXPOSE 8000

# Zdraví kontejneru (docker ps / compose depends_on): odpovídá /api/health?
HEALTHCHECK --interval=30s --timeout=5s --start-period=90s --retries=3 \
    CMD python -c "import sys, urllib.request; sys.exit(urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=4).status != 200)"

ENTRYPOINT ["/app/docker-entrypoint.sh"]
# --forwarded-allow-ips "*": backend není publikovaný ven (jen přes nginx
# kontejner), takže X-Forwarded-Proto/For od proxy může brát jako pravdivé -
# za HTTPS pak ví, že běží na https:// (bez toho uvicorn věří jen 127.0.0.1).
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers", "--forwarded-allow-ips", "*"]
