# 1. Použijeme moderní Miniforge3 (oficiální nástupce Mambaforge)
FROM condaforge/miniforge3:latest

WORKDIR /app

# 2. Instalace vědeckých balíčků přes mambu
# V novém obrazu mamba neudělá "sebevraždu" a korektně vše nainstaluje
RUN mamba install -y -c conda-forge \
    openmm \
    pdbfixer \
    numpy \
    rdkit \
    && mamba clean -afy

# 3. Instalace zbývajících Python závislostí
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Kopírování zbytku aplikace
COPY . .

EXPOSE 8000

# 5. Spuštění
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers"]