# Použijeme Conda/Mamba obraz, protože openmm a pdbfixer vyžadují specifické binární závislosti
FROM condaforge/mambaforge:latest

WORKDIR /app

# 1. Instalace vědeckých balíčků přes mamba (z conda-forge)
# Tohle zajistí, že OpenMM a PDBFixer budou fungovat správně na Linuxu
RUN mamba install -y -c conda-forge openmm pdbfixer numpy && \
    mamba clean -afy

# 2. Kopírování a instalace Python závislostí pro FastAPI
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 3. Kopírování zdrojového kódu
COPY . .

# Expozice portu 8000
EXPOSE 8000

# Spuštění aplikace pomocí uvicorn
# Bind na 0.0.0.0 je nutný pro přístup zvenčí (veřejná IP)
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers"]