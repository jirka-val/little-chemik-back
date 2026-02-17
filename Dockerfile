# Použijeme Mambaforge pro stabilní vědecké prostředí
FROM condaforge/mambaforge:latest

WORKDIR /app

# 1. Instalace vědeckých balíčků přes mamba
# Přidán rdkit a libgl1 (nutné pro rdkit/openmm v některých prostředích)
RUN mamba install -y -c conda-forge \
    openmm \
    pdbfixer \
    numpy \
    rdkit \
    && mamba clean -afy

# 2. Instalace zbývajících Python závislostí (FastAPI atd.)
# Ujisti se, že v requirements.txt je 'python-multipart'
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 3. Kopírování zbytku aplikace
COPY . .

EXPOSE 8000

# Spuštění
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers"]