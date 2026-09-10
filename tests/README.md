# Testy — Little Chemik backend

## Rozvržení

```
tests/
  conftest.py          sdílené fixtures (TestClient, workspace factory, PDB fixtures,
                        offline_forge_ff pro bezpečné testování FORGE builderu bez sítě)
  fixtures/pdb/         malé, deterministické PDB soubory používané testy
  unit/                 čistá logika, žádné I/O, žádná síť, žádný TestClient
  integration/          přes FastAPI TestClient, lokální workspace lifecycle
  network/               potřebuje skutečné RCSB/IDA API   (marker: network)
  performance/           dlouho běžící zátěžové testy       (marker: slow)
```

**Pravidlo umístění:** pokud test importuje a volá funkci/třídu přímo (žádný
`client.post(...)`), patří do `unit/`. Pokud jde přes `client` fixturu
(FastAPI `TestClient`), patří do `integration/` — pokud navíc nepotřebuje
skutečnou síť. Cokoliv, co stahuje z RCSB nebo IDA API, patří do `network/`.

## Jak spouštět

```bash
# výchozí běh - jen unit + integration, rychlé, žádná síť (~1-2s)
pytest

# konkrétní vrstva
pytest tests/unit -v
pytest tests/integration -v

# zapnout síťové testy (potřebuje internet)
pytest -m network

# zapnout výkonnostní testy
pytest -m slow

# úplně všechno
pytest -m "unit or integration or network or slow"
```

Markery jsou registrované v `pytest.ini` (`--strict-markers` - překlep
v markeru shodí test hned, ne až při čtení výstupu).

## Klíčová bezpečnostní zásada

**Žádný test nesmí zapisovat do `data/ff_cache/` ani `data/ff_cache_forge/`.**
Jsou to skutečná, draze stažená data z IDA API. Testy, které potřebují
silové pole pro FORGE builder (`app/builder`), použijí fixturu
`offline_forge_ff` z `conftest.py` — ta postaví izolovanou kopii v `tmp_path`
z toho, co je už lokálně nacachované, a pokud daný FF lokálně chybí, test se
elegantně přeskočí (`pytest.skip`) místo pádu nebo (horšího) přepsání reálné
cache prázdným/špatným obsahem.

Ze stejného důvodu `unit/test_forcefield_service.py` vždy monkeypatchuje
`ForceFieldService.CACHE_DIR`/`FORGE_CACHE_DIR` na `tmp_path` přes fixturu
`service` — nikdy nepracuje s opravdovými adresáři.

## Fixture PDB soubory (`fixtures/pdb/`)

| soubor | co testuje |
|---|---|
| `alanine_single.pdb` | nejmenší platný protein fragment |
| `protein_gap.pdb` | GLU83 → [84-88 chybí] → PHE89, syntetická obdoba 1JJ2 chain K - testuje terminalitu na okraji sekvenční díry |
| `rna_gap.pdb` | **reálný výřez ze struktury 1RNA** (chain A, residua 1-3 a 6-8, 4-5 vynechána) - musí mít skutečnou geometrii, ne vymyšlené souřadnice, protože FORGE builder na degenerovaných/kolineárních atomech spadne s `Cannot normalize near-zero central dihedral bond` |
| `altloc_sample.pdb` | SER42 se dvěma alternativními konformacemi OG (occupancy 0.6/0.4) |
| `two_model_nmr.pdb` | stejné reziduum ve 2 MODEL blocích (NMR ensemble) |

Když přidáváš nový fixture PDB, který půjde přes `ForgeStructureService`
(tzn. cokoliv v `integration/test_prepare_*` nebo `network/`), musí mít
reálnou, ne-degenerovanou geometrii - buď ručně dopočítanou, nebo (spolehlivěji)
vyříznutou ze skutečné stažené struktury, jako `rna_gap.pdb`.

## Solvatace a periodický box

`tests/network/test_solvation_box.py` je obdoba starého
`test_performance.py::test_03_solvation_performance` (PDBFixer `addSolvent`),
teď nad `ForgeStructureService`. Klíčový poznatek zjištěný při psaní: builder
hledá vodní silové pole striktně pod `mol_type="W3"` (ne obecné `"W"`, které
používá starší `TopologyService`/`pdb_service` pipeline pro AMBER topologii),
a i prázdný seznam solí spustí síťovou neutralizaci výchozími K+/Cl- ionty
(`mol_type="I1"`) - `ff_selections` proto pro solvataci vždy potřebuje i `I1`
položku, jinak spadne na `KeyError: Ion parameters missing for I1:K+`.

`TestSolvationCreatesBox` na 1RNA běží reálně a dokončí se (máme lokálně
OL3 + TIP3P + JC-TIP3P-I1). `TestLargeStructureSolvationPerformance` na 1JJ2
je parita se starým testem, ale **1JJ2 solvataci nikdy nedokončí** - narazí
na legitimní `missing_dof` přesně na `K:83 (CGLU:CG)` (GLU83, stejné reziduum
z úvodní diskuze o 1JJ2 gapu) po ~93 s, dřív než solvatace vůbec začne. To je
očekávané, správné chování (1JJ2 je reálně neúplná struktura), ne bug - test
proto měří čas do tohoto bodu, ne čas do úspěšné solvatace. Stejný nález a
stejný čas potvrzuje i `tests/performance/test_forge_performance.py::TestForgeBuildPerformance`.

## Známé, zdokumentované mezery (ne bugy v tomhle kódu)

- `tests/network/test_reference_structures.py::TestDnaReference::test_1bna_builds_successfully`
  je `xfail` - lokálně nacachované `FF99BSC0` nemá parametr pro `DT:H72`. Mezera
  v datech konkrétního FF, ne v kódu `ForgeStructureService`. Smaž `xfail`, až
  bude FF opravené/doplněné.
- `tests/unit/test_pdb_topology_dict.py::TestWithoutSidecar::test_protein_defaults_to_mol_type_r`
  dokumentuje pre-existující (mimo rozsah FORGE integrace) nedostatek:
  `parse_pdb_to_topology_dict` bez `forge_meta` sidecaru přiřadí proteinům
  `mol_type="R"` místo `"P"`.

## Přidávání nových testů

1. Vyber vrstvu podle pravidla výše.
2. Pokud test potřebuje FORGE builder (`ForgeStructureService`/`run_forge_workflow`),
   vždy použij `offline_forge_ff` fixturu - nikdy nevolej `ForceFieldService`
   napřímo bez monkeypatche v testu, který běží proti reálné `data/`.
3. Pokud test potřebuje PDB s reálnou geometrií, buď ho vyřízni ze skutečně
   stažené struktury (viz `rna_gap.pdb` výše), nebo ověř, že tvůj syntetický
   vstup neprochází přes builder (čistě tokenová/textová logika v `analysis_service`
   syntetické souřadnice snese bez problémů).
4. Nový marker přidávej i do `pytest.ini` (`--strict-markers` ho jinak odmítne).
