# Plán optimalizace: víc uživatelů najednou bez zpomalení

Cíl: stránka má zůstat rychlá, i když ji používá víc lidí současně. Kolik jich
bude, zatím nevíme, takže server musí **škálovat s počtem jader** a **férově
řadit** práci, místo aby se všichni zpomalili stejně.

Každá změna se měří stejným benchmarkem a porovnává s baseline:

```bash
.venv/Scripts/python.exe -m tests.load.benchmark run --label <změna>
.venv/Scripts/python.exe -m tests.load.benchmark compare tests/load/benchmarks/2026-11-09_1257_baseline_75a93ea.json tests/load/benchmarks/<nový>.json
```

Popis benchmarku a zátěžových testů: `tests/load/README.md`.

## Výchozí stav (baseline 2026-11-09)

Vývojové PC (24 jader), 1 worker uvicornu jako v Dockerfile, střední systém =
1CA2 (~30 tis. atomů). Plná data: `tests/load/benchmarks/2026-11-09_1257_baseline_75a93ea.md`.

| KPI | baseline |
|---|---|
| 1 uživatel, střední systém | 5,5 s (5,1 s CPU) |
| 1 uživatel, velký systém (124 tis. atomů) | 13,8 s |
| 4 / 8 / 16 uživatelů najednou, nejpomalejší čeká | 25 / 52 / 103 s |
| propustnost | ~9–11 příprav/min, **ať je jader kolik chce** |
| health check při 8 přípravách, p95 / max | 0,20 / ~2 s |
| uživatel každé 2 s, p95 zpomalení | 8,1× |
| RAM v klidu / při 16 uživatelích | 190 / 680 MB |

### Proč to je pomalé

1. **Server využije jen jedno jádro.** Jeden proces uvicornu a výpočty v Pythonu
   (GIL): přidání uživatele lineárně prodlouží čekání všem.
2. **Nefér řazení.** Souběžné přípravy se v threadpoolu střídají, takže při 16
   uživatelích skončí všichni až za ~100 s (p50 ≈ p95). Při frontě by první
   měl výsledek za 6 s a průměr by byl zhruba poloviční, bez jakéhokoli
   zrychlení výpočtu.
3. **Event loop občas čeká na GIL.** Lehké požadavky (health, sekvence) mají
   při zátěži špičky až 2 s, protože výpočetní vlákna drží GIL.
4. **Zbytečná práce v každé session** (profil 1CA2, podíl na CPU):

   | kde | podíl | poznámka |
   |---|---|---|
   | builder: plánování kroků (`plan_build_steps`) | ~27 % | 783 tis. vyhodnocení `resolve_ref`, `build_planner_dependency_index` |
   | builder: MM hledání volných rotorů | ~20 % | 192 rotorů, 2 144 zkušebních konformací |
   | builder: solvatace + ionty | ~9 % | |
   | validace PDB (`validate_pdb_content`) | ~13 % | 2× za session, pokaždé PDBFixer (načtení PDB přes OpenMM) |
   | topologie | ~13 % | z toho zápis prmtop ~4 %, parsování FF souborů při každém volání |
   | sekvence (`build_sequence_tokens`) | ~7 % | 8× za session, pokaždé znovu parsuje PDB |
   | serializace JSON odpovědí | ~7 % | 409 tis. objektů, **v event loopu** |

## Fáze 1: rychlé výhry bez změny architektury (1–2 dny)

| # | změna | očekávaný dopad | KPI |
|---|---|---|---|
| 1.1 | **Fronta těžkých úloh**: semafor na příprava / topologie / export (počet = počet worker procesů × 1), FIFO. Volitelně endpoint s pozicí ve frontě pro UI („jste 3. na řadě“). | průměrné čekání při souběhu ~½, první uživatelé hned; méně soupeření o GIL → lepší odezva lehkých požadavků | concurrency p50, arrivals slowdown p50, health p95 |
| 1.2 | **Validace jen jednou**: výsledek `validate_pdb_content` cachovat podle obsahu PDB (hash), PDBFixer nevytvářet, když se PDB nezměnilo | −10 až −13 % CPU na session | single.*.server_cpu_s |
| 1.3 | **Cache sekvence**: tokeny podle (workspace, mtime souboru) | −5 % CPU, rychlejší přepínání panelů | single.*.server_cpu_s, light_burst |
| 1.4 | **JSON mimo event loop**: velké odpovědi (příprava, sekvence, FF) serializovat přes `orjson` / `ORJSONResponse`, případně zmenšit payload | −5 % CPU, menší špičky lehkých požadavků | responsiveness |
| 1.5 | **gzip odpovědí** (nginx `gzip on` pro JSON/PDB/prmtop) | rychlejší přenos u pomalejšího připojení; lokálně se neprojeví | velikost odpovědí (doplnit do benchmarku) |
| 1.6 | **Zúžit odpověď silových polí** (nově zjištěno): `GET /api/forcefields` posílá celé obsahy FF souborů (7,1 MB) a frontend je posílá zpátky v `ff_selections` při přípravě i topologii. Frontend by měl posílat jen identifikátor FF, backend obsah vezme z katalogu (zároveň bezpečnější - klient teď může podstrčit libovolný obsah FF). | −7 MB stahování a −1–2 MB uploadu na session | velikost požadavků/odpovědí |

### Výsledky fáze 1 (2026-11-09, `2026-11-09_1329_phase1_feff387.json`)

Hotovo: 1.1 fronta (pořadí podle session), 1.2 validace bez vody v PDBFixeru,
1.3 cache sekvence, 1.4 velké JSON odpovědi mimo event loop, 1.5 gzip v nginx.

| KPI | baseline | fáze 1 | změna |
|---|---|---|---|
| 1 uživatel, střední / velký systém | 5,5 / 13,8 s | 5,3 / 12,4 s | −4 / −10 % |
| 1 uživatel, velký systém: RAM navíc | 299 MB | 168 MB | −44 % |
| 8 najednou: medián / nejpomalejší | 52 / 52 s | 37 / 45 s | −29 / −14 % |
| 16 najednou: medián / nejpomalejší | 103 / 103 s | 62 / 93 s | −40 / −10 % |
| 16 najednou: špička RAM | 684 MB | 493 MB | −28 % |
| health / sekvence při zátěži, p95 | 0,20 / 0,92 s | 0,14 / 0,65 s | −32 / −29 % |
| uživatel každé 2 s: medián session | 13,9 s | 22,4 s | **+62 %** |
| uživatel každé 2 s: p95 zpomalení | 8,1× | 19,7× | **+143 %** |

Gzip se lokálně neprojeví; spočítáno: silová pole 7,1 → 0,45 MB, připravené
PDB 2,8 → 0,4 MB, sekvence 3,6 → 0,07 MB (na 20 Mbit/s 3,0 s → 0,2 s u FF).

Poučení:
- Fronta po jednotlivých požadavcích nepomohla (session má víc těžkých
  požadavků, střídaly se dál). Pořadí podle session zkrátilo medián o 40 %.
- Za přetížení (uživatel každé 2 s je víc, než jedno jádro stihne) fronta
  znevýhodní malé úlohy: čekají za celými velkými. Dřív se o CPU dělily
  rovnoměrně. Na jednom jádře je to jen přerozdělení čekání - skutečné řešení
  je víc jader (fáze 2). Případně řadit malé úlohy přednostně (podle velikosti
  PDB), pokud to bude v provozu vadit.
- Samotné zrychlení výpočtu je zatím malé (−4 až −10 %), velký zbytek je v
  builderu (fáze 3).

### Výsledek bodu 1.6 (2026-11-09, `2026-11-09_1407_1.6-ff-by-id_00cf797.json`)

FF endpoint posílá FF bez obsahů souborů (`app/services/ff_selection.py`,
`public_ff`); příprava a topologie berou obsah z katalogu podle `id`
(`resolve_ff_selections`), FF mimo katalog odmítnou (400 `unknown_force_field`).
Hlavní soubor FF pro tlačítko Download: `GET /api/forcefields/catalog/{id}/main-file`.
Frontend beze změny formátu - posílá dál to, co dostal (jen už bez souborů).

| | před | po |
|---|---|---|
| seznam FF (stahování) | 7,1 MB | 61 kB |
| příprava + topologie (nahrávání) | 3,4 MB | 14 kB |
| špička RAM, 16 uživatelů | 493 MB | 384 MB |
| health / sekvence při zátěži, p95 | 0,14 / 0,65 s | 0,095 / 0,45 s |
| uživatel každé 2 s: p95 zpomalení | 19,7× | 15,7× |

Zároveň opravena bezpečnostní díra: obsah FF i `display_name` (z něj se
skládá cesta ke složce FF cache) dřív přicházely od klienta.

## Fáze 2: víc jader (2–4 dny)

| # | změna | očekávaný dopad | KPI |
|---|---|---|---|
| 2.1 | **Víc worker procesů** (`uvicorn --workers N`, N ≈ počet vCPU). Změřeno se 4 workery: propustnost 10 → 25 příprav/min, 16 uživatelů 95 → 38 s, health p95 0,05 s. | propustnost ~lineárně s jádry | concurrency.16.*, throughput |
| 2.2 | **Podmínka pro 2.1:** stav dostavování bočních řetězců (`sidechain_session_service`, start → update/optimize → commit) je v paměti procesu. Buď ho ukládat do workspace (serializace na disk), nebo v nginx směrovat podle workspace id na stejný worker (`hash $workspace_id consistent`). | víc workerů bez rozbití GUI dostavování | nový zátěžový test se systémem, který vyvolá `missing_dof` |
| 2.3 | Do benchmarku a zátěžových testů přidat systém s `missing_dof` (start → commit) | pokrytí cesty, která se s víc workery rozbije | — |

Alternativa k 2.1: process pool pro těžké kroky uvnitř jednoho serveru. Je
složitější kvůli stavu bočních řetězců (objekty builderu by musely žít v
jiném procesu), proto až druhá volba.

### Jak se fáze 2 nakonec udělala: pool výpočetních procesů (2026-11-09)

Místo 2.1 + 2.2 (víc procesů uvicornu + ukládání relací nebo směrování v
nginx) jeden webový proces a **pool výpočetních procesů**
(`app/core/process_pool.py`):

- Do poolu jdou bezstavové těžké výpočty: příprava (`/sidechains/start`),
  topologie, analýza solvatace (+ check-box), PDB z topologie a souřadnic.
- Interaktivní dostavování bočních řetězců zůstává ve webovém procesu: když
  výpočet v poolu narazí na `missing_dof`, vrátí `needs_session` a `/start`
  ho zopakuje ve webovém procesu, kde relace mezi požadavky přežije
  (vzácný případ, platí se jen tam). Výsledek přípravy se přenáší pickle
  (2,5 MB / 2 ms u středního systému, 10 MB / 10 ms u velkého).
- Proč ne víc procesů uvicornu: relace by se musely ukládat na disk při
  každém pohybu slideru, nebo by se muselo směrovat v nginx (víc instancí
  serveru, jiná konfigurace v Dockeru i lokálně). Pool funguje všude stejně
  a navíc drží výpočty mimo webový proces (žádné čekání event loopu na GIL).
- `HEAVY_JOB_WORKERS`: -1 = automaticky (jádra - 1, nejvýš 4), 0 = bez
  poolu (testy přes TestClient), N = N procesů. Fronta (`job_queue.py`)
  pouští najednou tolik úloh, kolik je procesů.
- Zprávy konzole z procesů poolu se přeposílají do webového procesu (konzole
  ve frontendu je dál vidí). Výjimky s HTTP kódem projdou beze změny.
- Procesy se při startu zahřejí (importy builderu), při vypnutí ukončí a
  když webový proces spadne natvrdo, ukončí se samy (hlídač rodiče).
- Katalog FF se v každém procesu načte znovu, když se změní soubor
  (refresh v jiném procesu).
- 2.3 hotovo: zátěžový test `test_interactive_sidechain_sessions_under_load`
  (6 uživatelů najednou start → commit → topologie, s poolem).

Ověřeno: golden testy s poolem (`HEAVY_JOB_WORKERS=2 pytest -m golden`)
dávají bajtově stejné výsledky jako bez něj; zátěžové testy 17 passed.

### Výsledky fáze 2 (`2026-11-09_1511_phase2-pool4_635f99b.json`, 4 procesy)

| KPI | baseline | fáze 1 + 1.6 | fáze 2 |
|---|---|---|---|
| 4 najednou: nejpomalejší | 25 s | 22 s | **8,8 s** |
| 8 najednou: medián / nejpomalejší | 52 / 52 s | 38 / 46 s | **12,8 / 16,2 s** |
| 16 najednou: medián / nejpomalejší | 103 / 103 s | 62 / 93 s | **23,9 / 31,1 s** |
| propustnost (16 najednou) | 9,3 /min | 10,4 /min | **29,7 /min** |
| sekvence při zátěži, p95 | 0,92 s | 0,45 s | **0,09 s** |
| uživatel každé 2 s: medián session / p95 zpomalení | 13,9 s / 8,1× | 19,6 s / 15,7× | **2,2 s / 1,5×** |
| RAM v klidu / špička při 16 | 188 / 684 MB | 188 / 384 MB | **753 / 1413 MB** |

Cena: každý proces poolu má v klidu ~140 MB, při práci podle velikosti
systému. Se 4 procesy ~0,75 GB v klidu, ~1,4 GB při plném vytížení
středními systémy - VM s 8 GB RAM má rezervu, velké systémy (300 tis.
atomů, ~0,8 GB na úlohu) ve 4 procesech najednou ~4 GB.

Přetížení z fáze 1 zmizelo: při uživateli každé 2 s server stíhá (medián
2,2 s, malé úlohy už nečekají za velkými). Jeden uživatel je pořád ~5 s -
to je výpočet sám, fáze 3.

## Fáze 3: zrychlení samotného výpočtu (průběžně)

| # | změna | očekávaný dopad |
|---|---|---|
| 3.1 | **Plánování builderu**: výsledky vyhodnocení pravidel (`resolve_ref`, `find_matching_torsion_group`) jsou pro stejný typ rezidua stejné → memoizace podle šablony rezidua | odhad −15 až −20 % CPU přípravy |
| 3.2 | **MM hledání volných rotorů**: vektorizovat skórování zkušebních konformací, hrubší první krok mřížky, předčasný konec | odhad −10 % |
| 3.3 | **Topologie**: parsované silové pole (FF_IDA) cachovat v paměti podle verze FF; rychlejší zápis prmtop (`write_int`/`write_dbl` po blocích) | −5 až −8 % |
| 3.4 | Velké systémy: profil pro 124–300 tis. atomů (solvatace, ionty, zápis) | velké systémy |

Builder (`app/builder/`) a `AMBER_topology.py` jsou citlivé: každá změna musí
projít golden testy (`pytest -m golden`), výstupy musí zůstat bajtově stejné.

## Fáze 4: když nebude stačit jeden server

- Fronta úloh mimo webový proces (Redis + worker kontejnery, např. arq/RQ),
  stav úlohy a průběh se dotazuje z frontendu. Workery se pak dají přidávat
  na další stroje.
- Limit souběžných úloh na uživatele / IP, aby jeden uživatel nezahltil frontu.
- Monitoring v produkci: doba odezvy po endpointech, délka fronty, CPU, RAM,
  disk, výstraha při p95 nad limitem.

## Cílové hodnoty

Na produkčním VM (4 vCPU, 4 workery) po fázích 1 a 2:

| KPI | baseline | cíl |
|---|---|---|
| 1 uživatel, střední systém | 5,5 s | ≤ 4 s |
| 8 uživatelů najednou, průměrné čekání | ~50 s | ≤ 15 s |
| propustnost | ~10/min | ≥ 35/min |
| health check při zátěži, p95 | 0,2 s | ≤ 0,05 s |
| lehké požadavky (sekvence) při zátěži, p95 | 0,9 s | ≤ 0,2 s |

Benchmark se musí zopakovat i na produkčním VM, jakmile bude: vývojové PC
má jiný procesor a čísla nejsou přímo přenosná (nejstabilnější je CPU čas na
přípravu).

## Známé problémy nalezené při měření (2026-11-09)

- Opraveno: přepisování sdílené FF cache při souběhu (chyby 500), generování
  topologie blokovalo celý server, nahrání prmtop + rst7 nad 9999 vod.
- Otevřené: rekonstrukce boxu u velkých RNA boxů bez CRYST1 (xfail v
  `tests/load/test_capacity.py`).
