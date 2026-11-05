# Stima differenziale dell'anisotropia di conduzione atriale

Tesi triennale di geometry processing. Su un atrio sinistro con fibre e
velocita' di conduzione (CV) note si simulano sei mappe di attivazione (LAT);
da queste si stima per faccia il tensore di conduzione (direzione di fibra,
CV_L, CV_T) con il gradiente per faccia, e si valida la stima contro la
verita' nota, anche a confronto con il metodo 2 di Roney et al. 2019
(*Computers in Biology and Medicine* 104, 278-290).

```
run_all.sh                     ricostruisce data/reference/ e output/ da data/
environment.yml                ambiente conda "tesi"
data/
  external/biatrial_model_completeFit.h5      SSM (incluso)
  augmenta/LA_endo_with_fiber.{pts,elem,lon}  LA con fibre, uscita di AugmentA (inclusa)
simulator/
  stage1/                      LA dall'SSM (facoltativo, vedi sotto), siti di stimolo
    opencarp/                  mesh in um, tag per zona, simulazioni eikonali
  convert/                     da openCARP a PLY: mappe LAT e verita' per faccia
estimator/
  src/                         Tesi0 (C++): gradiente per faccia, stima del tensore
  include/happly.h             lettore PLY vendorizzato (THIRD_PARTY.md)
  roney.py                     metodo 2 di Roney et al. 2019
analysis/analisi.py            errori pesati per area, figure, numeri.csv e confronto.csv
experiments/                   modello.conf e roney.conf, CV per zona (region_cv/),
                               subsets.py, mappe.py; formati in experiments/README.md
```

## Requisiti

- Linux; sviluppato su Ubuntu 26.04.1.
- conda, con l'ambiente di `environment.yml` (Python 3.13.15, numpy, scipy,
  pandas, matplotlib, h5py, meshio, pyvista 0.48.4, vtk 9.6.2: versioni
  esatte nel file).
- openCARP v19.0, dal pacchetto `.deb` standard (non ubuntu20):
  <https://opencarp.org/download/installation>. Serve `/usr/local/bin/openCARP`;
  `carputils` non serve.
- CMake >= 3.20 (testato 4.2.3), un compilatore C++20 (testato GCC 15.2.0) e
  glm, header-only (su Ubuntu `libglm-dev`, testato 1.0.1).
- Docker, solo per rigenerare facoltativamente l'uscita di AugmentA.

## Dati inclusi

Sono gli unici ingressi che non si rigenerano; tutto il resto di `data/` e
`output/` lo scrive `run_all.sh`.

**SSM** `data/external/biatrial_model_completeFit.h5` (23.9 MB). Modello
statistico di forma biatriale del KIT-IBT, Zenodo
<https://zenodo.org/records/4309958>, licenza CC-BY 4.0.
md5 `96d8bf64bc943dcba71e3702a552fa75`,
sha256 `d4b3a90f2b9c6db79bb225954a4e5eeaed59f44916e5a93935dde20458da63b8`.
Citazione: C. Nagel, S. Schuler, O. Dössel, A. Loewe, *A bi-atrial
statistical shape model for large-scale in silico studies of human atria:
model development and application to ECG simulations*, Medical Image
Analysis 74 (2021) 102210.

**Uscita di AugmentA** `data/augmenta/LA_endo_with_fiber.{pts,elem,lon}`
(3.7 MB): l'atrio sinistro medio dell'SSM, 21811 vertici e 43270
triangoli, con fibre e tag di regione per triangolo. Generata dall'SSM con
`make_mean.py`, `extract_la.py` e poi AugmentA con i flag
della sezione *Rigenerare la LA e le fibre*. Il commit di AugmentA usato
non e' noto. Licenza di AugmentA: Academic Public License. Citazione:
L. Azzolin et al., *AugmentA: patient-specific augmented atrial model
generation tool*, Computerized Medical Imaging and Graphics 108 (2023)
102265.

I dati clinici non sono nel repository.

## Riprodurre i risultati

```bash
git clone https://github.com/n0sum/cardiofiber
cd cardiofiber
conda env create -f environment.yml
conda activate tesi
```

Installare openCARP v19.0 (vedi *Requisiti*), poi dalla radice, con
l'ambiente `tesi` attivo:

```bash
./run_all.sh
```

Lo script converte la mesh in um, scrive i tag degli scenari, lancia le
simulazioni eikonali (sei siti di stimolo, `simulator/stage1/pacing_sites.csv`),
compila `estimator/` in `estimator/build/`, cancella `output/` e calcola
una sola griglia di esperimenti: per ogni scenario lo stimatore per faccia su
tutti i 42 sottoinsiemi di almeno 3 mappe e il metodo di Roney sulle 20
terne. Uscite in `output/`:

| cartella                                   | contenuto                                          |
|--------------------------------------------|----------------------------------------------------|
| `aniso_zone_N/`                            | `confs/sub_*.conf`; `faccia/sub_*/errors.csv` (Tesi0); `roney/sub_i-j-k/errors.csv` (`roney.py`) |
| `figure/`                                  | `numeri.csv`, `confronto.csv` e 16 grafici (`analisi.py`) |
| `mesh/`                                    | 12 PLY degli errori, uno per errore e caso, per MeshLab (`mappe.py`) |

Scenari: `aniso_zone_N`, la mesh divisa in N zone casuali (seme 0), ogni zona
con una delle 25 combinazioni di CV_L e r = CV_L/CV_T di
`region_cv/aniso_random.csv`, N = 1, 8, 70, 600: diametro medio di zona
d = sqrt(4 A / (pi N)) ~ 123 / sqrt(N) mm, cioe' 123, 44, 15 e 5.0 mm, passi
uguali in scala log. Con N = 1 le CV sono uniformi (e' uscito il tag 16:
CV_L 1.155, CV_T 0.722 mm/ms, r = 1.6); con N = 600 le zone sono grandi come
il disco del metodo di Roney (diametro 5 mm), che presuppone CV uniformi al
suo interno: con zone piu' piccole il confronto misurerebbe solo che il
presupposto e' violato. Dettagli in `experiments/README.md`.

## Risultati attesi

Errori pesati per area sulle facce stimabili: angolo di fibra in gradi,
errori relativi in valore assoluto, in %, su CV_L e sul rapporto r.

Le tabelle sono in `output/figure/`: `numeri.csv` per lo stimatore per faccia
(una riga per scenario e sottoinsieme), `confronto.csv` per il confronto con
Roney (una riga per scenario, terna e metodo, piu' le righe `sub` = `tutte`).
Le colonne hanno anche q99, che qui non e' riportato.

**aniso_zone_1**, sei mappe (`numeri.csv`, riga `sub_0-1-2-3-4-5`):

| area non stimabile | angolo q50 | angolo q90 | CV_L q50 | CV_L q90 | r q50 | r q90 |
|--------|--------|-------|--------|--------|--------|---------|
| 0.13 % | 0.907° | 4.27° | 0.67 % | 3.32 % | 1.30 % | 6.90 %  |

**Sottoinsiemi di mappe** su `aniso_zone_1` (`numeri.csv`, errori sulle
facce stimabili di ogni sottoinsieme): con 3 mappe, sulle 20 terne, area non
stimabile 4.0-9.7 % (mediana 5.3 %), angolo q50 1.17-1.77° (mediana
1.37°), angolo q90 9.2-15.0° (mediana 10.6°); con 6 mappe i valori di
`aniso_zone_1`.

**Zone**, sei mappe (`numeri.csv`, righe `sub_0-1-2-3-4-5`), d = diametro medio di zona:

| N   | d [mm] | area non stimabile | angolo q50 | angolo q90 | CV_L q50 | r q50  |
|-----|--------|--------------------|------------|------------|----------|--------|
| 600 | 5.03   | 1.28 %             | 1.78°      | 9.50°      | 1.99 %   | 4.69 % |
| 70  | 14.7   | 0.57 %             | 1.22°      | 6.50°      | 1.34 %   | 2.99 % |
| 8   | 43.6   | 0.42 %             | 0.882°     | 4.70°      | 1.00 %   | 2.16 % |
| 1   | 123    | 0.13 %             | 0.907°     | 4.27°      | 0.67 %   | 1.30 % |

**Roney contro stimatore per faccia**, terne di mappe (`confronto.csv`, righe
`sub` = `tutte`): area non stimabile di ogni metodo da solo, errori sulle
facce stimabili in entrambi (area appaiata), messe insieme sulle 20 terne:

| N   | area appaiata | metodo | area non stimabile | angolo q50 | angolo q90 | CV_L q50 | r q50  |
|-----|---------------|--------|--------------------|------------|------------|----------|--------|
| 1   | 75.7 %        | Roney  | 22.5 %             | 6.65°      | 21.2°      | 12.9 %   | 10.9 % |
|     |               | faccia | 5.8 %              | 1.31°      | 9.09°      | 0.82 %   | 2.13 % |
| 8   | 55.0 %        | Roney  | 42.4 %             | 9.16°      | 29.1°      | 13.2 %   | 14.6 % |
|     |               | faccia | 9.8 %              | 1.34°      | 10.0°      | 1.08 %   | 3.24 % |
| 70  | 54.8 %        | Roney  | 41.8 %             | 8.75°      | 30.2°      | 13.9 %   | 15.4 % |
|     |               | faccia | 12.5 %             | 1.63°      | 12.4°      | 1.27 %   | 3.84 % |
| 600 | 52.3 %        | Roney  | 42.6 %             | 8.69°      | 29.5°      | 15.0 %   | 17.8 % |
|     |               | faccia | 16.2 %             | 2.38°      | 17.0°      | 1.90 %   | 6.41 % |

## Unita' e formati

- Mesh e SSM in **mm**; openCARP vuole **um**: `mm_to_um.py` moltiplica per
  1000 le sole coordinate. LAT in ms, quindi CV in mm/ms. Gli errori di unita'
  sono silenziosi; controllo: l'estensione della LA e' circa [48, 71, 83] mm.
- `.pts`: prima riga il numero di vertici, poi `x y z` per riga.
- `.elem`: prima riga il numero di triangoli, poi `Tr i j k tag`, indici da 0,
  tag di regione 1-7 di AugmentA (`simulator/stage1/region_tags.txt`).
  AugmentA inverte il verso dei triangoli rispetto all'SSM; la connettivita'
  di riferimento e' sempre `LA_endo_with_fiber.elem`.
- `.lon`: prima riga `2`, poi per triangolo fibra e sheet, tre componenti ciascuna.
- `experiments/region_cv/aniso_random.csv`: `region_tag,CV_L,CV_T` in mm/ms;
  `run_sims.py` converte in mm/s per openCARP.
- `data/reference/LA_LAT_NOME.ply`: PLY binario, per vertice `x y z` (mm) e
  `lat_0 ... lat_5` (ms), `lat_k` dal sito k. `LA_GT_NOME.ply`: per faccia `fx fy fz` (fibra),
  `CV_L CV_T` (mm/ms), `elem_tag`.
- `output/aniso_zone_N/{faccia,roney}/sub_*/errors.csv`: una riga per faccia; colonne in `experiments/README.md`.

## Rigenerare la LA e le fibre dall'SSM (facoltativo)

Serve solo a verificare la provenienza di `data/augmenta/`; `run_all.sh` non
lo richiede. La procedura e' scritta dal sorgente di AugmentA e non e' stata
rieseguita da capo; il commit di AugmentA della run originale non e' noto.

```bash
python simulator/stage1/make_mean.py      # -> data/stage1/ssm_mean_surface.vtk
python simulator/stage1/extract_la.py     # -> data/stage1/ssm_mean_LA.vtk, 21811 vertici / 43270 triangoli
cp simulator/stage1/ssm_mean_LA_apexes.csv data/stage1/
```

`extract_la.py` riscrive `simulator/stage1/la_vertex_mask.npy`: `git status`
deve restare pulito. Poi AugmentA (<https://github.com/KIT-IBT/AugmentA>), con
il `docker/Docker-compose.yml` del suo repository e il volume dei dati puntato
su `data/stage1/`:

```yaml
    volumes:
      - /percorso/assoluto/di/cardiofiber/data/stage1:/headless/data
      - ../:/headless/augmenta
```

```bash
docker build -f Dockerfile-x64 -t augmenta_container_image .   # dalla cartella docker/ di AugmentA
docker compose -f Docker-compose.yml up
docker exec -u root CONTAINER pip install ruamel.yaml           # CONTAINER da docker ps
```

Nel desktop del container (VNC `localhost:5901`, o `localhost:6901/?password`):

```bash
cd /headless/augmenta
python main.py --mesh /headless/data/ssm_mean_LA.vtk --closed_surface 0 --open_orifices 0 --atrium LA --SSM_fitting 0 --resample_input 0 --normals_outside 1 --add_bridges 0 --MRI 0 --ofmt vtk --apex-file /headless/data/ssm_mean_LA_apexes.csv
```

`--resample_input 0` e' essenziale: tiene la mesh negli indici dell'SSM.
`--apex-file` sostituisce la selezione interattiva dell'apice
dell'auricola fatta nella run originale (2026-08-20): il punto piu' vicino
alle coordinate del file e' il vertice 1405, lo stesso. La config effettiva
della run originale e' `simulator/stage1/augmenta_config.yaml`. Le fibre
escono in `data/stage1/ssm_mean_LA_fibers/result_LA/LA_endo_with_fiber.{pts,elem,lon}`:
vanno confrontate con quelle di `data/augmenta/`, senza sovrascriverle.
