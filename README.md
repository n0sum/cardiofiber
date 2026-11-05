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
analysis/analisi.py            errori pesati per area, figure, numeri.csv
experiments/                   modello.conf e roney.conf, CV per zona (region_cv/),
                               subsets.py, zone.py, mappe.py; formati in experiments/README.md
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
compila `estimator/` in `estimator/build/` e calcola stime e confronti.
Uscite, per esperimento, in `output/`:

| cartella                                   | contenuto                                          |
|--------------------------------------------|----------------------------------------------------|
| `zone/`                                    | scenari `aniso_zone_N`: `errors.csv` per faccia, `figure/{zone.png,numeri.csv}`; `aniso_zone_1/figure/{errori.png,numeri.csv}` |
| `subsets/`                                 | sottoinsiemi di mappe su `aniso_zone_1`, `figure/numeri.csv` e grafici |
| `roney/`                                   | Roney contro stimatore per faccia, `figure/{zona_1.png,zone.png,numeri.csv,terne.csv}` |
| `mappe/`                                   | PLY degli errori per faccia su `aniso_zone_1` (tutte le mappe e terne scelte), uno per errore, per MeshLab |

Scenari: `aniso_zone_N`, la mesh divisa in N zone casuali (seme 0), ogni zona
con una delle 25 combinazioni di CV_L e r = CV_L/CV_T di
`region_cv/aniso_random.csv`, N = 1, 50, 500, 5000, 25000, 43270. Con N = 1
le CV sono uniformi (e' uscito il tag 16: CV_L 1.155, CV_T 0.722 mm/ms,
r = 1.6); con N = 43270 c'e' una zona per faccia. Dettagli in
`experiments/README.md`.

Il FIM di openCARP si ferma dopo 1000 iterazioni; le run che convergono ne
fanno 96-166. Una mappa che non converge non entra nel LAT PLY dello
scenario, e le altre tengono l'indice del loro sito. Mappe che non
convergono: la 5 in `aniso_zone_1` e `aniso_zone_5000`, la 0 in
`aniso_zone_50`, la 3 in `aniso_zone_43270`; negli altri scenari ci sono tutte e sei.

## Risultati attesi

Errori pesati per area sulle facce stimabili: angolo di fibra in gradi,
errori relativi in valore assoluto, in %, su CV_L e sul rapporto r.

**aniso_zone_1**, tutte le mappe, cioe' 0-4 (`output/zone/aniso_zone_1/figure/numeri.csv`):

| area non stimabile | angolo q50 | angolo q90 | CV_L q50 | CV_L q90 | r q50 | r q90 |
|--------|--------|-------|--------|--------|--------|---------|
| 0.33 % | 0.975° | 4.95° | 0.66 % | 3.69 % | 1.41 % | 7.81 %  |

**Sottoinsiemi di mappe** su `aniso_zone_1` (`output/subsets/figure/numeri.csv`,
confronto appaiato con tutte e cinque le mappe): con 3 mappe, sulle 10
terne, area non stimabile 4.0-9.6 % (mediana 5.6 %), angolo q50 1.17-1.76°
(mediana 1.37°), angolo q90 9.2-14.9° (mediana 10.5°); con 5 mappe i valori
di `aniso_zone_1`.

**Zone** (`output/zone/figure/numeri.csv`), d = diametro medio di zona:

| N     | d [mm] | mappe | area non stimabile | angolo q50 | angolo q90 |
|-------|--------|-------|--------------------|------------|------------|
| 43270 | 0.59   | 5     | 3.85 %             | 5.90°      | 22.1°      |
| 25000 | 0.78   | 6     | 1.58 %             | 4.71°      | 17.4°      |
| 5000  | 1.74   | 5     | 2.22 %             | 3.54°      | 16.4°      |
| 500   | 5.51   | 6     | 0.93 %             | 1.71°      | 8.89°      |
| 50    | 17.4   | 5     | 1.63 %             | 1.14°      | 7.23°      |
| 1     | 123    | 5     | 0.33 %             | 0.975°     | 4.95°      |

**Roney su aniso_zone_1**, terne di mappe (`output/roney/figure/numeri.csv`;
errori sulle facce stimabili in entrambi i metodi, 76.2 % dell'area):

| metodo  | area non stimabile | angolo q50 | angolo q90 | CV_L q50 | r q50  |
|---------|--------------------|------------|------------|----------|--------|
| Roney   | 21.9 %             | 6.60°      | 21.5°      | 13.5 %   | 11.1 % |
| faccia  | 5.8 %              | 1.30°      | 9.12°      | 0.82 %   | 2.06 % |

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
  `lat_k` (ms), una per sito k la cui mappa converge
  (`data/reference/sims/NOME/mappe.csv`). `LA_GT_NOME.ply`: per faccia `fx fy fz` (fibra),
  `CV_L CV_T` (mm/ms), `elem_tag`.
- `output/**/errors.csv`: una riga per faccia; colonne in `experiments/README.md`.

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
