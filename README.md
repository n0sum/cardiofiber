# Stima differenziale dell'anisotropia di conduzione atriale

Tesi triennale di geometry processing. Stima del tensore di conduzione
per via differenziale da mappe di attivazione (LAT), validata contro una
ground truth sintetica nota. La registrazione delle mappe tra geometrie
atriali distinte tramite functional maps e' uno stadio futuro, non
implementato: in `estimator/fmap/` ci sono solo prototipi.

Riferimento di confronto: Roney et al., *A technique for measuring
anisotropy in atrial conduction to estimate conduction velocity and atrial
fibre direction*, Computers in Biology and Medicine 104 (2019) 278-290.

---

## Percorso rapido

Dettagli e verifiche nei passi di *Riprodurre tutto da zero*.

```bash
git clone https://github.com/n0sum/cardiofiber
cd cardiofiber
conda env create -f environment.yml
conda activate tesi
```

1. openCARP v19.0, vedi *Requisiti*.
2. Scaricare l'SSM in `data/external/` (passo 1).
3. AugmentA, a mano nel container Docker, fino a `data/augmenta/` (passo 2).
4. Dalla radice, con l'ambiente `tesi` attivo:
   ```bash
   ./run_all.sh
   ```

## Risultati attesi

`output/aniso_const/figure/numeri.csv` della run di riferimento:

| grandezza                      | valore  |
|--------------------------------|---------|
| angolo di fibra, q50           | 0.747°  |
| angolo di fibra, q90           | 3.54°   |
| area non stimabile             | 0.53 %  |
| coppie (faccia, mappa) escluse | 2.35 %  |

Una versione diversa di AugmentA puo' spostare leggermente questi numeri
(le fibre sono l'ingresso di tutto il resto).

---

## Organizzazione

Le cartelle rispecchiano gli stadi della pipeline.

```
simulator/        stadi 1-2   produce geometrie e mappe LAT
estimator/        stadi 3-5   stima il tensore (registrazione: stadio futuro)
analysis/         stadio 6    valida contro la ground truth, confronta
experiments/                  cosa calcolare in una data run
data/             ignorato    ingressi e prodotti
output/           ignorato    risultati
papers/           ignorato    bibliografia, non redistribuibile
```

In git stanno il codice e la configurazione, non i dati. Chi clona non
riceve nessun risultato: li ricostruisce. Questo e' il test della
riproducibilita'.

L'unico ingresso che non si rigenera e' lo statistical shape model:
si scarica, e md5 e sha256 nel passo 1 ne fissano il contenuto esatto.

### Albero

```
simulator/
  stage1/             geometria di riferimento, fibre, mappe LAT
    opencarp/         mesh in um, tag casuali, lancio eikonale
  stage2/             forme SSM e remeshing  (da fare)
  convert/            da formati openCARP/AugmentA a PLY

estimator/
  src/                C++: gradiente per faccia, stima del tensore
  include/            happly vendorizzato (vedi THIRD_PARTY.md)
  fmap/               prototipi Python dello stadio 3 (functional maps),
                      non usati dalla pipeline

analysis/             validazione e figure

experiments/
  aniso_const.conf    un file per esperimento
  region_cv/          scenari: region_tag -> (CV_L, CV_T)

data/
  external/           SSM scaricato
  stage1/             intermedi dello stadio 1 e cartella di lavoro di AugmentA
  reference/          mesh di riferimento, fibre, mappe LAT
  augmenta/           output AugmentA, ingresso di run_all.sh
  clinical/           mesh cliniche, non distribuite e non usate

output/               CSV di errori per faccia, figure
```

`data/` e `output/` non esistono in un clone pulito: le creano gli script.

---

## Requisiti

- Linux. Sviluppato su Ubuntu 26.04.1, Python 3.13.
- **conda** con l'ambiente di `environment.yml`:
  ```bash
  conda env create -f environment.yml
  conda activate tesi
  ```
- **openCARP v19.0**, installato da `.deb` (pacchetto standard, non
  ubuntu20): <https://opencarp.org/download/installation>. Serve
  `/usr/local/bin/openCARP`. `carputils` non e' necessario e non ci sono
  file `.par`: i parametri stanno nella riga di comando di
  `simulator/stage1/opencarp/run_eik.sh`, piu' regioni e velocita' che
  aggiunge `run_sims.py`.
- **Docker**, per AugmentA.
- **CMake** >= 3.20, un compilatore C++20 e **glm** (header-only, l'unica
  dipendenza C++ esterna; su Ubuntu `libglm-dev`). Versioni testate:
  GCC 15.2.0, CMake 4.2.3, libglm-dev 1.0.1.
- ParaView, opzionale, per l'ispezione visiva.

---

## Riprodurre tutto da zero

> La catena ha **un anello non automatizzabile**: AugmentA gira solo in un
> container Docker con sessione desktop. Il passo 2 va eseguito a mano.
> I passi 3-6 li esegue `./run_all.sh`.

### 1. Statistical shape model

SSM biatriale del KIT-IBT (Nagel et al., CC-BY).

Zenodo: <https://zenodo.org/records/4309958>, file
`biatrial_model_completeFit.h5` (md5 `96d8bf64bc943dcba71e3702a552fa75`,
sha256 `d4b3a90f2b9c6db79bb225954a4e5eeaed59f44916e5a93935dde20458da63b8`).

Scaricarlo in `data/external/`, poi verificare dalla radice:

```bash
md5sum data/external/biatrial_model_completeFit.h5
sha256sum data/external/biatrial_model_completeFit.h5
```

### 2. Geometria di riferimento e fibre

Dalla radice del repository, con l'ambiente `tesi` attivo:

```bash
python simulator/stage1/make_mean.py      # -> data/stage1/ssm_mean_surface.vtk
python simulator/stage1/split_atria.py    # -> data/stage1/component_{0,1}.vtk
python simulator/stage1/extract_la.py     # -> data/stage1/ssm_mean_LA.vtk, 21811 vertici / 43270 triangoli
```

`extract_la.py` riscrive anche `simulator/stage1/la_vertex_mask.npy`
(gli indici SSM dei vertici della LA, tracciato in git): dopo il lancio
`git status` deve restare pulito. La LA e' la componente con le 4 vene
polmonari (`LA_REGION` in `extract_la.py`).

Poi le fibre, **a mano**, nel container AugmentA.

> Questa procedura e' scritta dal sorgente di AugmentA e non e' ancora
> stata rieseguita da capo.

Il commit di AugmentA della run originale non e' registrato. Commit da
usare: `<COMMIT DA FISSARE>`.

```bash
git clone https://github.com/KIT-IBT/AugmentA
cd AugmentA
git checkout <COMMIT DA FISSARE>
```

Si usa il `docker/Docker-compose.yml` del repository di AugmentA, che
monta il repository di AugmentA in `/headless/augmenta`. Il suo volume dei
dati va puntato su `data/stage1/`:

```yaml
    volumes:
      - /percorso/assoluto/di/cardiofiber/data/stage1:/headless/data
      - ../:/headless/augmenta
```

Immagine (una volta) e avvio, dalla cartella `docker/` di AugmentA:

```bash
docker build -f Dockerfile-x64 -t augmenta_container_image .
docker compose -f Docker-compose.yml up
```

Nel container manca `ruamel.yaml`, che AugmentA importa. Con il container
avviato (il nome lo da' `docker ps`):

```bash
docker exec -u root CONTAINER pip install ruamel.yaml
```

Dalla radice di questo repository si mette l'apice accanto alla mesh:

```bash
cp simulator/stage1/ssm_mean_LA_apexes.csv data/stage1/
```

Desktop del container: VNC su `localhost:5901`, oppure il browser su
`localhost:6901/?password`. In un terminale del container:

```bash
cd /headless/augmenta
python main.py --mesh /headless/data/ssm_mean_LA.vtk --closed_surface 0 --open_orifices 0 --atrium LA --SSM_fitting 0 --resample_input 0 --normals_outside 1 --add_bridges 0 --MRI 0 --ofmt vtk --apex-file /headless/data/ssm_mean_LA_apexes.csv
```

`--apex-file` evita la selezione interattiva dell'apice: AugmentA legge
le coordinate dal CSV e prende il punto della mesh piu' vicino. Il
percorso e' assoluto perche' la cwd e' `/headless/augmenta`. Il punto
piu' vicino alle coordinate del file e' il 1405, quello della run
originale (verificato).

AugmentA scrive dentro `data/stage1/`: il log
`ssm_mean_LAAugmentA.log`, la config effettiva
`_configs/augmenta_config.yaml`, gli anelli in `ssm_mean_LA_surf/` e le
fibre in `ssm_mean_LA_fibers/result_LA/`. In `data/augmenta/` vanno
copiati i file che legge `run_all.sh`, piu' log e config:

```bash
mkdir -p data/augmenta
cp data/stage1/ssm_mean_LA_fibers/result_LA/LA_endo_with_fiber.{pts,elem,lon} data/augmenta/
cp data/stage1/ssm_mean_LAAugmentA.log data/stage1/_configs/augmenta_config.yaml data/augmenta/
```

`data/augmenta/` contiene l'uscita gia' verificata: per controllare una
nuova run si confronta con quella, senza sovrascriverla.

La run originale e' del 2026-08-20, su un'altra macchina, con lo stesso
comando senza `--apex-file`: l'apice dell'auricola e' stato scelto a mano
con la selezione interattiva (vertice 1405), e le sue coordinate sono in
`simulator/stage1/ssm_mean_LA_apexes.csv`. Il suo ingresso `ssm_mean_LA.vtk` non e' stato conservato, ma si
rigenera identico: la mesh prodotta da `extract_la.py` coincide con
`data/augmenta/LA_endo_with_fiber.*` vertice per vertice (scarto
massimo 5e-9 mm, arrotondamento del `.pts`) e triangolo per triangolo,
nello stesso ordine.

`--resample_input 0` e' essenziale: mantiene l'output nello spazio di
indici dell'SSM, quindi non serve nessuna proiezione a valle.

**Verso dei triangoli.** AugmentA inverte il verso di tutti i triangoli
rispetto all'SSM: `elem == ssm_tri[:, ::-1]`. La connettivita' di
riferimento per tutte le forme e' `LA_endo_with_fiber.elem`, mai
`representer/cells`.

### 3. Mesh in micrometri e mappe LAT

Da qui in poi tutto e' scriptato: `./run_all.sh`, con l'ambiente `tesi`
attivo, esegue i passi 3-6 in ordine a partire da `data/augmenta/`. Se
`data/augmenta/` non c'e' si ferma: serve prima AugmentA.

```bash
python simulator/stage1/opencarp/mm_to_um.py   # -> data/reference/LA_um.{pts,elem,lon}
python simulator/stage1/opencarp/random_tags.py aniso_random 0
                                               # -> data/reference/aniso_random/LA_um.{pts,elem,lon}
python simulator/stage1/opencarp/run_sims.py experiments/region_cv/aniso_const.csv data/reference/LA_um
                                               # -> data/reference/sims/aniso_const/run_k/act.igb
```

`mm_to_um.py` moltiplica per 1000 le sole coordinate di
`LA_endo_with_fiber.pts` (vedi *Unita' di misura*); `.elem` e `.lon` si
copiano invariati.

`random_tags.py SCENARIO SEED` assegna a ogni triangolo di
`data/reference/LA_um` uno dei tag di `experiments/region_cv/SCENARIO.csv`,
a caso con la stessa probabilita' (`numpy.random.default_rng(SEED)`), e
scrive la mesh in `data/reference/SCENARIO/LA_um`: stessi punti, fibre e
triangoli, nuovo tag. Il seme (0) sta solo in `run_all.sh`.

`run_sims.py SCENARIO.csv MESH_BASE` simula lo scenario sulla mesh in um
`MESH_BASE` (senza estensione), i cui tag devono essere quelli del CSV:
`data/reference/LA_um` per `aniso_const` e `aniso_regional`,
`data/reference/aniso_random/LA_um` per `aniso_random`.

`run_sims.py` lancia `simulator/stage1/opencarp/run_eik.sh` una volta per
riga di `simulator/stage1/pacing_sites.csv`: sei run eikonali (FIM), un
sito di pacing ciascuno, ognuna in pochi decimi di secondo. Il CSV ha la
sola colonna `cell_id`: il punto di stimolo e' il baricentro di quella
cella sulla mesh in um, arrotondato all'intero. La riga k del CSV e' la
run `run_k` e la mappa `lat_k` in `LA_LAT_SCENARIO.ply`. Modalita' eikonale pura: le velocita'
si impongono direttamente con `gregion[].dream.vel_l / vel_t / vel_n`, in
**mm/s**.

Le CV hanno una sola fonte, il CSV dello scenario in
`experiments/region_cv/` (`region_tag,CV_L,CV_T`, in **mm/ms**): il nome
dello scenario e' il nome del CSV. `run_sims.py` ne costruisce una
gregion per tag, con `vel_l = CV_L*1000` e `vel_t = vel_n = CV_T*1000`;
e' l'unico punto della conversione. Dagli stessi tag costruisce l'unico
dominio eikonale (`phys_region`). `run_eik.sh` non contiene ne' velocita'
ne' tag. Lo stesso CSV da' le CV vere per faccia della ground truth
(passo 4). Gli scenari sono in *Ground truth*; `run_all.sh` li simula e
stima tutti e tre.

### 4. Conversione in PLY

```bash
python simulator/convert/multiLatToPly.py aniso_const                        # -> data/reference/LA_LAT_aniso_const.ply
python simulator/convert/fibersToPly.py aniso_const data/reference/LA_um     # -> data/reference/LA_GT_aniso_const.ply
```

`LA_GT_SCENARIO.ply` ha, per faccia, la fibra normalizzata (`fx fy fz`),
le CV vere dello scenario (`CV_L CV_T`, dal CSV per tag) e `elem_tag`. Lo
stimatore la legge cosi' com'e'. `fibersToPly.py SCENARIO MESH_BASE`
prende geometria e fibre da `data/augmenta/` e i tag per faccia da
`MESH_BASE.elem`, la stessa mesh passata a `run_sims.py`; si ferma se i
triangoli di `MESH_BASE.elem` non coincidono con quelli di
`data/augmenta/`.

Entrambi leggono le fonti canoniche (`.pts`, `.elem`, `.lon`, `act.igb`) e
**non** i `.vtk`: passare per pyvista puo' attraversare filtri che
riordinano le celle, e un disallineamento faccia-fibra non produce nessun
errore visibile, solo un campo di fibre sbagliato.

### 5. Stima e validazione

```bash
cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build
./estimator/build/Tesi0 experiments/aniso_const.conf
./estimator/build/Tesi0 experiments/aniso_const_nomask.conf
./estimator/build/Tesi0 experiments/aniso_regional.conf
./estimator/build/Tesi0 experiments/aniso_random.conf
```

L'eseguibile si lancia **dalla radice del repository**: i percorsi nei file
di configurazione sono relativi a essa.

### 6. Analisi

```bash
python analysis/analisi.py aniso_const
python analysis/analisi.py aniso_const_nomask
python analysis/analisi.py aniso_regional
python analysis/analisi.py aniso_random
```

L'argomento e' il nome della cartella in `output/`. Statistiche pesate per
area, sulle facce stimabili, di tre errori: angolo di fibra, errore
relativo con segno su CV_L e sul rapporto CV_L / CV_T. Scrive in
`output/ESPERIMENTO/figure/`:

- `errori.png`: per ogni errore la curva cumulativa (% di area stimabile
  con errore <= x), con le linee al 50% e al 90%; nel titolo l'area non
  stimabile e le coppie (faccia, mappa) escluse;
- `numeri.csv`: area non stimabile, coppie escluse, q50 e q90 di ogni
  errore (valore assoluto per i relativi) e mediana con segno dei relativi.

Le coppie escluse sono calcolate su K = righe di `pacing_sites.csv`: vale
per esperimenti con `maps = all`.

Poi i sottoinsiemi di mappe (vedi `experiments/README.md`):

```bash
python experiments/subsets.py
./estimator/build/Tesi0 output/subsets/confs/*.conf
python analysis/analisi.py --sottoinsiemi
```

---

## Esperimenti

Un file `experiments/*.conf` per esperimento. L'eseguibile ne riceve uno o
piu' e scrive un CSV di errori per faccia per esperimento, con
`experiment_id` come prima colonna. Formato dei `.conf` e colonne del CSV
in `experiments/README.md`.

Le soglie delle maschere e i percorsi stanno nel file di configurazione e
non nel codice, perche' **determinano i risultati**: sono dati, non
implementazione. Se stessero nel codice, un risultato non sarebbe piu'
attribuibile a una configurazione, e cambiare una soglia richiederebbe di
ricompilare.

E' anche cio' che rende possibile calcolare in una sola run anisotropia
costante, a tratti per tag e casuale per faccia: tre file di
configurazione invece di tre versioni del programma.

---

## Ground truth

| grandezza          | valore                                             |
|--------------------|----------------------------------------------------|
| CV                 | per tag, dallo scenario (tabella sotto)            |
| fibre              | AugmentA/LDRBM, per cella                          |
| corrispondenza     | identita' di indice (SSM)                          |

Le fibre sono ground truth **imposta**, non fibre reali di paziente:
servono note e controllabili, non anatomicamente corrette.

### Scenari

Le CV vere dipendono dallo scenario: `experiments/region_cv/SCENARIO.csv`
(`region_tag,CV_L,CV_T`, mm/ms), uno per `.conf` di `experiments/`.

| scenario         | tag                         | CV                                                        |
|------------------|-----------------------------|-----------------------------------------------------------|
| `aniso_const`    | 1-7 di AugmentA             | uniformi, CV_L 1.05, CV_T 0.525 mm/ms (r = 2)             |
| `aniso_regional` | 1-7 di AugmentA             | costanti a tratti per tag, costrutto di test (r da 1.25 a 2.5) |
| `aniso_random`   | 1-25, a caso per faccia     | CV_L = 1.05 x (0.8 ... 1.2), r = 2.0 x (0.8 ... 1.2), tutte le 25 combinazioni; tag 13 = `aniso_const` |

### Tag di regione

I tag di `LA_endo_with_fiber.elem` vanno da 1 a 7; si usano solo come
indici. Facce per tag (`aniso_const`, `aniso_regional`):

```
1: 30537   2: 1715   3: 968   4: 1194   5: 718   6: 989   7: 7149
```

In `aniso_random` i tag sono quelli di `random_tags.py`, circa 1/25
delle 43270 facce ciascuno (~1730).

### Unita' di misura

Mesh e SSM in **millimetri**. openCARP vuole **micrometri**: la conversione
x1000 riguarda le sole coordinate. I LAT sono in millisecondi, quindi i
gradienti sono in ms/mm e le velocita' stimate escono in mm/ms.

Gli errori di unita' qui sono **silenziosi**: nessun messaggio, risultati
sbagliati di tre ordini di grandezza. Controllo rapido: l'estensione della
LA di riferimento e' circa [48, 71, 83] mm.

### Dati clinici

Non servono per i risultati attuali e non sono distribuiti. Le mesh
cliniche (`2-LA.mesh`, `2-LA-FA.mesh`, `PRE_MESH.mesh`, `POST_MESH.mesh`)
sono previste per un confronto futuro conduzione normale contro
fibrillazione atriale, in `data/clinical/`; nessuno script attuale le legge.

---

## Limiti noti

- **Eikonale pura, senza termine di curvatura.** Il modello e' meno
  accurato dove il fronte e' molto curvo, cioe' vicino ai siti di stimolo.
  Questo non tocca la validazione: la ground truth sintetica e' generata
  dallo stesso modello, quindi l'errore di stima vicino ai siti e' di
  discretizzazione (il gradiente lineare per faccia non segue la curvatura
  del fronte), non di fisica mancante. La scelta dell'eikonale e'
  necessaria, non di comodo: nel monodominio la CV va calibrata e
  **dipende dalla risoluzione della mesh**, che e' uno degli assi
  sperimentali di questa tesi.
- **Le collisioni di fronte sono intrinseche all'equazione**, non artefatti
  della mesh. Dove due fronti si incontrano T ha una piega e il gradiente
  non esiste. Vanno mascherate; nessun raffinamento le elimina.
- **Il modello sintetico non e' un modello elettrofisiologico validato**:
  e' un generatore controllato di campi LAT con proprieta' note, che e'
  esattamente cio' che serve per validare uno stimatore.
- **Le forme generate dall'SSM non sono meccanica cardiaca**: sono un
  costrutto geometrico con non-isometria dosata.
- **I modi dell'SSM sono biatriali**: si deforma tutto e si tiene la sola
  parte LA.
- **Il confronto con Roney** va fatto implementando il loro metodo sui
  dati di questo repository, non confrontando numeri calcolati su dataset
  diversi.
