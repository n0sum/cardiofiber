# Stima differenziale dell'anisotropia di conduzione atriale

Tesi triennale di geometry processing. Registrazione di mappe di
attivazione (LAT) tra geometrie atriali distinte tramite functional maps,
e stima del tensore di conduzione per via differenziale, validate contro
una ground truth sintetica nota.

Riferimento di confronto: Roney et al., *A technique for measuring
anisotropy in atrial conduction to estimate conduction velocity and atrial
fibre direction*, Computers in Biology and Medicine 104 (2019) 278-290.

---

## Organizzazione

Le cartelle rispecchiano gli stadi della pipeline.

```
simulator/        stadi 1-2   produce geometrie e mappe LAT
estimator/        stadi 3-5   registra le mappe e stima il tensore
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
si scarica, e `checksums.sha256` ne fissa il contenuto esatto.

### Albero

```
simulator/
  stage1/             geometria di riferimento, fibre, mappe LAT
    opencarp/         file di parametri e script di lancio eikonale
  stage2/             forme SSM e remeshing  (da fare)
  convert/            da formati openCARP/AugmentA a PLY

estimator/
  src/                C++: gradiente per faccia, stima del tensore
  include/            happly vendorizzato (vedi THIRD_PARTY.md)
  fmap/               prototipi Python dello stadio 3, in attesa di
                      essere integrati nel C++

analysis/             validazione e figure

experiments/
  aniso_const.conf    un file per esperimento
  region_cv/          mappature region_tag -> (CV_L, CV_T)

data/
  external/           SSM scaricato
  stage1/             intermedi dello stadio 1 e cartella di lavoro di AugmentA
  reference/          mesh di riferimento, fibre, mappe LAT
  augmenta/           output AugmentA, ingresso di run_all.sh
  clinical/           mesh cliniche, NON redistribuibili

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
  ubuntu20). Serve `/usr/local/bin/openCARP`. `carputils` non e'
  necessario: i file di parametri si scrivono a mano.
- **Docker**, per AugmentA.
- **CMake** >= 4.1, un compilatore C++20 e **glm** (header-only, l'unica
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
`biatrial_model_completeFit.h5` (md5 `96d8bf64bc943dcba71e3702a552fa75`).

Scaricarlo in `data/external/`, poi verificare dalla radice:

```bash
md5sum data/external/biatrial_model_completeFit.h5
cd data/external && sha256sum -c ../../checksums.sha256
```

### 2. Geometria di riferimento e fibre

Dalla radice del repository, con l'ambiente `tesi` attivo:

```bash
python simulator/stage1/inspect_h5.py     # opzionale: struttura dell'h5
python simulator/stage1/make_mean.py      # -> data/stage1/ssm_mean_surface.vtk
python simulator/stage1/split_atria.py    # -> data/stage1/component_{0,1}.vtk
python simulator/stage1/extract_la.py     # -> data/stage1/ssm_mean_LA.vtk, 21811 vertici / 43270 triangoli
```

`extract_la.py` riscrive anche `simulator/stage1/la_vertex_mask.npy`
(gli indici SSM dei vertici della LA, tracciato in git): dopo il lancio
`git status` deve restare pulito. La LA e' la componente con le 4 vene
polmonari (`LA_REGION` in `extract_la.py`).

Poi le fibre, **a mano**, nel container AugmentA. Si usa il
`docker/Docker-compose.yml` del repository di AugmentA
(<https://github.com/KIT-IBT/AugmentA>), che monta il repository di
AugmentA in `/headless/augmenta`. Il suo volume dei dati va puntato su
`data/stage1/`:

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

La run originale e' del 2026-08-20, su un'altra macchina, con l'apice
dell'auricola scelto a mano con la selezione interattiva (punto 1405);
le sue coordinate sono in `simulator/stage1/ssm_mean_LA_apexes.csv`. Il
comando esatto di quella run e' in `simulator/stage1/run_augmenta.sh`.
Il suo ingresso `ssm_mean_LA.vtk` non e' stato conservato, ma si
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
python simulator/stage1/opencarp/run_sims.py experiments/region_cv/aniso_const.csv
                                               # -> data/reference/sims/aniso_const/run_NOME/act.igb
```

`mm_to_um.py` moltiplica per 1000 le sole coordinate di
`LA_endo_with_fiber.pts` (vedi *Unita' di misura*); `.elem` e `.lon` si
copiano invariati.

`run_sims.py` lancia `simulator/stage1/opencarp/run_eik.sh` una volta per
riga di `simulator/stage1/pacing_sites.csv`: sei run eikonali (FIM), un
sito di pacing ciascuno, ognuna in pochi decimi di secondo. Il punto di
stimolo e' il baricentro della cella `cell_id` sulla mesh in um,
arrotondato all'intero; le colonne di coordinate del CSV sono solo
informative. L'ordine delle righe del CSV e' l'ordine delle mappe
`lat_k` in `LA_LAT_SCENARIO.ply`. Modalita' eikonale pura: le velocita'
si impongono direttamente con `gregion[].dream.vel_l / vel_t / vel_n`, in
**mm/s**.

Le CV hanno una sola fonte, il CSV dello scenario in
`experiments/region_cv/` (`region_tag,CV_L,CV_T`, in **mm/ms**): il nome
dello scenario e' il nome del CSV. `run_sims.py` ne costruisce una
gregion per tag, con `vel_l = CV_L*1000` e `vel_t = vel_n = CV_T*1000`;
e' l'unico punto della conversione. `run_eik.sh` non contiene velocita'.
Lo stesso CSV va in `region_cv_csv` del `.conf`.

### 4. Conversione in PLY

```bash
python simulator/convert/multiLatToPly.py aniso_const   # -> data/reference/LA_LAT_aniso_const.ply
python simulator/convert/fibersToPly.py      # -> data/reference/LA_fibers_gt.ply
```

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
```

L'eseguibile si lancia **dalla radice del repository**: i percorsi nei file
di configurazione sono relativi a essa.

### 6. Analisi

```bash
python analysis/analisi.py aniso_const
```

L'argomento e' il nome della cartella in `output/`; le figure vanno in
`output/aniso_const/figures/`.

```bash
python analysis/regional.py aniso_regional aniso_const
```

Scenario con CV diverse per regione contro CV uniformi: stimato contro
imposto per regione, facce interne contro di confine (una faccia e' di
confine se condivide un vertice con una faccia di tag diverso), modello
nullo, maschera ai confini e confronto appaiato di `err_angle`. Solo
tabelle, nessuna figura.

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
costante, anisotropia variabile per regione e dati registrati via
functional maps: tre file di configurazione invece di tre versioni del
programma.

---

## Ground truth

| grandezza          | valore                                             |
|--------------------|----------------------------------------------------|
| CV longitudinale   | 1050 mm/s  (= 1.05 mm/ms)                          |
| CV trasversale     | 525 mm/s   (= 0.525 mm/ms)                         |
| rapporto           | 2.0                                                |
| fibre              | AugmentA/LDRBM, per cella                          |
| corrispondenza     | identita' di indice (SSM) o proiezione (remeshing)  |

Le fibre sono ground truth **imposta**, non fibre reali di paziente:
servono note e controllabili, non anatomicamente corrette.

### Tag di regione

```
1 = corpo LA (30537)      2 = anello mitralico (1715)
3 = LIPV (968)            4 = LSPV (1194)
5 = RIPV (718)            6 = RSPV (989)
7 = LAA (7149)
```

### Unita' di misura

Mesh e SSM in **millimetri**. openCARP vuole **micrometri**: la conversione
x1000 riguarda le sole coordinate. I LAT sono in millisecondi, quindi i
gradienti sono in ms/mm e le velocita' stimate escono in mm/ms.

Gli errori di unita' qui sono **silenziosi**: nessun messaggio, risultati
sbagliati di tre ordini di grandezza. Controllo rapido: l'estensione della
LA di riferimento e' circa [48, 71, 83] mm.

### Dati clinici

Le mesh cliniche (`2-LA.mesh`, `2-LA-FA.mesh`, `PRE_MESH.mesh`,
`POST_MESH.mesh`) non sono in questo repository e vanno rimesse in
`data/clinical/`. Entrano in gioco solo allo stadio 6, per il confronto
conduzione normale contro fibrillazione atriale, e non fanno parte dello
scheletro sintetico.

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
