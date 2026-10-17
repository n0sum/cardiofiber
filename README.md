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
  reference/          mesh di riferimento, fibre, mappe LAT
  augmenta/           output AugmentA
  clinical/           mesh cliniche, NON redistribuibili

output/               CSV di errori per faccia, figure
docs/                 note di progetto, dump dell'ambiente
```

`data/` e `output/` non esistono in un clone pulito: le creano gli script.

---

## Requisiti

- Linux. Sviluppato su Pop!_OS 24.04.
- **conda** con l'ambiente di `environment.yml`:
  ```bash
  conda env create -f environment.yml
  conda activate tesi
  ```
- **openCARP v19.0**, installato da `.deb` (pacchetto standard, non
  ubuntu20). Serve `/usr/local/bin/openCARP`. `carputils` non e'
  necessario: i file di parametri si scrivono a mano.
- **Docker**, per AugmentA.
- **CMake** >= 3.20 e un compilatore C++20.
- ParaView, opzionale, per l'ispezione visiva.

---

## Riprodurre tutto da zero

> La catena ha **un anello non automatizzabile**: AugmentA gira solo in un
> container Docker con sessione desktop. Il passo 2 va eseguito a mano.
> Tutto il resto e' scriptato.

### 1. Statistical shape model

SSM biatriale del KIT-IBT (Nagel et al., CC-BY).

<!-- DA COMPLETARE: URL esatto di download e versione del dataset -->

Scaricare `biatrial_model_completeFit.h5` in `data/external/`, poi
verificare:

```bash
cd data/external && sha256sum -c ../../checksums.sha256
```

### 2. Geometria di riferimento e fibre

```bash
conda activate tesi
python simulator/stage1/make_mean.py      # forma media dall'SSM
python simulator/stage1/split_atria.py    # separa LA da RA
python simulator/stage1/extract_la.py     # LA, 21811 vertici / 43270 triangoli
```

Poi le fibre, **a mano**, nel container AugmentA:

<!-- DA COMPLETARE: comando docker run esatto, e dove montare la mesh
     dentro il container -->

```bash
bash simulator/stage1/run_augmenta.sh
```

Flag usati: `--closed_surface 0 --open_orifices 0 --atrium LA
--SSM_fitting 0 --resample_input 0 --normals_outside 1 --add_bridges 0
--MRI 0 --ofmt vtk`, con `--apex-file ssm_mean_LA_apexes.csv` per evitare
la selezione interattiva dell'apice dell'auricola.

`--resample_input 0` e' essenziale: mantiene l'output nello spazio di
indici dell'SSM, quindi non serve nessuna proiezione a valle.

### 3. Mappe LAT

```bash
bash simulator/stage1/opencarp/run_eik.sh
```

Sei run eikonali (FIM), un sito di pacing ciascuno, definiti in
`simulator/stage1/pacing_sites.csv`. Modalita' eikonale pura: le velocita'
si impongono direttamente con `gregion[].dream.vel_l / vel_t / vel_n`,
in **mm/s**.

### 4. Conversione in PLY

```bash
python simulator/convert/multiLatToPly.py    # -> data/reference/LA_LAT.ply
python simulator/convert/fibersToPly.py      # -> data/reference/LA_fibers_gt.ply
```

Entrambi leggono le fonti canoniche (`.pts`, `.elem`, `.lon`, `.npz`) e
**non** i `.vtk`: passare per pyvista puo' attraversare filtri che
riordinano le celle, e un disallineamento faccia-fibra non produce nessun
errore visibile, solo un campo di fibre sbagliato.

### 5. Stima e validazione

```bash
cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build
./estimator/build/Tesi0 experiments/aniso_const.conf
```

L'eseguibile si lancia **dalla radice del repository**: i percorsi nei file
di configurazione sono relativi a essa.

### 6. Analisi

```bash
python analysis/analisi.py
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
  accurato dove il fronte e' molto curvo, cioe' vicino ai siti di stimolo:
  e' la ragione fisica, non solo numerica, della maschera di prossimita'.
  La scelta dell'eikonale e' necessaria, non di comodo: nel monodominio la
  CV va calibrata e **dipende dalla risoluzione della mesh**, che e' uno
  degli assi sperimentali di questa tesi.
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

---

## Convenzione di provenienza

Ogni cartella generata sotto `data/` porta un `provenance.json` che
dichiara da quale script, con quali parametri, da quale commit e da quali
ingressi e' stata prodotta. Template in `PROVENANCE_TEMPLATE.json`.

Senza, a distanza di mesi non e' piu' possibile sapere cosa ha prodotto
cosa, e un lettore esterno non puo' saperlo affatto.
