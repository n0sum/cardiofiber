# Esperimenti

Un file `.conf` per esperimento. L'eseguibile ne riceve uno o piu' e
scrive un CSV di errori per faccia per ciascuno:

```bash
./estimator/build/Tesi0 experiments/a.conf experiments/b.conf
```

Tutti i file vengono letti e validati prima di eseguire il primo.

## Perche' i parametri stanno qui e non nel codice

Le soglie delle maschere determinano i risultati: sono dati, non
implementazione. Se stanno nel codice, un risultato non e' piu'
attribuibile a una configurazione, e cambiare una soglia richiede di
ricompilare e di perdere la tracciabilita' di cosa ha prodotto cosa.

## Formato

`chiave = valore`, una per riga. Righe vuote e righe che iniziano con `#`
si saltano. Spazi attorno a `=` sono tollerati. La descrizione
dell'esperimento sta nei commenti in testa al file.

| chiave                      | significato                                   |
|-----------------------------|-----------------------------------------------|
| `experiment_id`             | identificativo, prima colonna del CSV          |
| `lat_ply`                   | mesh + K mappe LAT: l'ingresso dello stimatore |
| `ground_truth_ply`          | fibre, CV vere e tag di regione per faccia     |
| `maps`                      | `all`, oppure indici separati da virgola       |
| `mask.min_singular_value`   | soglia sul valore singolare minimo del 3x3     |
| `mask.collision_angle_deg`  | soglia del criterio di collisione, gradi       |
| `output_csv`                | dove scrivere gli errori per faccia            |

Tutte le 7 chiavi sono obbligatorie.

`ground_truth_ply` e' `data/reference/LA_GT_NOME.ply` dello scenario
che ha prodotto `lat_ply`: fibra, CV_L, CV_T e tag per faccia, letti cosi'
come sono. Le soglie `mask.*` non hanno un default nel codice: esistono
solo qui.

## Scenari

Uno scenario e' un nome `NOME`, un CSV `region_cv/*.csv`
(`region_tag,CV_L,CV_T`, mm/ms) e la mesh in um con quei tag, `MESH_BASE`:

| scenario       | `MESH_BASE`                         | `.conf`                                       |
|----------------|-------------------------------------|-----------------------------------------------|
| `aniso_const`  | `data/reference/LA_um`              | `aniso_const.conf`, `aniso_const_nomask.conf` |
| `aniso_zone_N` | `data/reference/aniso_zone_N/LA_um` | generati da `zone.py`; `aniso_zone_43270_nomask.conf` |

`aniso_const` usa `aniso_const.csv`, tutti gli `aniso_zone_N` usano
`aniso_random.csv`.

`aniso_random.csv` ha 25 tag, tutte le combinazioni di 5 livelli di
CV_L = 1.05 x (0.8, 0.9, 1.0, 1.1, 1.2) e 5 di r = CV_L/CV_T = 2.0 x (0.8,
0.9, 1.0, 1.1, 1.2), CV_L esterno e r interno (tag 13 = `aniso_const`).
La mesh di `aniso_zone_N` la scrive `simulator/stage1/opencarp/random_tags.py
aniso_zone_N experiments/region_cv/aniso_random.csv 0 N`: N zone
attorno a N facce scelte a caso (faccia al centro piu' vicino), un tag a
caso per zona con la stessa probabilita'; con N = 43270 ogni faccia e' una
zona. Semi e N stanno solo in `run_all.sh`.
`aniso_zone_43270_nomask.conf` e' il `.conf` di `aniso_zone_43270` con
`mask.collision_angle_deg = 181`, come `aniso_const_nomask.conf`.

Per ogni scenario, dalla radice:

```bash
python simulator/stage1/opencarp/run_sims.py NOME CSV MESH_BASE   # -> data/reference/sims/NOME/
python simulator/convert/multiLatToPly.py NOME                    # -> lat_ply
python simulator/convert/fibersToPly.py NOME CSV MESH_BASE        # -> ground_truth_ply
```

`fibersToPly.py` prende geometria e fibre da `data/augmenta/` e i tag da
`MESH_BASE.elem`, e si ferma se i triangoli non coincidono.

`maps = 0,2,4` tiene solo quelle mappe, nell'ordine dato; un indice
fuori range e' un errore.

Chiave sconosciuta (con `file:riga`) o mancante, file inesistente: errore,
nessun default silenzioso. Se una chiave compare piu' volte vale l'ultima.
La directory di `output_csv` viene creata se non esiste.

Tutti i percorsi sono relativi alla radice del repository, e l'eseguibile
si lancia da la'.

`experiment_id` finisce non quotato in ogni riga del CSV: non deve
contenere virgole.

## Maschere per coppia (faccia, mappa)

Prima della stima, ogni coppia (faccia `i`, mappa `k`) puo' essere esclusa
dal sistema lineare della faccia: si toglie la sola riga `k`, le altre
mappe della faccia restano. Le maschere si calcolano una volta sui
gradienti grezzi, solo da LAT, gradienti e geometria della mesh: niente
ground truth, CV vere, fibre o posizioni dei siti. In `main.cpp` vengono
calcolate prima che la ground truth sia caricata.

L'unico criterio e' la **collisione**: per ogni vicino `j` (adiacenza per
spigolo) il gradiente di `j` viene aperto nel piano di `i` attorno allo
spigolo condiviso, e si misura l'angolo pieno, in [0, 180] gradi, con
quello di `i`. La coppia e' esclusa se il massimo sui vicini supera
`mask.collision_angle_deg`. Gradiente nullo da' angolo 0, una faccia senza
vicini da' massimo 0.

Un criterio di stimolo (escludere le coppie attivate entro una soglia dal
minimo della mappa) e' stato provato e tolto: effetto trascurabile
sull'errore, e le righe a gradiente quasi nullo vicino al sito sono
innocue per i minimi quadrati, mentre il filtro toglieva righe utili
nella zona di tag 7.

Una faccia con meno di 3 righe rimaste e' non stimabile.

`mask.collision_angle_deg = 181` non esclude nulla, senza codice speciale:
l'angolo non supera mai 180. `aniso_const_nomask.conf` la usa come
regressione, `aniso_zone_43270_nomask.conf` per misurare l'effetto della
maschera sull'anisotropia casuale per faccia.

Le LAT negative (FIM non convergito, -1) vengono rifiutate a monte, da
`simulator/convert/multiLatToPly.py`, prima di scrivere `lat_ply`.

## CSV di uscita

Una riga per faccia, nessuna aggregazione. Colonne, nell'ordine:

| colonna         | significato                                            |
|-----------------|--------------------------------------------------------|
| `experiment_id` | dal `.conf`, uguale su tutte le righe                   |
| `face_id`       | indice della faccia nella mesh                          |
| `area`          | area della faccia, mm^2                                 |
| `region_tag`    | tag di regione dalla ground truth                       |
| `estimable`     | 1 se la faccia supera le maschere, 0 altrimenti         |
| `k_used`        | mappe rimaste dopo le maschere per coppia               |
| `cond`          | numero di condizionamento (vedi sotto); vuoto se non calcolato |
| `CV_L_est`      | CV longitudinale stimata, mm/ms; vuota se non stimabile |
| `CV_T_est`      | CV trasversale stimata, mm/ms; vuota se non stimabile   |
| `CV_L_true`     | CV longitudinale vera, mm/ms                            |
| `CV_T_true`     | CV trasversale vera, mm/ms                              |
| `err_angle`     | errore angolare di fibra, `acos(|a.b|)` in gradi; vuoto se non stimabile |

`cond` e' `sigma_max / sigma_min` del sistema della faccia, una riga per
mappa tenuta: `(g1^2, sqrt(2) g1 g2, g2^2)`, con `(g1, g2)` il gradiente
LAT nella base locale della faccia. Il `sqrt(2)` e' l'embedding di
Frobenius delle matrici simmetriche: una rotazione della base agisce come
rotazione 3x3, quindi `cond` e la soglia `mask.min_singular_value` non
dipendono dall'ordine dei vertici.

`experiment_id` sta su ogni riga perche' i CSV di esperimenti diversi si
possano concatenare e raggruppare senza perdere la provenienza:

```python
df = pd.concat(pd.read_csv(f) for f in Path("output").glob("*/errors.csv"))
df.groupby("experiment_id")["err_angle"].median()
```

`face_id` e' univoco solo dentro un esperimento: su dati concatenati la
chiave e' la coppia (`experiment_id`, `face_id`). `analysis/analisi.py`
lavora su un esperimento alla volta.

## Sottoinsiemi di mappe

`subsets.py TEMPLATE OUT_DIR [K]` scrive un `.conf` per ogni sottoinsieme di k >= 3 (o solo k = K) delle mappe
(le righe di `simulator/stage1/pacing_sites.csv`; con sei mappe:
20 + 15 + 6 + 1 = 42). Ogni file e'
`aniso_const.conf` con cambiati solo `experiment_id`, `maps` e
`output_csv`. L'ID e' `sub_` con gli indici uniti da trattino, es.
`sub_0-2-4`: niente virgole, perche' `experiment_id` finisce non quotato
nel CSV. `analysis/analisi.py --sottoinsiemi` etichetta i sottoinsiemi
con gli stessi indici, uniti da `+` (es. `0+2+4`).

Per rigenerare i 42 risultati e il confronto:

```bash
python experiments/subsets.py experiments/aniso_const.conf output/subsets
./estimator/build/Tesi0 output/subsets/confs/*.conf
python analysis/analisi.py --sottoinsiemi
```

Il confronto e' APPAIATO: gli errori di ogni sottoinsieme si calcolano
sulle facce stimabili sia in esso sia in `sub_0-1-2-3-4-5`; l'area non
stimabile e' quella del sottoinsieme. In `output/subsets/figure/`:
`area_non_stimabile.png`, `angolo.png`, `cv_l.png`, `rapporto.png` (q50 e
q90 del valore assoluto, un punto per sottoinsieme per k) e `numeri.csv`
con una riga per sottoinsieme. Tutto e' rigenerato da `run_all.sh`.

I `.conf` non sono in git: derivano in modo deterministico da
`aniso_const.conf` e da `pacing_sites.csv`, quindi la ricetta e' lo
script. Stanno in `output/subsets/confs/`, accanto ai risultati, e
`output/` e' ignorato. `sub_0-1-2-3-4-5` deve dare lo stesso
`errors.csv` di `aniso_const`, tolta `experiment_id`.

## Dimensione delle zone

Gli scenari `aniso_zone_N`, N = 43270, 10000, 2500, 600, 150, 40, seme 0,
hanno le 25 combinazioni di `aniso_random.csv` assegnate per zona (vedi
*Scenari*). `zone.py N ...` scrive un `.conf` per
ciascuno in `output/zone/confs/`: `aniso_const.conf` (maschera a 30
gradi) con cambiati solo `experiment_id` (`aniso_zone_N`), `lat_ply`,
`ground_truth_ply` e `output_csv` (`output/zone/aniso_zone_N/errors.csv`).
Come per i sottoinsiemi, i `.conf` non sono in git: la ricetta e' lo
script.

```bash
python experiments/zone.py 43270 10000 2500 600 150 40
./estimator/build/Tesi0 output/zone/confs/*.conf
python analysis/analisi.py --zone
```

`analisi.py --zone` calcola per ogni N il diametro medio di zona
d = sqrt(4 A / (pi N)) in mm, A = area totale della mesh, e scrive in
`output/zone/figure/`: `numeri.csv`, una riga per N con d, area non
stimabile, coppie escluse, q50 e q90 dei tre errori e mediane con segno
dei relativi; `zone.png`, quattro pannelli (area non stimabile, angolo,
CV_L, rapporto) in funzione di d in scala log, q50 e q90, con il valore di
`aniso_const` tratteggiato come riferimento.

## Confronto con Roney

Metodo 2 di Roney et al. 2019 su terne di mappe (`estimator/roney.py`,
`experiments/roney.conf`; comandi in fondo a `run_all.sh`):
`output/roney/SCENARIO/{roney,faccia}/sub_*/errors.csv`, confronto in
`output/roney/figure/` (`analisi.py --roney`).

## Perche' non JSON

Per leggere sette valori servirebbe vendorizzare un parser JSON. Questo
formato si legge in quindici righe senza dipendenze, ed e' coerente con
`region_cv/*.csv`, che e' gia' testo piatto. Una dipendenza in meno e' una
cosa in meno che puo' rompersi fra due anni.
