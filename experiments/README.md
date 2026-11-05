# Esperimenti

Un file `.conf` per esperimento. L'eseguibile ne riceve uno o piu' e
scrive un CSV di errori per faccia per ciascuno, in sequenza:

```bash
./estimator/build/Tesi0 output/zone/confs/a.conf output/zone/confs/b.conf
```

I `.conf` li generano `zone.py`, `subsets.py` e `mappe.py`, tutti con
`write_conf` di `subsets.py`: le righe di `modello.conf` (le soglie, uguali
per tutti) piu' `experiment_id`, `lat_ply`, `ground_truth_ply`, `maps`,
`output_csv` e `output_ply`. Non sono in git: la ricetta e' lo script.

## Perche' i parametri stanno qui e non nel codice

Le soglie determinano i risultati: sono dati, non
implementazione. Se stanno nel codice, un risultato non e' piu'
attribuibile a una configurazione, e cambiare una soglia richiede di
ricompilare e di perdere la tracciabilita' di cosa ha prodotto cosa.

## Formato

`chiave = valore`, una per riga. Righe vuote e righe che iniziano con `#`
si saltano. Spazi attorno a `=` sono tollerati.

| chiave                      | significato                                   |
|-----------------------------|-----------------------------------------------|
| `experiment_id`             | identificativo, prima colonna del CSV          |
| `lat_ply`                   | mesh + mappe LAT: l'ingresso dello stimatore   |
| `ground_truth_ply`          | fibre, CV vere e tag di regione per faccia     |
| `maps`                      | `all`, oppure indici di sito separati da virgola |
| `mask.min_singular_value`   | soglia sul valore singolare minimo del 3x3     |
| `output_csv`                | dove scrivere gli errori per faccia            |
| `output_ply`                | base dei PLY degli errori, oppure `none`       |

Tutte le 7 chiavi sono obbligatorie.

`ground_truth_ply` e' `data/reference/LA_GT_NOME.ply` dello scenario
che ha prodotto `lat_ply`: fibra, CV_L, CV_T e tag per faccia, letti cosi'
come sono. La soglia `mask.min_singular_value` non ha un default nel
codice: esiste solo in `modello.conf`.

## Scenari

Uno scenario e' `aniso_zone_N`: la mesh in um divisa in N zone, con
`MESH_BASE` = `data/reference/aniso_zone_N/LA_um`. Tutti usano
`region_cv/aniso_random.csv` (`region_tag,CV_L,CV_T`, mm/ms): 25 tag, tutte
le combinazioni di 5 livelli di CV_L = 1.05 x (0.8, 0.9, 1.0, 1.1, 1.2) e 5
di r = CV_L/CV_T = 2.0 x (0.8, 0.9, 1.0, 1.1, 1.2), CV_L esterno e r
interno. La mesh la scrive `simulator/stage1/opencarp/random_tags.py
aniso_zone_N experiments/region_cv/aniso_random.csv 0 N`: N zone attorno a
N facce scelte a caso (faccia al centro piu' vicino), un tag a caso per
zona con la stessa probabilita'. Con N = 1 le CV sono uniformi, con
N = 43270 ogni faccia e' una zona. Semi e N stanno solo in `run_all.sh`:
N = 1, 50, 500, 5000, 25000, 43270, seme 0.

Per ogni scenario, dalla radice:

```bash
python simulator/stage1/opencarp/run_sims.py NOME CSV MESH_BASE   # -> data/reference/sims/NOME/
python simulator/convert/toPly.py NOME CSV MESH_BASE              # -> lat_ply e ground_truth_ply
```

`run_sims.py` lancia una run eikonale per sito di stimolo, con al massimo
1000 iterazioni del FIM (le run che convergono ne fanno un centinaio), e
legge le iterazioni da `eik_stats.dat`. In `sims/NOME/mappe.csv` scrive,
per sito, `mappa,iterazioni,converge`: una mappa che ha raggiunto 1000
iterazioni non converge. `toPly.py` mette in `lat_ply` solo le mappe che
convergono, come `lat_k` con k l'indice del sito (la riga di
`simulator/stage1/pacing_sites.csv`): gli indici possono saltare. Prende
coordinate e fibre da `data/augmenta/` (mm), triangoli e tag da
`MESH_BASE.elem`.

Tutto cio' che enumera le mappe usa solo quelle di `lat_ply`: `maps = all`,
i sottoinsiemi di `subsets.py`, le terne di `roney.py` e quindi la scelta
delle terne di `mappe.py`. `maps = 0,2,4` tiene solo le mappe di quei
siti, nell'ordine dato; un sito che non c'e' in `lat_ply` e' un errore.

Una chiave mancante, o un file inesistente, ferma Tesi0; non ci sono
default. Le chiavi sconosciute si ignorano; se una chiave compare piu'
volte vale l'ultima. La directory di `output_csv` viene creata se non
esiste.

Tutti i percorsi sono relativi alla radice del repository, e l'eseguibile
si lancia da la'.

`experiment_id` finisce non quotato in ogni riga del CSV: non deve
contenere virgole.

## Facce non stimabili

Ogni faccia usa tutte le mappe di `maps`, una riga per mappa nel suo sistema lineare.
E' non stimabile se ha meno di 3 mappe, se il valore singolare minimo del
sistema non supera `mask.min_singular_value` (fronti quasi paralleli) o se
l'autovalore minore del tensore stimato (CV_T^2) non e' positivo.

Le mappe che non convergono non entrano in `lat_ply` (vedi *Scenari*); una
LAT negativa (-1, FIM non convergito) in una mappa che converge ferma
`simulator/convert/toPly.py`.

## CSV di uscita

Una riga per faccia, nessuna aggregazione. Colonne, nell'ordine:

| colonna         | significato                                            |
|-----------------|--------------------------------------------------------|
| `experiment_id` | dal `.conf`, uguale su tutte le righe                   |
| `face_id`       | indice della faccia nella mesh                          |
| `area`          | area della faccia, mm^2                                 |
| `region_tag`    | tag di regione dalla ground truth                       |
| `estimable`     | 1 se la faccia e' stimabile, 0 altrimenti               |
| `k_used`        | mappe nel sistema della faccia: tutte quelle di `maps`  |
| `cond`          | numero di condizionamento (vedi sotto); vuoto se non calcolato |
| `CV_L_est`      | CV longitudinale stimata, mm/ms; vuota se non stimabile |
| `CV_T_est`      | CV trasversale stimata, mm/ms; vuota se non stimabile   |
| `CV_L_true`     | CV longitudinale vera, mm/ms                            |
| `CV_T_true`     | CV trasversale vera, mm/ms                              |
| `err_angle`     | errore angolare di fibra, `acos(|a.b|)` in gradi; vuoto se non stimabile |

`cond` e' `sigma_max / sigma_min` del sistema della faccia, una riga per
mappa: `(g1^2, sqrt(2) g1 g2, g2^2)`, con `(g1, g2)` il gradiente
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

## PLY degli errori

Per MeshLab, che colora le facce con un solo campo, `quality`. Se
`output_ply` = `BASE` non e' `none`, dopo il CSV Tesi0 scrive tre PLY
binari little endian, uno per errore, con le definizioni di `load()` in
`analysis/analisi.py`:

| file                | `quality` per faccia                         |
|---------------------|----------------------------------------------|
| `BASE_angolo.ply`   | errore angolare di fibra, gradi              |
| `BASE_cv_l.ply`     | errore relativo su CV_L, %, con segno        |
| `BASE_rapporto.ply` | errore relativo su r = CV_L / CV_T, %, con segno |

Ogni file ha solo i vertici `x y z` (mm) e tutte le facce, nell'ordine di
`face_id`, con `property float quality`; sulle facce non stimabili e' NaN.
I PLY si calcolano rileggendo il CSV appena scritto, quindi i valori sono
quelli del CSV arrotondati a float. Li scrive happly, che aggiunge una riga
`comment`.

Per un `errors.csv` scritto da altri (es. `estimator/roney.py`):

```bash
./estimator/build/Tesi0 --ply MESH.ply ERRORS.csv BASE
```

con `MESH.ply` un `lat_ply` della stessa mesh. Nei `.conf` di
`subsets.py` e `zone.py` `output_ply` e' `none`; i PLY li chiede solo
`mappe.py` (vedi *Mappe degli errori*).

## Sottoinsiemi di mappe

`subsets.py SCENARIO OUT_DIR [K]` scrive un `.conf` per ogni sottoinsieme di
k >= 3 (o solo k = K) delle mappe di SCENARIO che convergono (da
`data/reference/sims/SCENARIO/mappe.csv`; con sei mappe 20 + 15 + 6 + 1 =
42, con cinque 10 + 5 + 1 = 16), con `output_ply` = `none`. L'ID e' `sub_`
con gli indici di sito uniti da trattino, es. `sub_0-2-4`: niente virgole,
perche' `experiment_id` finisce non quotato nel CSV.
`analysis/analisi.py --sottoinsiemi` etichetta i sottoinsiemi con gli
stessi indici, uniti da `+` (es. `0+2+4`).

Sottoinsiemi su `aniso_zone_1` e confronto:

```bash
python experiments/subsets.py aniso_zone_1 output/subsets
./estimator/build/Tesi0 output/subsets/confs/*.conf
python analysis/analisi.py --sottoinsiemi
```

Il confronto e' APPAIATO: gli errori di ogni sottoinsieme si calcolano
sulle facce stimabili sia in esso sia nel sottoinsieme con tutte le mappe;
l'area non stimabile e' quella del sottoinsieme. In `output/subsets/figure/`:
`area_non_stimabile.png`, `angolo.png`, `cv_l.png`, `rapporto.png` (q50 e
q90 del valore assoluto, un punto per sottoinsieme per k) e `numeri.csv`
con una riga per sottoinsieme. Il sottoinsieme con tutte le mappe deve dare
lo stesso `errors.csv` di `output/zone/aniso_zone_1/`, tolta `experiment_id`.

## Dimensione delle zone

`zone.py N ...` scrive in `output/zone/confs/` un `.conf` per scenario
`aniso_zone_N`, con `maps` = `all` e `output_csv` =
`output/zone/aniso_zone_N/errors.csv`.

```bash
python experiments/zone.py 1 50 500 5000 25000 43270
./estimator/build/Tesi0 output/zone/confs/*.conf
python analysis/analisi.py zone/aniso_zone_1
python analysis/analisi.py --zone
```

`analisi.py zone/aniso_zone_1` scrive `numeri.csv` ed `errori.png` in
`output/zone/aniso_zone_1/figure/`. `analisi.py --zone` calcola per ogni N
il diametro medio di zona d = sqrt(4 A / (pi N)) in mm, A = area totale
della mesh (con N = 1, d = 123 mm), e scrive in `output/zone/figure/`:
`numeri.csv`, una riga per N con d, numero di mappe usate, area non
stimabile, q50 e q90 dei tre errori e mediane con segno dei relativi;
`zone.png`, quattro pannelli (area non stimabile, angolo, CV_L, rapporto)
in funzione di d in scala log, q50 e q90.

## Confronto con Roney

Metodo 2 di Roney et al. 2019 su terne di mappe (`estimator/roney.py`,
`experiments/roney.conf`; comandi in fondo a `run_all.sh`), per ogni
scenario sulle terne delle mappe che convergono:
`output/roney/SCENARIO/{roney,faccia}/sub_*/errors.csv`, confronto in
`output/roney/figure/` (`analisi.py --roney`): `numeri.csv` per scenario e
braccio, con il numero di mappe dello scenario; `terne.csv` per terna;
`zona_1.png`, un punto per terna su `aniso_zone_1`; `zone.png`, in
funzione di d.

## Mappe degli errori

`mappe.py` sceglie, per braccio, le terne di `aniso_zone_1` con angolo q50
minimo e massimo in `output/roney/figure/terne.csv` (facce stimabili nei
due bracci, i punti di `zona_1.png`) e scrive in `output/mappe/confs/` i
`.conf` dello stimatore per faccia: tutte le mappe e le sue due terne, con
`output_ply` in `output/mappe/`. Per le due terne di Roney stampa
`ERRORS_CSV BASE`, da passare a `Tesi0 --ply`. Comandi in fondo a
`run_all.sh`. Uscite, tre PLY per caso (vedi *PLY degli errori*):
`output/mappe/CASO_{angolo,cv_l,rapporto}.ply`, con CASO
`faccia_tutte_le_mappe` e `{faccia,roney}_{migliore,peggiore}_sub_i-j-k`.

## Perche' non JSON

Per leggere sei valori servirebbe vendorizzare un parser JSON. Questo
formato si legge in quindici righe senza dipendenze, ed e' coerente con
`region_cv/*.csv`, che e' gia' testo piatto. Una dipendenza in meno e' una
cosa in meno che puo' rompersi fra due anni.
