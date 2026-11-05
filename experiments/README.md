# Esperimenti

Un file `.conf` per esperimento. L'eseguibile ne riceve uno o piu' e
scrive un CSV di errori per faccia per ciascuno, in sequenza:

```bash
./estimator/build/Tesi0 output/aniso_zone_1/confs/sub_0-1-2.conf output/aniso_zone_1/confs/sub_0-1-3.conf
```

I `.conf` li genera `subsets.py`, con `write_conf`: le righe di
`modello.conf` (le soglie, uguali per tutti) piu' `experiment_id`,
`lat_ply`, `ground_truth_ply`, `maps` e `output_csv`. Non sono in git: la
ricetta e' lo script.

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

Tutte le 6 chiavi sono obbligatorie.

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
zona con la stessa probabilita'. Semi e N stanno solo in `run_all.sh`:
N = 1, 8, 70, 600, seme 0, cioe' diametro medio di zona d = sqrt(4 A / (pi N))
~ 123 / sqrt(N) = 123, 44, 15, 5.0 mm, passi uguali in scala log. Con N = 1
le CV sono uniformi; con N = 600 le zone sono grandi come il disco del metodo
di Roney (2 `patch.radius_mm` di `roney.conf`, 5 mm), che presuppone CV
uniformi al suo interno: zone piu' piccole misurerebbero solo che il
presupposto e' violato.

Per ogni scenario, dalla radice:

```bash
python simulator/stage1/opencarp/run_sims.py NOME CSV MESH_BASE   # -> data/reference/sims/NOME/
python simulator/convert/toPly.py NOME CSV MESH_BASE              # -> lat_ply e ground_truth_ply
```

`run_sims.py` lancia una run eikonale per sito di stimolo. `toPly.py` mette
in `lat_ply` le sei mappe come `lat_k`, con k l'indice del sito (la riga di
`simulator/stage1/pacing_sites.csv`); prende coordinate e fibre da
`data/augmenta/` (mm), triangoli e tag da `MESH_BASE.elem`.

`maps = 0,2,4` tiene solo le mappe di quei siti, nell'ordine dato; un sito
che non c'e' in `lat_ply` e' un errore.

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

Le LAT negative (FIM non convergito, -1) vengono rifiutate a monte, da
`simulator/convert/toPly.py`, prima di scrivere `lat_ply`.

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
df = pd.concat(pd.read_csv(f) for f in Path("output/aniso_zone_1").glob("*/sub_*/errors.csv"))
df.groupby("experiment_id")["err_angle"].median()
```

`face_id` e' univoco solo dentro un esperimento: su dati concatenati la
chiave e' la coppia (`experiment_id`, `face_id`). `experiment_id` non dice
lo scenario: lo dice la cartella, `output/SCENARIO/`. `analysis/analisi.py`
legge un `errors.csv` alla volta.

## La griglia di esperimenti

Per ogni scenario `aniso_zone_N` (la lista `ZONE` in `run_all.sh`):

```
output/
  aniso_zone_N/
    confs/sub_….conf             subsets.py
    faccia/sub_…/errors.csv      Tesi0, 42 sottoinsiemi
    roney/sub_i-j-k/errors.csv   estimator/roney.py, 20 terne
  figure/                        analysis/analisi.py
  mesh/                          Tesi0 --ply, da mappe.py
```

`run_all.sh` cancella `output/` prima di ricalcolarla, quindi nelle figure
non entrano risultati vecchi.

`subsets.py SCENARIO` scrive in `output/SCENARIO/confs/` un `.conf` per ogni
sottoinsieme di almeno 3 mappe (le righe di
`simulator/stage1/pacing_sites.csv`; con sei mappe 20 + 15 + 6 + 1 = 42),
con `output_csv` = `output/SCENARIO/faccia/ID/errors.csv`. L'ID e' `sub_`
con gli indici di sito uniti da trattino, es. `sub_0-2-4`: niente virgole,
perche' `experiment_id` finisce non quotato nel CSV. Le sei mappe sono
`sub_0-1-2-3-4-5`.

```bash
python experiments/subsets.py aniso_zone_1
./estimator/build/Tesi0 output/aniso_zone_1/confs/*.conf
```

`estimator/roney.py experiments/roney.conf SCENARIO ...` applica il metodo 2
di Roney et al. 2019 a tutte le terne di mappe di ogni scenario e scrive
`output/SCENARIO/roney/sub_i-j-k/errors.csv`, con le stesse colonne. Le
chiavi di `roney.conf` sono documentate nel file.

## Analisi

`python analysis/analisi.py`, senza argomenti, legge tutta la griglia e
scrive in `output/figure/`. Statistiche pesate per area, errori come in
*PLY degli errori*; q50, q90 e q99 sul valore assoluto, la mediana con
segno di cv_l e rapporto solo nelle tabelle. d = sqrt(4 A / (pi N)) e' il
diametro medio di zona in mm, A l'area totale della mesh (con N = 1,
d = 123 mm).

- `numeri.csv`: una riga per scenario e sottoinsieme dello stimatore per
  faccia: area non stimabile ed errori sulle facce stimabili del
  sottoinsieme.
- `confronto.csv`: una riga per scenario, terna e metodo (`faccia`,
  `roney`), con l'area non stimabile del metodo da solo, l'area appaiata
  (facce stimabili da entrambi i metodi con quella terna) e gli errori sulle
  facce appaiate. Le righe con `sub` = `tutte` mettono insieme le facce
  appaiate delle 20 terne dello scenario.

Quattro grafici, ognuno in quattro PNG: `_area` (area non stimabile, un
pannello) e `_angolo`, `_cv_l`, `_rapporto`. Nei grafici 1 e 2 questi
hanno un solo pannello con q50, q90 e q99 insieme (cerchio, quadrato,
triangolo, colori Okabe-Ito), solo punti senza linee, con l'asse
verticale in scala log ed etichette in numeri semplici (1, 2, 5, 10...);
nei grafici 3 e 4 tre pannelli q50, q90 e q99. Gli altri assi verticali,
compresi i file `_area`, sono lineari e partono da 0 (nel grafico 4 anche
quello orizzontale, con gli stessi limiti).

| file                       | contenuto                                              |
|----------------------------|--------------------------------------------------------|
| `grafico1_mappe_*`         | `aniso_zone_1`, in funzione di k: un punto per sottoinsieme e quantile, sfalsati dentro ogni k |
| `grafico2_zone_*`          | sei mappe, in funzione di d (log): un punto per scenario e quantile |
| `grafico3_zone_mappe_*`    | in funzione di d (log): per ogni k mediana dei sottoinsiemi e fascia dal minimo al massimo |
| `grafico4_roney_*`         | un punto per terna e scenario (un colore e un simbolo per scenario): Roney in orizzontale, per faccia in verticale, con y = x |

Mediana, minimo e massimo fra sottoinsiemi contano ogni sottoinsieme uno;
la pesatura per area vale dentro ogni sottoinsieme. Nei grafici 2 e 3 la
linea tratteggiata e' il diametro del disco di Roney, 2 `patch.radius_mm`.

## PLY degli errori

Per MeshLab, che colora le facce con un solo campo, `quality`:

```bash
./estimator/build/Tesi0 --ply MESH.ply ERRORS.csv PREFISSO
```

con `MESH.ply` un `lat_ply` della stessa mesh e `ERRORS.csv` un CSV di
Tesi0 o di `estimator/roney.py`. Scrive tre PLY binari little endian, uno
per errore, con le definizioni di `load()` in `analysis/analisi.py`:

| file                    | `quality` per faccia                         |
|-------------------------|----------------------------------------------|
| `PREFISSO_angolo.ply`   | errore angolare di fibra, gradi              |
| `PREFISSO_cv_l.ply`     | errore relativo su CV_L, %, con segno        |
| `PREFISSO_rapporto.ply` | errore relativo su r = CV_L / CV_T, %, con segno |

Ogni file ha solo i vertici `x y z` (mm) e tutte le facce, nell'ordine di
`face_id`, con `property float quality`; sulle facce non stimabili e' NaN.
I valori sono quelli del CSV arrotondati a float. Li scrive happly, che
aggiunge una riga `comment`.

`mappe.py` sceglie i casi e stampa per ognuno una riga
`MESH ERRORS_CSV PREFISSO`, che `run_all.sh` passa a `Tesi0 --ply`; crea
`output/mesh/` se serve. I PLY sono
`output/mesh/METODO_SCENARIO_SUB_{angolo,cv_l,rapporto}.ply`, 12 in tutto:

1. per faccia, sei mappe, `aniso_zone_1`;
2. per faccia, sei mappe, `aniso_zone_600`;
3. e 4. per faccia e Roney sulla terna di `aniso_zone_600` con l'area
   appaiata piu' grande (facce stimabili da entrambi i metodi, pesate per
   area, dai due `errors.csv` della terna).

Gli scenari sono costanti in cima a `mappe.py`.

## Perche' non JSON

Per leggere sei valori servirebbe vendorizzare un parser JSON. Questo
formato si legge in quindici righe senza dipendenze, ed e' coerente con
`region_cv/*.csv`, che e' gia' testo piatto. Una dipendenza in meno e' una
cosa in meno che puo' rompersi fra due anni.
