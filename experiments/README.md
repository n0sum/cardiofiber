# Esperimenti

Un file `.conf` per esperimento. L'eseguibile ne riceve uno o piu' e
scrive un CSV di errori per faccia per ciascuno:

```bash
./estimator/build/Tesi0 experiments/a.conf experiments/b.conf
```

Tutti i file vengono letti e validati prima di eseguire il primo.

Cosi' "anisotropia costante + anisotropia variabile + fmap in una sola
run" e' una lista di tre file, non tre versioni del programma.

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
| `ground_truth_ply`          | fibre e tag di regione per faccia              |
| `region_cv_csv`             | region_tag -> (CV_L, CV_T), in mm/ms           |
| `maps`                      | `all`, oppure indici separati da virgola       |
| `mask.min_singular_value`   | soglia sul valore singolare minimo del 3x3     |
| `mask.collision_angle_deg`  | soglia del criterio di collisione, gradi       |
| `output_csv`                | dove scrivere gli errori per faccia            |
| `output_mask_csv`           | dove scrivere le maschere per coppia           |

Tutte le chiavi sono obbligatorie. Le soglie `mask.*` non hanno un
default nel codice: esistono solo qui.

`maps = 0,2,4` tiene solo quelle mappe, nell'ordine dato; un indice
fuori range e' un errore.

Chiave sconosciuta (con `file:riga`) o mancante, file inesistente: errore,
nessun default silenzioso. Se una chiave compare piu' volte vale l'ultima.
La directory di `output_csv` viene creata se non esiste.

La directory di `output_mask_csv` idem.

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
nella LAA.

Una faccia con meno di 3 righe rimaste e' non stimabile.

`mask.collision_angle_deg = 181` non esclude nulla, senza codice speciale:
l'angolo non supera mai 180. `aniso_const_nomask.conf` la usa come
regressione: il suo `errors.csv`, tolte `experiment_id` e `k_used`, e'
identico byte per byte a quello di prima delle maschere.

Le LAT negative (FIM non convergito, -1) vengono rifiutate a monte, da
`simulator/convert/multiLatToPly.py`, prima di scrivere `lat_ply`.

I nomi delle mappe vengono dai commenti `lat_k NOME` dell'header di
`lat_ply`, che deve averne uno per ogni mappa. Come `experiment_id`
finiscono non quotati nel CSV delle maschere.

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
| `cond`          | numero di condizionamento; vuoto se non calcolato       |
| `CV_L_est`      | CV longitudinale stimata, mm/ms; vuota se non stimabile |
| `CV_T_est`      | CV trasversale stimata, mm/ms; vuota se non stimabile   |
| `CV_L_true`     | CV longitudinale vera, mm/ms                            |
| `CV_T_true`     | CV trasversale vera, mm/ms                              |
| `err_angle`     | errore angolare di fibra, gradi; vuoto se non stimabile |

Gli errori relativi sulle velocita' non stanno nel CSV: li calcola
`analysis/analisi.py` da `CV_*_est` e `CV_*_true`, con segno,
`(stimato - vero) / vero`.

`experiment_id` sta su ogni riga perche' i CSV di esperimenti diversi si
possano concatenare e raggruppare senza perdere la provenienza:

```python
df = pd.concat(pd.read_csv(f) for f in Path("output").glob("*/errors.csv"))
df.groupby("experiment_id")["err_angle"].median()
```

`face_id` e' univoco solo dentro un esperimento: su dati concatenati la
chiave e' la coppia (`experiment_id`, `face_id`). `analysis/analisi.py`
lavora su un esperimento alla volta e si ferma se ne trova piu' d'uno.

## CSV delle maschere

`output_mask_csv`, una riga per coppia (faccia, mappa) di ogni mappa
selezionata:

| colonna         | significato                                                   |
|-----------------|---------------------------------------------------------------|
| `experiment_id` | dal `.conf`                                                   |
| `face_id`       | indice della faccia, come in `errors.csv`                     |
| `map_name`      | nome del sito di stimolo della mappa                          |
| `masked`        | 0 tenuta, 1 esclusa per collisione                            |
| `r`             | `g^T M_true g - 1`, 6 cifre significative                     |

La mappa si identifica per nome, non per indice: `maps` rinumera le
mappe, e in un esperimento su un sottoinsieme l'indice indicherebbe la
mappa sbagliata senza dare errore.

`r` usa `M_true = CV_L^2 f f^T + CV_T^2 t t^T`, con `f` la fibra vera
proiettata sul piano della faccia e normalizzata, `t = n x f`, CV in mm/ms
da `region_cv_csv`, `g` il gradiente LAT della faccia in ms/mm. Vale 0
quando la coppia e' coerente con un fronte piano nel tessuto vero. Si
calcola per tutte le coppie, escluse comprese, solo nell'export: e' una
metrica di valutazione e non rientra mai nella stima.

## Perche' non JSON

Per leggere nove valori servirebbe vendorizzare un parser JSON. Questo
formato si legge in quindici righe senza dipendenze, ed e' coerente con
`region_cv/*.csv`, che e' gia' testo piatto. Una dipendenza in meno e' una
cosa in meno che puo' rompersi fra due anni.
