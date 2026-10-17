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
si saltano. Spazi attorno a `=` e fine riga CRLF sono tollerati.

| chiave                      | significato                                   |
|-----------------------------|-----------------------------------------------|
| `experiment_id`             | identificativo, prima colonna del CSV; solo `[A-Za-z0-9_.-]` |
| `description`               | testo libero, finisce nei metadati del run     |
| `lat_ply`                   | mesh + K mappe LAT: l'ingresso dello stimatore |
| `ground_truth_ply`          | fibre e tag di regione per faccia              |
| `region_cv_csv`             | region_tag -> (CV_L, CV_T), in mm/ms           |
| `maps`                      | `all`, oppure indici separati da virgola       |
| `mask.min_singular_value`   | soglia sul valore singolare minimo del 3x3     |
| `mask.min_area_vector_norm` | soglia di faccia degenere                      |
| `output_csv`                | dove scrivere gli errori per faccia            |

Tutte le chiavi sono obbligatorie tranne `description`. Le soglie
`mask.*` non hanno un default nel codice: esistono solo qui.

`maps = 0,2,4` tiene solo quelle mappe, nell'ordine dato; un indice
ripetuto o fuori range e' un errore.

Chiave sconosciuta, ripetuta o mancante, valore vuoto o non numerico,
file inesistente: errore con `file:riga`, nessun default silenzioso.
La directory di `output_csv` viene creata se non esiste.

Tutti i percorsi sono relativi alla radice del repository, e l'eseguibile
si lancia da la'.

`experiment_id` finisce non quotato in ogni riga del CSV, per questo
ammette solo lettere, cifre, `_`, `-` e `.`: una virgola sposterebbe tutte
le colonne senza nessun errore.

## CSV di uscita

Una riga per faccia, nessuna aggregazione. Colonne, nell'ordine:

| colonna         | significato                                            |
|-----------------|--------------------------------------------------------|
| `experiment_id` | dal `.conf`, uguale su tutte le righe                   |
| `face_id`       | indice della faccia nella mesh                          |
| `area`          | area della faccia, mm^2                                 |
| `region_tag`    | tag di regione dalla ground truth                       |
| `estimable`     | 1 se la faccia supera le maschere, 0 altrimenti         |
| `cond`          | numero di condizionamento; vuoto se non calcolato       |
| `CV_L_est`      | CV longitudinale stimata, mm/ms; vuota se non stimabile |
| `CV_T_est`      | CV trasversale stimata, mm/ms; vuota se non stimabile   |
| `CV_L_true`     | CV longitudinale vera, mm/ms                            |
| `CV_T_true`     | CV trasversale vera, mm/ms                              |
| `err_angle`     | errore angolare di fibra, gradi                         |
| `err_CV_L_rel`  | errore relativo su CV_L                                 |
| `err_CV_T_rel`  | errore relativo su CV_T                                 |
| `err_ratio_rel` | errore relativo sul rapporto CV_L/CV_T                  |

Gli errori relativi sono con segno, `(stimato - vero) / vero`. Le ultime
quattro colonne sono vuote sulle facce non stimabili.

`experiment_id` sta su ogni riga perche' i CSV di esperimenti diversi si
possano concatenare e raggruppare senza perdere la provenienza:

```python
df = pd.concat(pd.read_csv(f) for f in Path("output").glob("*/errors.csv"))
df.groupby("experiment_id")["err_angle"].median()
```

`face_id` e' univoco solo dentro un esperimento: su dati concatenati la
chiave e' la coppia (`experiment_id`, `face_id`). `analysis/analisi.py`
lavora su un esperimento alla volta e si ferma se ne trova piu' d'uno.

## Perche' non JSON

Per leggere otto valori servirebbe vendorizzare un parser JSON. Questo
formato si legge in quindici righe senza dipendenze, ed e' coerente con
`region_cv/*.csv`, che e' gia' testo piatto. Una dipendenza in meno e' una
cosa in meno che puo' rompersi fra due anni.
