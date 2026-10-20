# Note sullo stadio 3 (functional maps)

Prototipi Python, oggi non usati dalla pipeline. Le dipendenze (geomfum,
geomstats, trame) non sono nell'ambiente: i pin esatti e i motivi stanno
nel blocco commentato in fondo a `environment.yml`.

## Stato al 2026-07-10

- La functional map tra le due mesh **flippa** (traccia della
  sottomatrice iniziale negativa, una sola coordinata correla). Ipotesi:
  i landmark usati (anello mitralico, ostii delle vene, auricola) sono
  quasi collineari lungo l'asse dell'atrio e lasciano libera una
  simmetria speculare.
- `DEFAULT_LABELS` in `landmarks.py` e' passato da 6 a 9 punti, con tre
  fuori asse: parete posteriore, tetto tra le vene, setto.
  `check_landmark_spread` avvisa se i landmark sono quasi collineari.
- Bug corretto: la cache dei landmark non veniva confrontata con il
  numero di label, e la diagnostica precedente girava in silenzio con 4
  punti. Ora una cache disallineata e' un errore che suggerisce `--repick`.
- **Da fare:** `--repick` interattivo su entrambe le mesh (9 punti
  ciascuna), poi rilanciare `diag_stages.py` per vedere se il flip
  sparisce, e verificare il peso `--n-times-landmarks` (default 3).
- `diag_fmap.py` sembra il predecessore di `diag_stages.py`: non usa i
  landmark e non ha ricevuto gli ultimi aggiornamenti.

## Note pratiche

- **Picking dei landmark senza X11** (Wayland): `landmarks.py` avvia un
  server trame e stampa un URL `http://127.0.0.1:PORT/`. Nel browser:
  click sinistro aggiunge un punto, stessi punti anatomici nello stesso
  ordine sulle due mesh, poi "Salva". La cache `<mesh>.landmarks.txt` +
  `.meta` sta accanto alla mesh: conviene copiarla altrove, perche'
  rifarla richiede di nuovo il picking a mano.
- trame logga a livello INFO per ogni token; per silenziarlo
  `TRAME_LOG_LEVEL=WARNING`.
- In `diag_fmap.py` e `diag_stages.py` ZoomOut chiede piu' autovettori di
  quelli usati all'inizio: serve `--k >= --use-k + --zo-nit * --zo-step`,
  altrimenti fallisce con `ValueError: Not enough eigenvectors`. Con i
  default (120 >= 30 + 6*5) va bene.
