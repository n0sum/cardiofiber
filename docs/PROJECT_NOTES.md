# Note di progetto (Tesi0) — stato al 2026-07-10

Questo file è uno snapshot di quello che so sul progetto in questo momento, per dare
contesto a chiunque (io in una sessione futura, o altri) riprenda il lavoro. Non è
generato/aggiornato automaticamente: potrebbe invecchiare, verificare contro il
codice prima di fidarsi ciecamente di un dettaglio qui sotto.

## Obiettivo della tesi

Due mesh dell'atrio sinistro dello stesso paziente, una a **ritmo normale** e una in
**fibrillazione atriale (FA)**, ciascuna con un segnale scalare **LAT** (Local
Activation Time — tempo di attivazione del fronte elettrico) per vertice. Obiettivo
finale: mappare la LAT di una mesh sull'altra (via corrispondenza geometrica), così da
avere **due LAT per vertice** sulla stessa mesh → calcolarne i **gradienti** (direzione
di propagazione) → derivare la **conduzione anisotropica** (ellissi lungo le fibre
cardiache, i cui assi rappresentano il modulo dell'anisotropia nelle due condizioni).
Il metodo di corrispondenza scelto è la **functional map** (spettrale, LBO), rinforzata
con **landmark anatomici** perché i soli descrittori spettrali (HKS/WKS) non bastano a
disambiguare la quasi-simmetria della geometria atriale.

La parte di **ellissi/anisotropia non è ancora implementata**: al momento la pipeline
C++ arriva fino al calcolo dei gradienti (tangenti) del LAT, non oltre.

## Pipeline end-to-end (chi chiama cosa)

Orchestrata da **`main.cpp`** (C++, usa `Mesh`/`Face`/`Half` in `Mesh.cpp`/`Mesh.h`),
compilato come eseguibile `Tesi0` (CMake, C++20, dipende solo da `glm`). Gira da
`cmake-build-debug/` come working directory. Fasi:

1. **Import mesh grezze**: `Mesh::importMesh()` legge `.mesh` da `assets/rawMesh/`
   (`2-LA.mesh` = normale, `2-LA-FA.mesh` = FA — gli altri due file in quella cartella,
   `POST_MESH.mesh`/`PRE_MESH.mesh`, non sono referenziati da nessuna parte, sembrano
   residui). Calcola le normali (`updateVertexNormals`) ed esporta in PLY ASCII:
   `mesh1elaborata.ply` (normale) e `mesh2elaborata.ply` (FA), dentro
   `cmake-build-debug/`.
2. **Chiusura buchi**: `main.cpp` invoca `python3 ../fix_holes.py mesh1elaborata.ply
   mesh2elaborata.ply` (sovrascrive i PLY in place). Chiude solo i loop di bordo
   piccoli con fan-triangolazione, **non aggiunge vertici nuovi**.
3. **Functional map + transfer LAT**: `main.cpp` invoca `python3 ../lat_transfer.py
   mesh1elaborata.ply mesh2elaborata.ply mesh1_dual_lat.ply`. Nessun flag passato:
   usa i default hardcoded nello script. Produce `mesh1_dual_lat.ply` con le property
   `quality` (LAT originale di mesh1/FA) e `quality_transferred` (LAT di mesh2/normale
   trasferito via pullback p2p), più `oor_flag` (fuori range). Scrive anche un PLY
   diagnostico `*_oor.ply` con l'overshoot con segno.
4. **Gradienti**: `main.cpp` reimporta `mesh1_dual_lat.ply`
   (`Mesh::importPlyDualLat`), ricalcola le normali, ed esporta due PLY diagnostici
   separati (`diag_lat_originale.ply`, `diag_lat_trasferito.ply`, uno per canale LAT,
   utili per confronto visivo in MeshLab). Poi calcola i **gradienti tangenziali** del
   campo LAT (`Mesh::updateTangets1AsQuality1Gradient` per `quality`,
   `updateTangets2AsQuality2Gradient` per `quality_transferred`) ed esporta come OBJ
   (`tangenti_lat1.obj`, `tangenti_lat1_NotPiInv.obj`, `tangenti_lat2.obj`) — vettori
   tangenti alla superficie, proiezione del gradiente. `Mesh.h` ha già le funzioni di
   base (`facesGradients`, `distributeGradientOnVertex`, con flag `piInvariant` e
   `projectOnTangentPlane`) — quindi lo stadio "gradienti" esiste già in C++.
   **Le ellissi di anisotropia (assi = modulo dell'anisotropia) non esistono ancora**
   da nessuna parte nel codice (nessun grep positivo per "ellip"/"anisotrop"/"fiber").

## Script Python (tutti a livello di repo root, eseguiti con cwd
`cmake-build-debug/` quando lanciati da `main.cpp`, quindi i path relativi tipo
`../fix_holes.py` sono relativi a quella cwd)

- **`fix_holes.py`** (497 righe): chiude loop di bordo piccoli nei PLY, fan
  triangulation, nessun vertice nuovo. `python3 fix_holes.py mesh1.ply mesh2.ply`,
  sovrascrive in place.
- **`landmarks.py`** (444 righe): picking interattivo di landmark anatomici (pyvista),
  condiviso da `lat_transfer.py` e `diag_stages.py`. Vedi sezione dedicata sotto.
- **`lat_transfer.py`** (265 righe): script di **produzione**, quello che
  effettivamente fa il transfer LAT usato da `main.cpp`. Legge le due mesh, acquisisce
  i landmark (`get_matched_landmarks`), costruisce `TriangleMesh` geomfum, calcola
  spettro LBO (K=120 default), descrittori HKS+WKS+LandmarkHKS+LandmarkWKS, risolve la
  functional map (`FunctionalMap` con `SDPFactorBuilder`+`LBCFactorBuilder`+
  `MultFactorBuilder`), rileva il flip (`is_flipped`: traccia della sottomatrice
  leading < 0) e se serve ricalcola con termine orientation-preserving
  (`compute_fmap_orient`, L-BFGS-B esplicito), poi ICP (opzionale) + ZoomOut, poi p2p
  via Sinkhorn (fallback nearest-neighbour), poi pullback del LAT (`lat2[p2p]`, nessuna
  ricostruzione spettrale — non può "overshoot-are" oltre il range sorgente per
  costruzione del pullback, ma lo si può comunque ottenere se il p2p mappa a un
  sottoinsieme non rappresentativo, da cui i diagnostici di overshoot). Scrive
  `mesh1_dual_lat.ply` + un PLY diagnostico di overshoot.
- **`diag_stages.py`** (416 righe): script **diagnostico**, non tocca file di
  produzione. Confronta più varianti di pipeline fmap (stessa convenzione
  source/target di `lat_transfer.py`) su: `tr` (traccia, segno = flip), diagonalità,
  correlazione coordinate xyz trasferite vs native, correlazione autovettori LBO 1-3
  trasferiti vs nativi. Scrive `diag_out/report_stages.txt`. Vedi sezione "Stato
  landmark/fmap" sotto per le varianti attuali (a)-(g) e la logica di scelta finale.
- **`diag_fmap.py`** (535 righe): script diagnostico **più vecchio/diverso** da
  `diag_stages.py` — non usa ancora i landmark, produce PNG della fmap
  (matplotlib, opzionale) e vari PLY con funzioni geometriche trasferite (coordinate,
  autovettori LBO, distanza da un riferimento) per ispezione visiva in MeshLab. Scrive
  `diag_out/report.txt` (nome diverso da `report_stages.txt`). Sembra il predecessore
  di `diag_stages.py`: non risulta più invocato da `main.cpp`.

Tutti gli script Python (tranne `fix_holes.py`) importano `geomfum` — installato senza
torch da un fork locale patchato (`pipeline.py`+`convert.py`), richiede `geomstats`;
vedi memoria `project_geomfum_setup` per i dettagli d'installazione (fuori scope di
questo file, che riguarda solo lo stato del progetto/pipeline).

## Convenzione naming (identica in tutti gli script Python)

```
mesh1elaborata.ply = 2-LA.mesh    = NORMALE = TARGET = b
mesh2elaborata.ply = 2-LA-FA.mesh = FA      = SOURCE = a
mesh_src = TriangleMesh(V2, F2)   (a)
mesh_tgt = TriangleMesh(V1, F1)   (b)
p2p = converter(fmap, basis_src, basis_tgt)   shape=[n_tgt], valori in [0, n_src)
lat_transferred = lat_src[p2p]
```
Cioè: il LAT viaggia da mesh2 (FA, sorgente) verso mesh1 (normale, target), che è un
po' controintuitivo dato che concettualmente si vorrebbe "vedere come si comporterebbe
il LAT normale sulla geometria FA" — ma è la convenzione stabilita e usata ovunque,
compreso `main.cpp`.

## Stato landmark / qualità fmap (lavoro più recente, 2026-07-10)

Diagnosi (basata su `diag_out/report_stages.txt`, generato però — si è scoperto dopo —
con solo **4 landmark reali** invece dei 6 nominali per via di un bug, vedi sotto):
l'init grezzo della fmap **flippa** (`tr=-1.651`) e il termine orientation-preserving
non lo corregge a nessun peso testato (0.5/1/5). Nel flip un solo asse correla forte
(`corr_y≈+0.56`), gli altri due ~0: sintomo di landmark che vincolano una sola
direzione. Ipotesi: i landmark usati (anulus mitralico, ostii delle 4 vene polmonari,
auricola) sono quasi collineari lungo l'asse dominante dell'atrio, lasciando libertà
rotazionale/speculare di cui la fmap approfitta.

Interventi fatti oggi in `landmarks.py`:
- **`DEFAULT_LABELS`** esteso da 6 a **9 punti**, aggiungendo 3 punti fuori asse:
  *parete posteriore*, *tetto (roof) tra le vene*, *setto interatriale*.
- **`check_landmark_spread(V, indices, mesh_path, ratio_threshold=0.1)`**: nuova
  guardia di collinearità, PCA sulle posizioni 3D dei landmark, stampa un avviso
  esplicito se il rapporto tra il 2°/3° autovalore e il 1° è sotto soglia (landmark
  quasi collineari, atteso flip). Chiamata automaticamente da `get_matched_landmarks`
  per entrambe le mesh, prima della fmap. È solo un avviso (non blocca l'esecuzione).
- **Bug corretto**: `get_landmarks_for_mesh()` non verificava mai che la cache su
  disco (`<mesh>.landmarks.txt`) avesse lo stesso numero di punti delle label
  correnti — le cache reali (`cmake-build-debug/mesh1elaborata.landmarks.txt`,
  `mesh2elaborata.landmarks.txt`) avevano **solo 4 righe** anche quando
  `DEFAULT_LABELS` ne aveva già 6, quindi la diagnostica precedente girava
  silenziosamente su 4 landmark insufficienti. Ora solleva `RuntimeError` esplicito se
  `len(cache) != len(labels)`, consigliando `--repick`.

Interventi fatti oggi in `diag_stages.py`:
- Nuova variante **`(g) landmark+ZO(step=1)`**: solo ZoomOut (step=1, nit=30)
  dall'init grezzo (già landmark-aware visto che i descrittori landmark sono sempre
  attivi salvo `--no-landmarks`), senza orientation-preserving né ICP — pensata per
  essere confrontabile "a parità di step" con `(e) reinforced+ZO(step=1)`, isolando
  l'effetto del solo rinforzo orientation-preserving.
- La **selezione finale** (PASSO 4) ora considera tutte le varianti `(a)`-`(g)` +
  init rinforzato, non solo `(e)`/`(f)` come prima: sceglie quella con `tr>0` che
  massimizza il punteggio aggregato (`score()`, somma dei valori assoluti delle
  correlazioni), con un tie-break verso la pipeline più semplice (meno passaggi) se i
  punteggi sono vicini (entro 0.02).

**Ancora da fare, non fattibile da un agente**: le cache dei landmark vanno
rigenerate interattivamente dall'utente (`--repick` su entrambe le mesh, servono 9
click guidati per mesh — un agente non ha accesso visivo per identificare le
strutture anatomiche). Finché non viene fatto, `diag_stages.py`/`lat_transfer.py`
falliscono subito con l'errore di cache disallineata (verificato: comportamento
corretto, non un bug). Dopo il repicking, va rilanciato `diag_stages.py` per
verificare se il flip si risolve con i landmark ben distribuiti, e se serve
correggere il peso dei descrittori landmark (`--n-times-landmarks`, default 3).

## Vincoli ambiente (perché il picking a volte serve un fallback browser)

Desktop dell'utente: COSMIC (System76), **Wayland puro, nessun X11/XWayland
raggiungibile**. Le wheel pip di pyvista/VTK renderizzano interattivamente solo via
X11/GLX. `landmarks.py` gestisce questo con due livelli: finestra nativa se
`DISPLAY` è settato (`_has_display()`), altrimenti fallback a un server locale via
`trame`/`trame-vtk`/`trame-vuetify` (`pick_landmarks_browser`, rendering server-side
off-screen, l'utente interagisce da un tab browser qualsiasi). Se nessuno dei due è
disponibile e manca la cache, l'errore è esplicito con istruzioni.

## File di dati / output rilevanti

- `assets/rawMesh/2-LA.mesh`, `2-LA-FA.mesh` — mesh grezze di partenza (formato
  `.mesh`, non PLY). `POST_MESH.mesh`/`PRE_MESH.mesh` nella stessa cartella: non usati.
- `assets/landmarks/` — esiste ma non ancora ispezionato in questa sessione.
- `cmake-build-debug/mesh1elaborata.ply`, `mesh2elaborata.ply` — PLY post-import
  C++ e post fix_holes, working copy usata da tutti gli script Python.
- `cmake-build-debug/mesh1elaborata.landmarks.txt`, `mesh2elaborata.landmarks.txt` —
  cache landmark, **attualmente stale (4 punti, servono 9)**, da rigenerare.
- `cmake-build-debug/mesh1_dual_lat.ply` — output di `lat_transfer.py`, LAT originale +
  trasferito.
- `cmake-build-debug/diag_lat_originale.ply`, `diag_lat_trasferito.ply`,
  `diag_lat_source.ply` — PLY diagnostici per confronto visivo in MeshLab.
- `diag_out/report_stages.txt` — referto di `diag_stages.py` (**da rigenerare** con i
  9 landmark, quello attuale riflette la baseline a 4 punti/bug cache).
- `diag_out/report.txt` — referto del più vecchio `diag_fmap.py`.
- `diag_out_landmarks/` — cartella vuota, presumibilmente predisposta per un run
  futuro con `--outdir` diverso, non ancora popolata.

## Cose che NON so / da verificare

- Non ho ispezionato `assets/landmarks/` né `fix_holes.py`/`diag_fmap.py` in
  dettaglio (solo l'header/docstring) — se serve lavorarci, vanno letti per intero.
- Non è chiaro se `diag_fmap.py` sia ancora mantenuto o solo un residuo storico
  pre-landmark (non riceve gli aggiornamenti fatti a `diag_stages.py`/`landmarks.py`).
- La correlazione tra "gradienti già calcolati in C++" e il piano futuro di
  "anisotropia/ellissi" menzionato come obiettivo di tesi non ha ancora nessun
  codice: è la prossima fase logica ma non è stata scoping-ata né iniziata.
