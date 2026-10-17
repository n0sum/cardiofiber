# SETUP_SERVER.md — ambiente di build/esecuzione di Tesi0

Documentazione riproducibile dell'ambiente per compilare ed eseguire `Tesi0` su
**Ubuntu 26.04.1 LTS (resolute)**. Requisito di tesi: ogni pacchetto qui sotto ha
versione esatta e comando esatto usato per installarlo.

- **Data setup:** 2026-09-10
- **Macchina:** Ubuntu 26.04.1 LTS (resolute), x86_64, 32 core, GPU NVIDIA presente
  (`torch.cuda.is_available() == True`)
- **Shell:** fish
- **Toolchain già presente (non installata da noi):** GCC/G++ 15.2.0
  (`Ubuntu 15.2.0-16ubuntu1`), GNU Make 4.4.1, GNU ld 2.46, pkg-config 2.5.1,
  git, OpenMP 4.5 (via GCC / libgomp)

---

## 1. Pacchetti di sistema (apt)

Repo: `http://it.archive.ubuntu.com/ubuntu resolute`.

| Pacchetto | Versione installata | Comando |
|---|---|---|
| `cmake` | **4.2.3-2ubuntu2** | `sudo apt install -y cmake` |
| `ninja-build` | **1.13.2-1** | `sudo apt install -y ninja-build` |
| `libglm-dev` | **1.0.1+ds-3** | `sudo apt install -y libglm-dev` |

Comando unico eseguito (dall'utente, `sudo` interattivo):
```bash
sudo apt update
sudo apt install -y cmake ninja-build libglm-dev
```
`glmConfig.cmake` → `/usr/share/cmake/glm/glmConfig.cmake`; header → `/usr/include/glm/`.

`glm` è l'**unica** dipendenza C++ esterna (header-only, `find_package(glm REQUIRED)`
→ target `glm::glm`). `happly.h` (parser PLY) è già vendorizzato in `include/`.
**Non** eseguire `git submodule update --init` (gitlink orfano `include/happly`,
non serve e fallisce).

---

## 2. Miniforge (gestore ambienti Python, installazione utente)

| Voce | Valore |
|---|---|
| Installer | `Miniforge3-Linux-x86_64.sh` |
| URL | `https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` |
| Versione installer | Miniforge3 **26.7.2-0** (`constructor 3.16.1`, `MD5: 38c0a8c6a2297f7bf1d29e96da9688cd`) |
| SHA256 (file scaricato) | `281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05` |
| Prefix | `~/miniforge3` |
| `conda --version` | conda **26.7.2** |
| Python della base | 3.14.7 (non usato dal progetto — vedi env `tesi`) |

```bash
cd /tmp
curl -fL -o Miniforge3-Linux-x86_64.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
$HOME/miniforge3/bin/conda init fish bash
```

---

## 3. Ambiente conda `tesi` (Python 3.13)

**Python installato:** 3.13.15 (conda-forge).

### 3.1 Pacchetti core da conda-forge

```bash
conda create -y -n tesi -c conda-forge python=3.13 numpy scipy pyvista vtk matplotlib
```

| Pacchetto | Versione |
|---|---|
| python | 3.13.15 |
| numpy | 2.5.3 |
| scipy | 1.18.0 |
| pyvista | 0.48.4 |
| vtk | 9.6.2 |
| matplotlib | 3.11.1 |

### 3.2 Pacchetti da pip DENTRO l'env `tesi` (mai nel Python di sistema)

```bash
conda activate tesi
pip install pot robust_laplacian
pip install --force-reinstall --no-deps "geomstats @ git+https://github.com/geomstats/geomstats.git@main"
pip install "geomfum @ git+https://github.com/3diglab/geomfum.git@main"
# trame per il picking landmark headless — VERSIONI PINNATE (vedi DEVIAZIONE 3):
pip install "trame<4" "trame-server<4" "trame-client<4" "trame-vuetify<3" "trame-vtk==2.11.8"
```

> **DEVIAZIONE 1 — geomstats da git, non da PyPI.**
> `PROJECT_NOTES.md` / il messaggio d'errore in `lat_transfer.py` suggeriscono
> `pip install geomstats` (PyPI). La 2.8.0 di PyPI **non funziona con numpy ≥ 2**
> (`ImportError: cannot import name 'trapz' from 'numpy'`). La `main` di git
> (stessa version string 2.8.0, commit sotto) ha il fix `trapz → trapezoid`.
>
> **DEVIAZIONE 2 — geomfum upstream ora richiede `torch`.**
> `PROJECT_NOTES.md` diceva "installato senza torch da un fork locale patchato".
> Il fork/patch è andato perso; si è usato l'upstream `3diglab/geomfum@main`, che
> ha `torch` come dipendenza **core** (non opzionale) e tira dentro anche
> `pyfmaps`, `potpourri3d`, `gsops`, `meshio` e l'intero stack CUDA di torch
> (~3 GB). Installazione riuscita, GPU rilevata. Tutti gli import usati dagli
> script (`lat_transfer.py`, `diag_stages.py`, `diag_fmap.py`) sono verificati OK.
>
> **DEVIAZIONE 3 — trame va pinnato per pyvista 0.48.4.**
> Con `trame` 4.x (default) o `trame-vtk` < 2.11, il picking landmark headless di
> `landmarks.py` crasha con `ValueError: Protocol does not exist yet` in
> `plotter_ui()`. Combinazione funzionante trovata:
> `trame` 3.13.2 / `trame-server` 3.15.0 / `trame-client` 3.13.6 /
> `trame-vuetify` 2.7.2 / **`trame-vtk` 2.11.8** (pyvista 0.48.4 pin: `trame-vtk<2.11.9`).
> Verificato: il server browser del picker si avvia
> (`[landmarks] apri nel browser: http://127.0.0.1:PORT/`).

**Commit esatti installati:**

| Pacchetto | Sorgente / commit | Versione |
|---|---|---|
| geomfum | `git+https://github.com/3diglab/geomfum.git@63737102c23a701d03d7c6ae00985a38c904bdb7` | 1.0.0 |
| geomstats | `git+https://github.com/geomstats/geomstats.git@9941df0a37d900d540b7869d8c2898fca46266c8` | 2.8.0 |

**Versioni pip rilevanti (env `tesi`):**

| Pacchetto | Versione | Note |
|---|---|---|
| torch | 2.14.0+cu130 | dep core di geomfum; CUDA 13.0 |
| pot (POT) | 0.9.7.post1 | Sinkhorn |
| robust-laplacian | 1.1.0 | LBO |
| pyfmaps | 1.3.0 | dep core di geomfum |
| potpourri3d | 1.4.0 | dep core di geomfum |
| gsops | 0.1.0 | dep core di geomfum |
| meshio | 5.3.5 | dep core di geomfum |
| scikit-learn | 1.9.0 | dep di geomfum/geomstats |
| pandas | 3.0.5 | dep di geomstats |
| trame / trame-server / trame-client | 3.13.2 / 3.15.0 / 3.13.6 | picking landmark headless — pin (DEVIAZIONE 3) |
| trame-vtk / trame-vuetify | 2.11.8 / 2.7.2 | picking landmark headless — pin (DEVIAZIONE 3) |

Snapshot completi salvati fuori dal repo (scratchpad di sessione):
`pip_freeze_tesi.txt`, `conda_list_tesi.txt`, `environment_tesi.yml`.
Rigenerabili con `conda list -n tesi`, `pip freeze` (env attivo),
`conda env export -n tesi --no-builds`.

> Nota: `h5py` e `meshio` che avevi citato per il vecchio env "tesi" non servono
> agli script (PLY I/O scritto a mano); `meshio` è comunque entrato come dep di
> geomfum.

---

## 4. Build C++ (Release, da zero) — ✅ ESEGUITO

```bash
cd /home/gian/cardiofiber/Tesi0
rm -rf cmake-build-debug build          # nessun riuso di build vecchie
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel "$(nproc)"   # 32 core  (fish: (nproc))
```

- Generatore: **Ninja 1.13.2**
- Compilatore: **GNU C++ 15.2.0** (`/usr/bin/c++`)
- `find_package(glm REQUIRED)` → risolto via `libglm-dev` senza errori
- Build: **exit 0**, prodotto `build/Tesi0` (~258 KB)
- Un solo warning (pre-esistente, non toccato):
  `Mesh.cpp:184: warning: no return statement in function returning non-void` in
  `Mesh::populateHalfs()`

---

## 5. Esecuzione

- L'eseguibile `Tesi0` va lanciato con **cwd = `build/`** (path relativi `../` in
  `main.cpp`): `conda activate tesi; cd build; ./Tesi0`.
- L'env conda `tesi` deve essere **attivo** quando si lancia `Tesi0` (invoca
  `python3 ../fix_holes.py` e `python3 ../lat_transfer.py` dal `PATH`).
- `lat_transfer.py` richiede una cache landmark valida
  (`build/mesh1elaborata.landmarks.txt` + `.meta`, idem `mesh2elaborata`);
  su una build nuova non esiste → picking interattivo una tantum, da `build/`:
  ```bash
  conda activate tesi
  cd build
  python ../lat_transfer.py mesh1elaborata.ply mesh2elaborata.ply mesh1_dual_lat.ply --repick
  ```
  Senza `$DISPLAY` (sessione Wayland/tty) parte un server browser: aprire l'URL
  stampato (`http://127.0.0.1:PORT/`), click sinistro = aggiungi landmark, gli
  stessi punti anatomici nello stesso ordine sulle due mesh, poi "Salva". La
  cache resta in `build/` — **se si rifà `rm -rf build` va rifatto il picking**
  (conviene copiare i 4 file `*.landmarks.*` altrove come backup).
- Dopo il picking, `cd build && ./Tesi0` esegue la pipeline completa fino agli
  OBJ delle tangenti.
- Il picker trame logga in modo molto verboso (INFO per token) — è il default di
  `trame-server`, innocuo; per silenziarlo `export TRAME_LOG_LEVEL=WARNING` o
  `logging.getLogger('trame').setLevel(logging.WARNING)` (non fatto: sarebbe una
  modifica agli script).

---

## 6. Verifiche eseguite

| # | Test | Esito |
|---|---|---|
| 1 | Build C++ (`cmake --build build --parallel`) | ✅ exit 0 → `build/Tesi0` |
| 2 | `Tesi0` fase 1 da `build/` (env `tesi` attivo): import `.mesh` → `mesh1elaborata.ply` / `mesh2elaborata.ply` → `python3 ../fix_holes.py` (`[2] Buchi chiusi`) → `diag_lat_source.ply` | ✅ le 3 fasi OK; poi si ferma su `lat_transfer.py` per cache landmark mancante (atteso) |
| 3 | Smoke-test import: tutti i simboli `geomfum` usati da `lat_transfer.py`, `diag_stages.py`, `diag_fmap.py` | ✅ OK |
| 4 | `diag_fmap.py` end-to-end su `diag_out/tgt_hemi_x_xfer.ply` + `src_hemi_x.ply` (`--k 120 --use-k 30`): LBO (`robust_laplacian`) + HKS/WKS + solve `FunctionalMap` + ICP + `ZoomOut` + P2P + Sinkhorn (`POT`) + round-trip | ✅ exit 0 in ~2m10s; PNG e PLY prodotti |
| 4b | `lat_transfer.py --repick`: legge le mesh, costruisce la UI trame, avvia il server browser del picker (`http://127.0.0.1:PORT/`) | ✅ arriva al prompt del picker (l'interazione poi è manuale) |
| 5 | Pipeline completa via `Tesi0` (dopo picking landmark) | ⏳ **manuale** — richiede il picking interattivo di 2× N landmark nel browser |

Nota su #4: il referto segnala "flip di simmetria" / "descrittori deboli" — è la
*diagnosi sulla qualità della mappa* (problema di ricerca già noto nel progetto, e
qui su mesh-emisfero senza landmark), **non** un errore d'ambiente. Lo script
esce 0 ed esegue tutta la catena numerica.

> Con `--k 40` (valore usato in una prima prova) `ZoomOut` fallisce con
> `ValueError: Not enough eigenvectors` perché lo schedule di default
> (`--zo-nit 6 --zo-step 5`) chiede più autovettori di quanti calcolati: usare
> `--k` ≥ `use_k + zo_nit*zo_step`. Con i default (`--k 120`) non si presenta.

## 7. Log delle modifiche

| Data | Modifica |
|---|---|
| 2026-09-10 | Setup completo: Miniforge 26.7.2, env `tesi` (Python 3.13.15), stack Python + geomfum(+torch) + geomstats(git). apt `cmake`/`ninja-build`/`libglm-dev`. Build C++ Release OK (`build/Tesi0`). Verifiche #1–#4b superate. Pin trame per fixare il picker headless (DEVIAZIONE 3). Resta manuale solo il picking landmark (#5). |
