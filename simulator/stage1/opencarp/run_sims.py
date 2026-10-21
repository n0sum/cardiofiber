#!/usr/bin/env python3
"""
Una simulazione eikonale openCARP per ogni riga di pacing_sites.csv.

Il punto di stimolo e' il baricentro della cella cell_id sulla mesh in um,
arrotondato all'intero piu' vicino. I siti si definiscono per indice di
cella, non per coordinate: le colonne x_um, y_um, z_um del CSV non si
leggono (allo stadio 2 gli indici restano, le coordinate cambiano).

Le run finiscono in data/reference/sims/run_NOME/. openCARP non ha un'opzione
per sovrascrivere una run esistente, quindi la cartella si cancella prima.
"""

import csv
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
HERE     = Path(__file__).resolve().parent
ROOT     = HERE.parents[2]                                    # radice del repository
SITES    = HERE.parent / "pacing_sites.csv"
MESH_UM  = ROOT / "data" / "reference" / "LA_um"              # senza estensione
SIMS_DIR = ROOT / "data" / "reference" / "sims"
# --------------------------------------------------------------------------


def main():
    pts = np.loadtxt(f"{MESH_UM}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{MESH_UM}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    SIMS_DIR.mkdir(parents=True, exist_ok=True)

    for row in csv.DictReader(open(SITES)):
        name = row["name"]
        p = np.round(pts[tri[int(row["cell_id"])]].mean(0)).astype(np.int64)
        print(f"{name}: cella {row['cell_id']}  stimolo {p} um", flush=True)
        shutil.rmtree(SIMS_DIR / f"run_{name}", ignore_errors=True)
        t0 = time.time()
        subprocess.run([HERE / "run_eik.sh", str(MESH_UM), name, *map(str, p)],
                       cwd=SIMS_DIR, check=True)
        print(f"{name}: {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
