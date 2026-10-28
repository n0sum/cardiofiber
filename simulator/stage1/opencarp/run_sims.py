#!/usr/bin/env python3
"""
Una simulazione eikonale openCARP per ogni riga di pacing_sites.csv.

Il punto di stimolo e' il baricentro della cella cell_id sulla mesh in um,
arrotondato all'intero piu' vicino. I siti si definiscono per indice di
cella, non per coordinate (allo stadio 2 gli indici restano, le coordinate
cambiano). La riga k del CSV e' la mappa k.

Le CV vengono dal CSV dello scenario (region_tag,CV_L,CV_T in mm/ms), unica
fonte: qui si convertono in mm/s per openCARP, una gregion per tag. Il nome
dello scenario e' il nome del CSV.

Le run finiscono in data/reference/sims/SCENARIO/run_k/. openCARP non ha
un'opzione per sovrascrivere una run esistente, quindi la cartella si
cancella prima.

Uso, dalla radice del repository:
    python simulator/stage1/opencarp/run_sims.py experiments/region_cv/aniso_const.csv
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
HERE     = Path(__file__).resolve().parent
ROOT     = HERE.parents[2]                                    # radice del repository
SITES    = HERE.parent / "pacing_sites.csv"
MESH_UM  = ROOT / "data" / "reference" / "LA_um"              # senza estensione
SIMS_DIR = ROOT / "data" / "reference" / "sims"
# --------------------------------------------------------------------------


def gregion_args(cv_csv):
    """Argomenti gregion di openCARP dal CSV delle CV. Unico punto della
    conversione mm/ms -> mm/s."""
    rows = list(csv.DictReader(l for l in open(cv_csv) if not l.startswith("#")))
    args = ["-num_gregions", str(len(rows))]
    for i, row in enumerate(rows):
        vel_l = float(row["CV_L"]) * 1000                     # mm/ms -> mm/s
        vel_t = float(row["CV_T"]) * 1000
        args += [f"-gregion[{i}].num_IDs", "1",
                 f"-gregion[{i}].ID[0]", row["region_tag"],
                 f"-gregion[{i}].dream.vel_l", f"{vel_l:g}",
                 f"-gregion[{i}].dream.vel_t", f"{vel_t:g}",
                 f"-gregion[{i}].dream.vel_n", f"{vel_t:g}"]
    return args


def main():
    cv_csv = Path(sys.argv[1]).resolve()
    sims_dir = SIMS_DIR / cv_csv.stem
    gregions = gregion_args(cv_csv)
    pts = np.loadtxt(f"{MESH_UM}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{MESH_UM}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    sims_dir.mkdir(parents=True, exist_ok=True)

    for k, row in enumerate(csv.DictReader(open(SITES))):
        p = np.round(pts[tri[int(row["cell_id"])]].mean(0)).astype(np.int64)
        shutil.rmtree(sims_dir / f"run_{k}", ignore_errors=True)
        subprocess.run([HERE / "run_eik.sh", str(MESH_UM), str(k), *map(str, p), *gregions],
                       cwd=sims_dir, check=True)


if __name__ == "__main__":
    main()
