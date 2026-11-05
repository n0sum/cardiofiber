#!/usr/bin/env python3
"""
Una simulazione eikonale openCARP per ogni riga di pacing_sites.csv.

Il punto di stimolo e' il baricentro della cella cell_id sulla mesh in um,
arrotondato all'intero piu' vicino. I siti si definiscono per indice di
cella, non per coordinate. La riga k del CSV e' la mappa k.

Le CV vengono dal CSV dello scenario (region_tag,CV_L,CV_T in mm/ms), unica
fonte: qui si convertono in mm/s per openCARP, una gregion per tag. Dagli
stessi tag viene l'unico dominio eikonale (phys_region).

MESH_BASE e' la mesh in um senza estensione: i suoi tag per faccia devono
essere quelli del CSV.

Le run finiscono in data/reference/sims/NOME/run_k/. openCARP non ha
un'opzione per sovrascrivere una run esistente, quindi la cartella si
cancella prima.

Il FIM si ferma dopo MAX_ITER iterazioni; le run che convergono ne fanno
97-154. Dopo ogni run si leggono le iterazioni da eik_stats.dat: se ha
raggiunto MAX_ITER, la mappa non converge. Scrive
data/reference/sims/NOME/mappe.csv (mappa,iterazioni,converge): toPly.py
mette nel LAT PLY solo le mappe che convergono.

Uso, dalla radice del repository:
    python simulator/stage1/opencarp/run_sims.py NOME CSV MESH_BASE
    python simulator/stage1/opencarp/run_sims.py aniso_zone_50 experiments/region_cv/aniso_random.csv data/reference/aniso_zone_50/LA_um
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- percorsi
HERE     = Path(__file__).resolve().parent
ROOT     = HERE.parents[2]                                    # radice del repository
SITES    = HERE.parent / "pacing_sites.csv"
SIMS_DIR = ROOT / "data" / "reference" / "sims"
MAX_ITER = 1000                                               # iterazioni massime del FIM
# --------------------------------------------------------------------------


def region_args(cv_csv):
    """Argomenti phys_region e gregion di openCARP dal CSV delle CV. Unico
    punto della conversione mm/ms -> mm/s."""
    rows = list(csv.DictReader(l for l in open(cv_csv) if not l.startswith("#")))
    args = ["-num_phys_regions", "1",
            "-phys_region[0].name", "Eikonal domain",
            "-phys_region[0].ptype", "2",
            "-phys_region[0].num_IDs", str(len(rows))]
    for i, row in enumerate(rows):
        args += [f"-phys_region[0].ID[{i}]", row["region_tag"]]
    args += ["-num_gregions", str(len(rows))]
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
    name = sys.argv[1]
    cv_csv = Path(sys.argv[2]).resolve()
    mesh_um = Path(sys.argv[3]).resolve()
    sims_dir = SIMS_DIR / name
    regions = region_args(cv_csv)
    pts = np.loadtxt(f"{mesh_um}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{mesh_um}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    sims_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for k, row in enumerate(csv.DictReader(open(SITES))):
        p = np.round(pts[tri[int(row["cell_id"])]].mean(0)).astype(np.int64)
        shutil.rmtree(sims_dir / f"run_{k}", ignore_errors=True)
        subprocess.run([HERE / "run_eik.sh", str(mesh_um), str(k), *map(str, p),
                        "-dream.fim.max_iter", str(MAX_ITER), *regions],
                       cwd=sims_dir, check=True)
        # terza riga di eik_stats.dat: la quarta colonna sono le iterazioni del FIM
        its = int(open(sims_dir / f"run_{k}" / "eik_stats.dat").readlines()[2].split()[3])
        rows.append({"mappa": k, "iterazioni": its, "converge": int(its < MAX_ITER)})
    pd.DataFrame(rows).to_csv(sims_dir / "mappe.csv", index=False)


if __name__ == "__main__":
    main()
