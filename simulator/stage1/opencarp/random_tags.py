#!/usr/bin/env python3
"""
Tag casuali per zona per uno scenario.

Legge data/reference/LA_um.{pts,elem,lon}, sceglie N facce a caso come centri
e assegna ogni faccia al centro piu' vicino (distanza euclidea fra baricentri):
sono le N zone. A ogni zona va uno dei tag di CSV, a caso con la stessa
probabilita'. Un solo generatore, numpy.random.default_rng(SEED), per centri e
tag. N = numero di facce da' una zona per faccia.

Scrive data/reference/NOME/LA_um: .pts e .lon copiati, .elem con gli stessi
triangoli e il nuovo tag.

Uso, dalla radice del repository:
    python simulator/stage1/opencarp/random_tags.py NOME CSV SEED N
    python simulator/stage1/opencarp/random_tags.py aniso_random experiments/region_cv/aniso_random.csv 0 43270
"""

import csv
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

# ---------------------------------------------------------------- percorsi
ROOT   = Path(__file__).resolve().parents[3]                     # radice del repository
SRC    = ROOT / "data" / "reference" / "LA_um"                   # senza estensione
# --------------------------------------------------------------------------


def main():
    name, cv_csv, seed, n = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    dst = ROOT / "data" / "reference" / name / "LA_um"
    tags = [int(r["region_tag"]) for r in
            csv.DictReader(l for l in open(cv_csv) if not l.startswith("#"))]
    pts = np.loadtxt(f"{SRC}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{SRC}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    bary = pts[tri].mean(1)

    rng = np.random.default_rng(seed)
    centers = rng.choice(len(tri), size=n, replace=False)
    zone = cKDTree(bary[centers]).query(bary)[1]
    tag = rng.choice(tags, size=n)[zone]

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(f"{SRC}.pts", f"{dst}.pts")
    shutil.copyfile(f"{SRC}.lon", f"{dst}.lon")
    with open(f"{dst}.elem", "w") as f:
        f.write(f"{len(tri)}\n")
        for (a, b, c), t in zip(tri, tag):
            f.write(f"Tr {a} {b} {c} {t}\n")


if __name__ == "__main__":
    main()
