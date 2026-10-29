#!/usr/bin/env python3
"""
Tag casuali per faccia per uno scenario.

Legge data/reference/LA_um.{pts,elem,lon} e assegna a ogni triangolo uno dei
tag di experiments/region_cv/SCENARIO.csv, a caso con la stessa probabilita'
(numpy.random.default_rng(SEED)). Scrive data/reference/SCENARIO/LA_um: .pts
e .lon copiati, .elem con gli stessi triangoli e il nuovo tag.

Uso, dalla radice del repository:
    python simulator/stage1/opencarp/random_tags.py aniso_random 0
"""

import csv
import shutil
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
ROOT   = Path(__file__).resolve().parents[3]                     # radice del repository
SRC    = ROOT / "data" / "reference" / "LA_um"                   # senza estensione
CV_DIR = ROOT / "experiments" / "region_cv"                      # SCENARIO.csv
# --------------------------------------------------------------------------


def main():
    scenario, seed = sys.argv[1], int(sys.argv[2])
    dst = ROOT / "data" / "reference" / scenario / "LA_um"
    tags = [int(r["region_tag"]) for r in
            csv.DictReader(l for l in open(CV_DIR / f"{scenario}.csv") if not l.startswith("#"))]
    tri = np.loadtxt(f"{SRC}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    tag = np.random.default_rng(seed).choice(tags, size=len(tri))

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(f"{SRC}.pts", f"{dst}.pts")
    shutil.copyfile(f"{SRC}.lon", f"{dst}.lon")
    with open(f"{dst}.elem", "w") as f:
        f.write(f"{len(tri)}\n")
        for (a, b, c), t in zip(tri, tag):
            f.write(f"Tr {a} {b} {c} {t}\n")


if __name__ == "__main__":
    main()
