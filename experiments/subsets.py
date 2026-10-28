#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT (passo 2).

Modello: experiments/aniso_const.conf. Si cambiano solo experiment_id,
maps e output_csv; ingressi e soglie restano quelli.
Il numero di mappe e' il numero di righe di pacing_sites.csv (la riga k e'
la mappa k). Nessuna ground truth.

experiment_id = sub_ + indici uniti da trattino (es. sub_0-2-4): finisce
non quotato nel CSV, quindi niente virgole.

Uso, dalla radice del repository:
    python experiments/subsets.py
I .conf vanno in output/subsets/confs/, le uscite in output/subsets/<ID>/.
"""

import csv
import itertools
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"
SITES = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
CONF_DIR = ROOT / "output" / "subsets" / "confs"
MIN_MAPS = 3                                              # kMinGradients in Mesh.cpp


def main():
    lines = TEMPLATE.read_text().splitlines()
    n_maps = len(list(csv.DictReader(open(SITES))))
    CONF_DIR.mkdir(parents=True, exist_ok=True)

    for k in range(MIN_MAPS, n_maps + 1):
        for subset in itertools.combinations(range(n_maps), k):
            exp_id = "sub_" + "-".join(map(str, subset))
            new = {"experiment_id": exp_id,
                   "maps": ",".join(map(str, subset)),
                   "output_csv": f"output/subsets/{exp_id}/errors.csv"}
            out = []
            for line in lines:
                key, sep, _ = line.partition("=")
                if not line.startswith("#") and key.strip() in new:
                    line = f"{key}{sep} {new[key.strip()]}"
                out.append(line)
            (CONF_DIR / f"{exp_id}.conf").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
