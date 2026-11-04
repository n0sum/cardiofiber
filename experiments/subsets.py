#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT (passo 2), oppure
solo per i sottoinsiemi di K mappe.

Modello: il .conf TEMPLATE. Si cambiano solo experiment_id, maps,
output_csv e output_ply (none); ingressi e soglie restano quelli.
Il numero di mappe e' il numero di righe di pacing_sites.csv (la riga k e'
la mappa k). Nessuna ground truth.

experiment_id = sub_ + indici uniti da trattino (es. sub_0-2-4): finisce
non quotato nel CSV, quindi niente virgole.

Uso, dalla radice del repository:
    python experiments/subsets.py TEMPLATE OUT_DIR [K]
    python experiments/subsets.py experiments/aniso_const.conf output/subsets
I .conf vanno in OUT_DIR/confs/, le uscite in OUT_DIR/<ID>/errors.csv.
"""

import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
SITES = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
MIN_MAPS = 3                                              # kMinGradients in Mesh.cpp


def main():
    template, out_dir = Path(sys.argv[1]), sys.argv[2]
    lines = template.read_text().splitlines()
    n_maps = len(list(csv.DictReader(open(SITES))))
    conf_dir = Path(out_dir) / "confs"
    conf_dir.mkdir(parents=True, exist_ok=True)

    ks = [int(sys.argv[3])] if len(sys.argv) > 3 else range(MIN_MAPS, n_maps + 1)
    for k in ks:
        for subset in itertools.combinations(range(n_maps), k):
            exp_id = "sub_" + "-".join(map(str, subset))
            new = {"experiment_id": exp_id,
                   "maps": ",".join(map(str, subset)),
                   "output_csv": f"{out_dir}/{exp_id}/errors.csv",
                   "output_ply": "none"}
            out = []
            for line in lines:
                key, sep, _ = line.partition("=")
                if not line.startswith("#") and key.strip() in new:
                    line = f"{key}{sep} {new[key.strip()]}"
                out.append(line)
            (conf_dir / f"{exp_id}.conf").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
