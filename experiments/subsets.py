#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT, oppure solo per i
sottoinsiemi di K mappe.

Ogni .conf e' il modello TEMPLATE con cambiati experiment_id, maps,
output_csv e output_ply (none). Il numero di mappe e' il numero di righe di
pacing_sites.csv (la riga k e' la mappa k).

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


def write_conf(template, path, new):
    """Scrive in path il .conf template senza commenti, con i valori di new al posto dei suoi.
    La usano anche zone.py e mappe.py."""
    out = []
    for line in Path(template).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            key, value = (s.strip() for s in line.split("=", 1))
            out.append(f"{key} = {new.get(key, value)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n")


def main():
    template, out_dir = sys.argv[1], sys.argv[2]
    n_maps = len(list(csv.DictReader(open(SITES))))
    ks = [int(sys.argv[3])] if len(sys.argv) > 3 else range(MIN_MAPS, n_maps + 1)
    for k in ks:
        for subset in itertools.combinations(range(n_maps), k):
            exp_id = "sub_" + "-".join(map(str, subset))
            write_conf(template, Path(out_dir) / "confs" / f"{exp_id}.conf",
                       {"experiment_id": exp_id,
                        "maps": ",".join(map(str, subset)),
                        "output_csv": f"{out_dir}/{exp_id}/errors.csv",
                        "output_ply": "none"})


if __name__ == "__main__":
    main()
