#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT di uno scenario,
oppure solo per i sottoinsiemi di K mappe. Il numero di mappe e' il numero
di righe di pacing_sites.csv (la riga k e' la mappa k).

experiment_id = sub_ + indici uniti da trattino (es. sub_0-2-4): finisce
non quotato nel CSV, quindi niente virgole.

Uso, dalla radice del repository:
    python experiments/subsets.py SCENARIO OUT_DIR [K]
    python experiments/subsets.py aniso_zone_1 output/subsets
I .conf vanno in OUT_DIR/confs/, le uscite in OUT_DIR/<ID>/errors.csv.
"""

import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
MODEL = ROOT / "experiments" / "modello.conf"
SITES = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
MIN_MAPS = 3                                              # kMinGradients in Mesh.cpp


def write_conf(path, scenario, experiment_id, maps, output_csv, output_ply="none"):
    """Scrive in path un .conf: le soglie di modello.conf piu' scenario, mappe e uscite.
    La usano anche zone.py e mappe.py."""
    keys = {"experiment_id": experiment_id,
            "lat_ply": f"data/reference/LA_LAT_{scenario}.ply",
            "ground_truth_ply": f"data/reference/LA_GT_{scenario}.ply",
            "maps": maps,
            "output_csv": output_csv,
            "output_ply": output_ply}
    lines = [l for l in MODEL.read_text().splitlines() if l.strip() and not l.startswith("#")]
    lines += [f"{k} = {v}" for k, v in keys.items()]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main():
    scenario, out_dir = sys.argv[1], sys.argv[2]
    maps = range(len(list(csv.DictReader(open(SITES)))))
    ks = [int(sys.argv[3])] if len(sys.argv) > 3 else range(MIN_MAPS, len(maps) + 1)
    for k in ks:
        for subset in itertools.combinations(maps, k):
            exp_id = "sub_" + "-".join(map(str, subset))
            write_conf(Path(out_dir) / "confs" / f"{exp_id}.conf", scenario, exp_id,
                       ",".join(map(str, subset)), f"{out_dir}/{exp_id}/errors.csv")


if __name__ == "__main__":
    main()
