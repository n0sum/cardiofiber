#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT di uno scenario
(con sei mappe 20 + 15 + 6 + 1 = 42). Il numero di mappe e' il numero di
righe di pacing_sites.csv (la riga k e' la mappa k).

experiment_id = sub_ + indici uniti da trattino (es. sub_0-2-4): finisce
non quotato nel CSV, quindi niente virgole.

Uso, dalla radice del repository:
    python experiments/subsets.py SCENARIO
I .conf vanno in output/SCENARIO/confs/, le uscite di Tesi0 in
output/SCENARIO/faccia/<ID>/errors.csv.
"""

import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
MODEL = ROOT / "experiments" / "modello.conf"
SITES = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
MIN_MAPS = 3                                              # kMinGradients in Mesh.cpp


def write_conf(path, scenario, experiment_id, maps, output_csv):
    """Scrive in path un .conf: le soglie di modello.conf piu' scenario, mappe e uscita."""
    keys = {"experiment_id": experiment_id,
            "lat_ply": f"data/reference/LA_LAT_{scenario}.ply",
            "ground_truth_ply": f"data/reference/LA_GT_{scenario}.ply",
            "maps": maps,
            "output_csv": output_csv}
    lines = [l for l in MODEL.read_text().splitlines() if l.strip() and not l.startswith("#")]
    lines += [f"{k} = {v}" for k, v in keys.items()]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main():
    scenario = sys.argv[1]
    out_dir = f"output/{scenario}"
    maps = range(len(list(csv.DictReader(open(SITES)))))
    for k in range(MIN_MAPS, len(maps) + 1):
        for subset in itertools.combinations(maps, k):
            exp_id = "sub_" + "-".join(map(str, subset))
            write_conf(ROOT / out_dir / "confs" / f"{exp_id}.conf", scenario, exp_id,
                       ",".join(map(str, subset)), f"{out_dir}/faccia/{exp_id}/errors.csv")


if __name__ == "__main__":
    main()
