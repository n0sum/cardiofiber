#!/usr/bin/env python3
"""
Genera un .conf per ogni scenario aniso_zone_N.

Modello: experiments/aniso_const.conf (tutte le mappe, maschera a 30 gradi).
Si cambiano solo experiment_id, lat_ply, ground_truth_ply e output_csv.

Uso, dalla radice del repository:
    python experiments/zone.py N [N ...]
I .conf vanno in output/zone/confs/, le uscite in output/zone/aniso_zone_N/.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"
CONF_DIR = ROOT / "output" / "zone" / "confs"


def main():
    lines = TEMPLATE.read_text().splitlines()
    CONF_DIR.mkdir(parents=True, exist_ok=True)

    for n in sys.argv[1:]:
        exp_id = f"aniso_zone_{n}"
        new = {"experiment_id": exp_id,
               "lat_ply": f"data/reference/LA_LAT_{exp_id}.ply",
               "ground_truth_ply": f"data/reference/LA_GT_{exp_id}.ply",
               "output_csv": f"output/zone/{exp_id}/errors.csv"}
        out = []
        for line in lines:
            key, sep, _ = line.partition("=")
            if not line.startswith("#") and key.strip() in new:
                line = f"{key}{sep} {new[key.strip()]}"
            out.append(line)
        (CONF_DIR / f"{exp_id}.conf").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
