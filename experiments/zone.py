#!/usr/bin/env python3
"""
Genera un .conf per ogni scenario aniso_zone_N: experiments/aniso_const.conf
con cambiati experiment_id, lat_ply, ground_truth_ply, output_csv e
output_ply (none).

Uso, dalla radice del repository:
    python experiments/zone.py N [N ...]
I .conf vanno in output/zone/confs/, le uscite in output/zone/aniso_zone_N/.
"""

import sys
from pathlib import Path

from subsets import write_conf

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"


def main():
    for n in sys.argv[1:]:
        exp_id = f"aniso_zone_{n}"
        write_conf(TEMPLATE, ROOT / "output" / "zone" / "confs" / f"{exp_id}.conf",
                   {"experiment_id": exp_id,
                    "lat_ply": f"data/reference/LA_LAT_{exp_id}.ply",
                    "ground_truth_ply": f"data/reference/LA_GT_{exp_id}.ply",
                    "output_csv": f"output/zone/{exp_id}/errors.csv",
                    "output_ply": "none"})


if __name__ == "__main__":
    main()
