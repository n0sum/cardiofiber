#!/usr/bin/env python3
"""
Genera un .conf per ogni scenario aniso_zone_N, con tutte le mappe che
convergono (maps = all).

Uso, dalla radice del repository:
    python experiments/zone.py N [N ...]
I .conf vanno in output/zone/confs/, le uscite in output/zone/aniso_zone_N/.
"""

import sys
from pathlib import Path

from subsets import write_conf

ROOT = Path(__file__).resolve().parents[1]                # radice del repository


def main():
    for n in sys.argv[1:]:
        scenario = f"aniso_zone_{n}"
        write_conf(ROOT / "output" / "zone" / "confs" / f"{scenario}.conf", scenario, scenario,
                   "all", f"output/zone/{scenario}/errors.csv")


if __name__ == "__main__":
    main()
