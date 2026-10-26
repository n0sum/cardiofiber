#!/usr/bin/env python3
"""
Genera un .conf per ogni sottoinsieme di k >= 3 mappe LAT (passo 2).

Modello: experiments/aniso_const.conf. Si cambiano solo experiment_id,
maps e output_csv; ingressi e soglie restano quelli.
Le mappe e il loro ordine vengono dai commenti "lat_k NOME" dell'header
di lat_ply. Nessuna ground truth.

experiment_id = sub_ + indici uniti da trattino (es. sub_0-2-4): finisce
non quotato nel CSV, quindi niente virgole.

Uso, dalla radice del repository:
    python experiments/subsets.py
I .conf vanno in output/subsets/confs/, le uscite in output/subsets/<ID>/.
"""

import itertools
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"
CONF_DIR = ROOT / "output" / "subsets" / "confs"
MIN_MAPS = 3                                              # kMinGradients in Mesh.cpp


def conf_value(lines, key):
    for line in lines:
        k, _, v = line.partition("=")
        if not line.startswith("#") and k.strip() == key:
            return v.strip()


def map_names(ply):
    """Nomi delle mappe dai commenti 'lat_k NOME' dell'header, in ordine di k."""
    names = {}
    with open(ply, "rb") as f:
        for raw in f:
            line = raw.decode("ascii").strip()
            if line == "end_header":
                break
            parts = line.split()
            if len(parts) == 3 and parts[0] == "comment" and parts[1].startswith("lat_"):
                names[int(parts[1][4:])] = parts[2]
    return [names[k] for k in range(len(names))]


def main():
    lines = TEMPLATE.read_text().splitlines()
    names = map_names(ROOT / conf_value(lines, "lat_ply"))
    CONF_DIR.mkdir(parents=True, exist_ok=True)

    n = 0
    for k in range(MIN_MAPS, len(names) + 1):
        for subset in itertools.combinations(range(len(names)), k):
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
            n += 1

    print(f"mappe: {', '.join(f'{i}={s}' for i, s in enumerate(names))}")
    print(f"{n} .conf scritti in {CONF_DIR.relative_to(ROOT)}/")
    print(f"\n./estimator/build/Tesi0 {CONF_DIR.relative_to(ROOT)}/*.conf")


if __name__ == "__main__":
    main()
