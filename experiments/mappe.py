#!/usr/bin/env python3
"""
Mappe degli errori per faccia in PLY (Tesi0), su aniso_const, in output/mappe/:
    faccia_6_mappe.ply             stimatore per faccia, sei mappe
    faccia_migliore_sub_i-j-k.ply  stimatore per faccia, terna con angolo q50 minimo
    faccia_peggiore_sub_i-j-k.ply  stimatore per faccia, terna con angolo q50 massimo
    roney_migliore_sub_i-j-k.ply   Roney, idem
    roney_peggiore_sub_i-j-k.ply
Angolo q50 di ogni terna da output/roney/figure/terne.csv (analisi.py --roney):
facce stimabili nei due bracci, i punti di const.png.

Scrive in output/mappe/confs/ i .conf dello stimatore per faccia:
experiments/aniso_const.conf con cambiati solo experiment_id, maps, output_csv
e output_ply. Le stime di Roney sono gia' in output/roney/aniso_const/roney/:
per ognuna stampa una riga "ERRORS_CSV OUT_PLY", da passare a Tesi0 --ply.

Uso, dalla radice del repository:
    python experiments/mappe.py
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"
TERNE = ROOT / "output" / "roney" / "figure" / "terne.csv"
OUT = "output/mappe"


def main():
    t = pd.read_csv(TERNE)
    t = t[t["scenario"] == "aniso_const"].set_index("sub")
    lines = TEMPLATE.read_text().splitlines()
    (ROOT / OUT / "confs").mkdir(parents=True, exist_ok=True)

    mappe = {"faccia_6_mappe": "all"}
    for arm in ("faccia", "roney"):
        q50 = t.loc[t["braccio"] == arm, "angolo_q50"]
        for nome, sub in (("migliore", q50.idxmin()), ("peggiore", q50.idxmax())):
            if arm == "faccia":
                mappe[f"faccia_{nome}_{sub}"] = sub.removeprefix("sub_").replace("-", ",")
            else:
                print(f"output/roney/aniso_const/roney/{sub}/errors.csv {OUT}/roney_{nome}_{sub}.ply")

    for exp_id, maps in mappe.items():
        new = {"experiment_id": exp_id, "maps": maps,
               "output_csv": f"{OUT}/{exp_id}/errors.csv", "output_ply": f"{OUT}/{exp_id}.ply"}
        out = []
        for line in lines:
            key, sep, _ = line.partition("=")
            if not line.startswith("#") and key.strip() in new:
                line = f"{key}{sep} {new[key.strip()]}"
            out.append(line)
        (ROOT / OUT / "confs" / f"{exp_id}.conf").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
