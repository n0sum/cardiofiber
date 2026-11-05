#!/usr/bin/env python3
"""
Mappe degli errori per faccia in PLY (Tesi0), su aniso_const, in output/mappe/.
Per ogni caso tre PLY, uno per errore: NOME_angolo.ply, NOME_cv_l.ply,
NOME_rapporto.ply. I casi:
    faccia_6_mappe             stimatore per faccia, sei mappe
    faccia_migliore_sub_i-j-k  stimatore per faccia, terna con angolo q50 minimo
    faccia_peggiore_sub_i-j-k  stimatore per faccia, terna con angolo q50 massimo
    roney_migliore_sub_i-j-k   Roney, idem
    roney_peggiore_sub_i-j-k
Angolo q50 di ogni terna da output/roney/figure/terne.csv (analisi.py --roney):
facce stimabili nei due bracci, i punti di const.png.

Scrive in output/mappe/confs/ i .conf dello stimatore per faccia:
experiments/aniso_const.conf con cambiati experiment_id, maps, output_csv e
output_ply. Le stime di Roney sono gia' in output/roney/aniso_const/roney/:
per ognuna, dopo aver scritto i .conf, stampa una riga "ERRORS_CSV BASE", da
passare a Tesi0 --ply.

Uso, dalla radice del repository:
    python experiments/mappe.py
"""

from pathlib import Path

import pandas as pd

from subsets import write_conf

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
TEMPLATE = ROOT / "experiments" / "aniso_const.conf"
TERNE = ROOT / "output" / "roney" / "figure" / "terne.csv"
OUT = "output/mappe"


def main():
    t = pd.read_csv(TERNE)
    t = t[t["scenario"] == "aniso_const"].set_index("sub")

    faccia = {"faccia_6_mappe": "all"}                     # caso -> maps
    roney = {}                                             # caso -> errors.csv
    for arm in ("faccia", "roney"):
        q50 = t.loc[t["braccio"] == arm, "angolo_q50"]
        for nome, sub in (("migliore", q50.idxmin()), ("peggiore", q50.idxmax())):
            if arm == "faccia":
                faccia[f"faccia_{nome}_{sub}"] = sub.removeprefix("sub_").replace("-", ",")
            else:
                roney[f"roney_{nome}_{sub}"] = f"output/roney/aniso_const/roney/{sub}/errors.csv"

    for exp_id, maps in faccia.items():
        write_conf(TEMPLATE, ROOT / OUT / "confs" / f"{exp_id}.conf",
                   {"experiment_id": exp_id, "maps": maps,
                    "output_csv": f"{OUT}/{exp_id}/errors.csv", "output_ply": f"{OUT}/{exp_id}"})
    for caso, errors_csv in roney.items():
        print(f"{errors_csv} {OUT}/{caso}")


if __name__ == "__main__":
    main()
