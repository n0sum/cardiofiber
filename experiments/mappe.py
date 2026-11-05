#!/usr/bin/env python3
"""
Mappe degli errori per faccia in PLY (Tesi0), su aniso_zone_1, in output/mappe/.
Per ogni caso tre PLY, uno per errore: CASO_angolo.ply, CASO_cv_l.ply,
CASO_rapporto.ply. I casi:
    faccia_tutte_le_mappe      stimatore per faccia, tutte le mappe che convergono
    faccia_migliore_sub_i-j-k  stimatore per faccia, terna con angolo q50 minimo
    faccia_peggiore_sub_i-j-k  stimatore per faccia, terna con angolo q50 massimo
    roney_migliore_sub_i-j-k   Roney, idem
    roney_peggiore_sub_i-j-k
Angolo q50 di ogni terna da output/roney/figure/terne.csv (analisi.py --roney):
facce stimabili nei due bracci, i punti di zona_1.png. Le terne sono solo
quelle delle mappe che convergono.

Scrive in output/mappe/confs/ i .conf dello stimatore per faccia. Le stime di
Roney sono gia' in output/roney/aniso_zone_1/roney/: per ognuna, dopo aver
scritto i .conf, stampa una riga "ERRORS_CSV BASE", da passare a Tesi0 --ply.

Uso, dalla radice del repository:
    python experiments/mappe.py
"""

from pathlib import Path

import pandas as pd

from subsets import write_conf

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
SCENARIO = "aniso_zone_1"
TERNE = ROOT / "output" / "roney" / "figure" / "terne.csv"
OUT = "output/mappe"


def main():
    t = pd.read_csv(TERNE)
    t = t[t["scenario"] == SCENARIO].set_index("sub")

    faccia = {"faccia_tutte_le_mappe": "all"}              # caso -> maps
    roney = {}                                             # caso -> errors.csv
    for arm in ("faccia", "roney"):
        q50 = t.loc[t["braccio"] == arm, "angolo_q50"]
        for nome, sub in (("migliore", q50.idxmin()), ("peggiore", q50.idxmax())):
            if arm == "faccia":
                faccia[f"faccia_{nome}_{sub}"] = sub.removeprefix("sub_").replace("-", ",")
            else:
                roney[f"roney_{nome}_{sub}"] = f"output/roney/{SCENARIO}/roney/{sub}/errors.csv"

    for caso, maps in faccia.items():
        write_conf(ROOT / OUT / "confs" / f"{caso}.conf", SCENARIO, caso, maps,
                   f"{OUT}/{caso}/errors.csv", f"{OUT}/{caso}")
    for caso, errors_csv in roney.items():
        print(f"{errors_csv} {OUT}/{caso}")


if __name__ == "__main__":
    main()
