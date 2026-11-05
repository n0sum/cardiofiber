#!/usr/bin/env python3
"""
Mappe degli errori in PLY per MeshLab, dagli errors.csv della griglia.
Per ogni caso stampa una riga "MESH ERRORS_CSV PREFISSO", da passare a
Tesi0 --ply, che scrive PREFISSO_{angolo,cv_l,rapporto}.ply. MESH e' il LAT
PLY dello scenario; PREFISSO e' output/mesh/METODO_SCENARIO_SUB. I casi:
    per faccia, tutte le mappe, su SCENARIO_TUTTE e su SCENARIO_TERNA
    per faccia e Roney sulla terna di SCENARIO_TERNA con l'area appaiata
    piu' grande: facce stimabili da entrambi i metodi, pesate per area

Uso, dalla radice del repository:
    python experiments/mappe.py | while read -r mesh csv prefisso; do
        ./estimator/build/Tesi0 --ply "$mesh" "$csv" "$prefisso"
    done
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
SCENARIO_TUTTE = "aniso_zone_1"
SCENARIO_TERNA = "aniso_zone_5000"


def line(method, scenario, sub):
    errors = f"output/{scenario}/{method}/{sub}/errors.csv"
    print(f"data/reference/LA_LAT_{scenario}.ply {errors} output/mesh/{method}_{scenario}_{sub}")


def main():
    (ROOT / "output" / "mesh").mkdir(parents=True, exist_ok=True)

    for scenario in (SCENARIO_TUTTE, SCENARIO_TERNA):
        subs = [d.name for d in (ROOT / "output" / scenario / "faccia").glob("sub_*")]
        line("faccia", scenario, max(subs, key=lambda s: s.count("-")))   # tutte le mappe

    paired = {}                                           # terna -> area appaiata
    for d in (ROOT / "output" / SCENARIO_TERNA / "roney").glob("sub_*"):
        roney = pd.read_csv(d / "errors.csv")
        faccia = pd.read_csv(ROOT / "output" / SCENARIO_TERNA / "faccia" / d.name / "errors.csv")
        both = (roney["estimable"] == 1) & (faccia["estimable"] == 1)
        paired[d.name] = roney["area"][both].sum() / roney["area"].sum()
    sub = max(paired, key=paired.get)
    line("faccia", SCENARIO_TERNA, sub)
    line("roney", SCENARIO_TERNA, sub)


if __name__ == "__main__":
    main()
