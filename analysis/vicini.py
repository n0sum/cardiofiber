#!/usr/bin/env python3
"""
Errore della stima contro la verita' della faccia e contro la verita' del vicinato.

Vicinato di una faccia: le facce che toccano almeno uno dei suoi tre vertici, lei
compresa. Verita' del vicinato: media pesata per area di CV_L_true e di CV_T_true sul
vicinato; r del vicinato = rapporto delle due medie.

Sulle facce stimabili, pesato per area, q50 e q90 del valore assoluto dell'errore
relativo su CV_L e su r = CV_L / CV_T, in tre versioni:
  faccia   stima contro verita' della faccia (come numeri.csv)
  vicinato stima contro verita' del vicinato
  nullo    costante = mediana pesata della verita' su tutta la mesh, contro verita'
           della faccia

Uso, dalla radice del repository:
    python analysis/vicini.py aniso_random
        -> output/aniso_random/figure/vicini.csv
"""

import argparse

import numpy as np
import pandas as pd

from analisi import ROOT, load, wquantile

ELEM = ROOT / "data" / "augmenta" / "LA_endo_with_fiber.elem"   # stesso ordine dei face_id


def neighbourhood_mean(tri, w, values):
    """Media pesata con w di ogni colonna di values sul vicinato di ogni faccia."""
    fv = pd.DataFrame({"f": np.repeat(np.arange(len(tri)), 3), "v": tri.ravel()})
    pairs = fv.merge(fv, on="v").drop_duplicates(["f_x", "f_y"])
    src = pairs["f_y"].values
    num = pd.DataFrame(values[src] * w[src, None]).groupby(pairs["f_x"].values).sum()
    den = pd.Series(w[src]).groupby(pairs["f_x"].values).sum()
    return num.values / den.values[:, None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment", help="nome della cartella in output/")
    name = ap.parse_args().experiment

    df = load(ROOT / "output" / name / "errors.csv").sort_index()
    tri = np.loadtxt(ELEM, skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    if len(tri) != len(df):
        raise SystemExit(f"{ELEM}: {len(tri)} triangoli, errors.csv ne ha {len(df)}")

    w = df["area"].values
    l_true, t_true = df["CV_L_true"].values, df["CV_T_true"].values
    l_vic, t_vic = neighbourhood_mean(tri, w, np.c_[l_true, t_true]).T
    l_est, r_est = df["CV_L_est"].values, (df["CV_L_est"] / df["CV_T_est"]).values
    r_true = l_true / t_true

    def rel(est, true):
        return 100 * np.abs(est - true) / true

    versions = {
        "faccia": (rel(l_est, l_true), rel(r_est, r_true)),
        "vicinato": (rel(l_est, l_vic), rel(r_est, l_vic / t_vic)),
        "nullo": (rel(wquantile(l_true, w, 0.5), l_true), rel(wquantile(r_true, w, 0.5), r_true)),
    }
    est = df["estimable"].values == 1
    rows = []
    for v, (cv_l, rapporto) in versions.items():
        rows.append({"versione": v,
                     "cv_l_q50": wquantile(cv_l[est], w[est], 0.50),
                     "cv_l_q90": wquantile(cv_l[est], w[est], 0.90),
                     "rapporto_q50": wquantile(rapporto[est], w[est], 0.50),
                     "rapporto_q90": wquantile(rapporto[est], w[est], 0.90)})
    out_dir = ROOT / "output" / name / "figure"
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_dir / "vicini.csv", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
