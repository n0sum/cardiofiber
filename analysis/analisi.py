#!/usr/bin/env python3
# analisi.py si usa solo su esperimenti con maps = all
"""
Tabelle di errore di un esperimento, da output/<esperimento>/errors.csv:
  a) riepilogo
  b) per tag di regione
  c) errore angolare per fasce di condizionamento

Tutte le statistiche sono PESATE PER AREA, non per conteggio di facce.

Uso, dalla radice del repository:
    python analysis/analisi.py aniso_const
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
K = len(pd.read_csv(ROOT / "simulator" / "stage1" / "pacing_sites.csv"))   # mappe dell'esperimento


def wquantile(values, weights, q):
    """Quantile pesato. Ignora i NaN."""
    ok = np.isfinite(values) & np.isfinite(weights)
    v, w = np.asarray(values)[ok], np.asarray(weights)[ok]
    if v.size == 0:
        return np.nan
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    return np.interp(q * cw[-1], cw - 0.5 * w, v)


def rule(title):
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


def summary(g):
    """Colonne della tabella a). err_angle e' vuoto sulle facce non stimabili."""
    w = g["area"]
    est = g["estimable"] == 1
    return {
        "q50": wquantile(g["err_angle"], w, 0.50),
        "q90": wquantile(g["err_angle"], w, 0.90),
        "<20° %": 100 * w[g["err_angle"] < 20].sum() / w[est].sum(),   # sull'area stimabile
        "non stim. %": 100 * w[~est].sum() / w.sum(),
        "coppie escl. %": 100 * (w * (K - g["k_used"])).sum() / (K * w.sum()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment", help="nome della cartella in output/")
    exp = ap.parse_args().experiment
    df = pd.read_csv(ROOT / "output" / exp / "errors.csv")
    fmt = lambda x: f"{x:8.3f}"

    rule(f"a) RIEPILOGO  -  {exp}")
    print(pd.DataFrame([summary(df)]).to_string(index=False, float_format=fmt))

    rule("b) PER REGIONE")
    rows = []
    for tag, g in df.groupby("region_tag"):
        w = g["area"]
        rows.append({
            "tag": tag,
            "CV_L est": wquantile(g["CV_L_est"], w, 0.50),
            "CV_L vera": wquantile(g["CV_L_true"], w, 0.50),
            "CV_T est": wquantile(g["CV_T_est"], w, 0.50),
            "CV_T vera": wquantile(g["CV_T_true"], w, 0.50),
            "q50": wquantile(g["err_angle"], w, 0.50),
            "q90": wquantile(g["err_angle"], w, 0.90),
            "non stim. %": 100 * w[g["estimable"] == 0].sum() / w.sum(),
        })
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))

    rule("c) ERRORE ANGOLARE PER FASCE DI COND  (facce stimabili, fasce con almeno 20 facce)")
    est = df[df["estimable"] == 1]
    c = est["cond"].replace([np.inf, -np.inf], np.nan).dropna()
    bins = np.logspace(np.log10(max(c.min(), 1.0)), np.log10(c.max()), 13)
    idx = np.digitize(est["cond"], bins)
    print("  fascia di cond          facce      q50      q90")
    for b in range(1, len(bins)):
        g = est[idx == b]
        if len(g) < 20:
            continue
        w = g["area"]
        print(f"  [{bins[b-1]:8.1f},{bins[b]:9.1f})  {len(g):7d} "
              f"{wquantile(g['err_angle'], w, 0.50):8.3f} {wquantile(g['err_angle'], w, 0.90):8.3f}")


if __name__ == "__main__":
    main()
