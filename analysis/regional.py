#!/usr/bin/env python3
"""
Valutazione dello scenario con CV per regione (aniso_regional) contro lo
scenario a CV uniformi (aniso_const). Solo lettura: errors.csv,
mask_pairs.csv e la mesh (tag per faccia da LA_endo_with_fiber.elem).

Una faccia e' DI CONFINE se condivide un vertice con una faccia di tag
diverso, altrimenti e' INTERNA. Le viste:
  1. per regione: CV_L, CV_T e rapporto stimati (mediane) contro imposti,
     err_angle q50 e q90, area non stimabile
  2. le stesse statistiche separate per facce interne e di confine, e
     quanta area di ogni regione e' di confine
  3. modello nullo: ovunque la media globale (pesata per area) delle CV
     vere. Errore relativo mediano |stima - vero| / vero del nullo contro
     quello dello stimatore, sulle facce interne stimabili
  4. maschera ai confini: frazione di coppie escluse e mediana di |r|
     delle escluse, interne contro di confine, nei due scenari sulle
     stesse facce
  5. confronto appaiato con il riferimento sulle facce stimabili in
     entrambi: dq50 e dq90 di err_angle (scenario - riferimento)

Tutte le statistiche sono PESATE PER AREA, con wquantile di analisi.py.

Uso, dalla radice del repository:
    python analysis/regional.py aniso_regional aniso_const
"""

import argparse

import numpy as np
import pandas as pd

from analisi import REGION_NAMES, ROOT, rule, wquantile

ELEM = ROOT / "data" / "augmenta" / "LA_endo_with_fiber.elem"
ZONES = ("tutte", "interne", "confine")


def border_faces():
    """True per le facce che condividono un vertice con una faccia di tag diverso."""
    el = np.loadtxt(ELEM, skiprows=1, usecols=(1, 2, 3, 4), dtype=np.int64)
    tri, tag = el[:, :3], el[:, 3]
    n = tri.max() + 1
    vmin, vmax = np.full(n, tag.max()), np.full(n, tag.min())
    for j in range(3):
        np.minimum.at(vmin, tri[:, j], tag)
        np.maximum.at(vmax, tri[:, j], tag)
    return (vmin[tri] != vmax[tri]).any(1), tag


def load(exp, border, tag):
    df = pd.read_csv(ROOT / "output" / exp / "errors.csv")
    if not (df["region_tag"].to_numpy() == tag).all():
        raise SystemExit(f"{exp}: region_tag di errors.csv diverso dai tag di {ELEM}")
    df["border"] = border
    mp = pd.read_csv(ROOT / "output" / exp / "mask_pairs.csv")
    mp = mp.merge(df[["face_id", "area", "region_tag", "border"]], on="face_id")
    return df, mp


def zone(df, z):
    return df if z == "tutte" else df[df["border"] == (z == "confine")]


def groups(df):
    """(regione, zona, sottoinsieme) per ogni regione e per la mesh intera."""
    for t in [*REGION_NAMES, 0]:
        sub = df if t == 0 else df[df["region_tag"] == t]
        for z in ZONES:
            yield ("tutte" if t == 0 else REGION_NAMES[t]), z, zone(sub, z)


def wmed(v, w):
    return wquantile(v.to_numpy(), w.to_numpy(), 0.5)


def show(rows):
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))


def main():
    ap = argparse.ArgumentParser(description="CV per regione contro CV uniformi.")
    ap.add_argument("scenario")
    ap.add_argument("reference")
    args = ap.parse_args()

    border, tag = border_faces()
    df, mp = load(args.scenario, border, tag)
    ref, mpr = load(args.reference, border, tag)

    rule(f"1-2. STIMATO CONTRO IMPOSTO  -  {args.scenario}")
    rows = []
    for name, z, s in groups(df):
        e = s[s["estimable"] == 1]
        rows.append({"regione": name, "zona": z,
                     "area_%": 100 * s["area"].sum() / df["area"].sum(),
                     "L_imp": wmed(s["CV_L_true"], s["area"]), "L_est": wmed(e["CV_L_est"], e["area"]),
                     "T_imp": wmed(s["CV_T_true"], s["area"]), "T_est": wmed(e["CV_T_est"], e["area"]),
                     "rap_imp": wmed(s["CV_L_true"] / s["CV_T_true"], s["area"]),
                     "rap_est": wmed(e["CV_L_est"] / e["CV_T_est"], e["area"]),
                     "ang_q50": wmed(e["err_angle"], e["area"]),
                     "ang_q90": wquantile(e["err_angle"].to_numpy(), e["area"].to_numpy(), 0.9),
                     "nonstim_%": 100 * s.loc[s["estimable"] == 0, "area"].sum() / s["area"].sum()})
    show(rows)
    print("\narea_%: quota dell'area totale; L, T in mm/ms; rap = CV_L / CV_T per faccia.")
    print("Quota di confine per regione = area_% confine / area_% tutte.")

    rule("3. MODELLO NULLO (media globale delle CV vere) CONTRO STIMATORE  -  facce interne stimabili")
    w = df["area"] / df["area"].sum()
    null_L, null_T = (w * df["CV_L_true"]).sum(), (w * df["CV_T_true"]).sum()
    print(f"nullo: CV_L = {null_L:.4f}   CV_T = {null_T:.4f} mm/ms")
    rows = []
    for name, z, s in groups(df):
        if z != "interne":
            continue
        e = s[s["estimable"] == 1]
        rel = lambda v, t: (v - t).abs() / t
        rows.append({"regione": name,
                     "L_nullo_%": 100 * wmed(rel(null_L, e["CV_L_true"]), e["area"]),
                     "L_stima_%": 100 * wmed(rel(e["CV_L_est"], e["CV_L_true"]), e["area"]),
                     "T_nullo_%": 100 * wmed(rel(null_T, e["CV_T_true"]), e["area"]),
                     "T_stima_%": 100 * wmed(rel(e["CV_T_est"], e["CV_T_true"]), e["area"])})
    show(rows)

    rule(f"4. MASCHERA: coppie escluse e mediana |r| delle escluse  -  {args.scenario} contro {args.reference}")
    rows = []
    for (name, z, s), (_, _, sr) in zip(groups(mp), groups(mpr)):
        row = {"regione": name, "zona": z}
        for lab, x in (("scen", s), ("ref", sr)):
            ex = x[x["masked"] == 1]
            row[f"escl_%_{lab}"] = 100 * (x["area"] * x["masked"]).sum() / x["area"].sum()
            row[f"|r|_{lab}"] = wmed(ex["r"].abs(), ex["area"])
        rows.append(row)
    show(rows)
    print(f"\nscen = {args.scenario}, ref = {args.reference}, stesse facce e stesse mappe.")
    print("escl_%: area delle coppie escluse / area di tutte le coppie.")

    rule(f"5. CONFRONTO APPAIATO err_angle: {args.scenario} - {args.reference}, facce stimabili in entrambi")
    both = (df["estimable"] == 1) & (ref["estimable"] == 1)
    rows = []
    for (name, z, s), (_, _, sr) in zip(groups(df[both]), groups(ref[both])):
        a, ar = s["area"].to_numpy(), sr["area"].to_numpy()
        q = lambda x, aa, p: wquantile(x["err_angle"].to_numpy(), aa, p)
        rows.append({"regione": name, "zona": z, "area_%": 100 * a.sum() / df["area"].sum(),
                     "q50": q(s, a, 0.5), "q50_ref": q(sr, ar, 0.5), "dq50": q(s, a, 0.5) - q(sr, ar, 0.5),
                     "q90": q(s, a, 0.9), "q90_ref": q(sr, ar, 0.9), "dq90": q(s, a, 0.9) - q(sr, ar, 0.9)})
    show(rows)
    print("\nangoli in gradi; dq > 0: lo scenario sbaglia di piu' del riferimento sulle stesse facce.")


if __name__ == "__main__":
    main()
