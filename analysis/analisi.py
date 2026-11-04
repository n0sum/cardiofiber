#!/usr/bin/env python3
"""
Errori di stima dagli errors.csv. Tutte le statistiche sono PESATE PER AREA.

Per ogni faccia stimabile:
  angolo    err_angle, gradi
  cv_l      (CV_L_est - CV_L_true) / CV_L_true, con segno
  rapporto  (r_est - r_true) / r_true, con segno, r = CV_L / CV_T

Uso, dalla radice del repository:
    python analysis/analisi.py aniso_const
        -> output/aniso_const/figure/{errori.png, numeri.csv}
    python analysis/analisi.py --sottoinsiemi
        -> output/subsets/figure/{area_non_stimabile,angolo,cv_l,rapporto}.png, numeri.csv
           da output/subsets/sub_*/errors.csv (generati da experiments/subsets.py)
    python analysis/analisi.py --zone
        -> output/zone/figure/{zone.png, numeri.csv}
           da output/zone/aniso_zone_*/errors.csv (generati da experiments/zone.py)
    python analysis/analisi.py --roney
        -> output/roney/figure/{const.png, zone.png, numeri.csv, terne.csv}
           da output/roney/SCENARIO/{roney,faccia}/sub_*/errors.csv (estimator/roney.py)

Sottoinsiemi: confronto APPAIATO, gli errori di ogni sottoinsieme si calcolano
sulle facce stimabili sia in esso sia in sub_0-1-2-3-4-5.

Zone: per ogni aniso_zone_N il diametro medio di zona d = sqrt(4 A / (pi N)),
A = area totale in mm^2; riferimento orizzontale: aniso_const.

Roney: per scenario e terna, confronto APPAIATO fra i due bracci (facce
stimabili in entrambi); l'area non stimabile e' di ogni braccio da solo.
numeri.csv ha una riga per scenario e braccio, terne.csv una per terna.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
K = len(pd.read_csv(ROOT / "simulator" / "stage1" / "pacing_sites.csv"))   # mappe totali
SUBSETS = ROOT / "output" / "subsets"
ZONE = ROOT / "output" / "zone"
RONEY = ROOT / "output" / "roney"
ARMS = {"roney": "#eb6834", "faccia": "#2a78d6"}
ERRORS = {"angolo": "errore angolare [gradi]",
          "cv_l": "errore relativo su CV_L [%]",
          "rapporto": "errore relativo su CV_L / CV_T [%]"}
XLIM = {"angolo": (1e-2, 90), "cv_l": (-20, 20), "rapporto": (-40, 40)}


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


def load(path):
    """errors.csv con le colonne degli errori; relativi in percentuale, NaN sulle non stimabili."""
    df = pd.read_csv(path).set_index("face_id")
    df["angolo"] = df["err_angle"]
    df["cv_l"] = 100 * (df["CV_L_est"] - df["CV_L_true"]) / df["CV_L_true"]
    r_est, r_true = df["CV_L_est"] / df["CV_T_est"], df["CV_L_true"] / df["CV_T_true"]
    df["rapporto"] = 100 * (r_est - r_true) / r_true
    return df


def nonstim(df):
    """% di area non stimabile."""
    w = df["area"]
    return 100 * w[df["estimable"] == 0].sum() / w.sum()


def quantiles(df, prefix=""):
    """q50 e q90 di ogni errore sulle facce stimabili (valore assoluto per i relativi), piu' la
    mediana con segno dei relativi."""
    est = df[df["estimable"] == 1]
    w = est["area"]
    out = {}
    for e in ERRORS:
        a = est[e].abs()
        out[f"{e}_q50"] = wquantile(a, w, 0.50)
        out[f"{e}_q90"] = wquantile(a, w, 0.90)
        if e != "angolo":
            out[f"{e}_mediana"] = wquantile(est[e], w, 0.50)
    return out


def experiment(name):
    df = load(ROOT / "output" / name / "errors.csv")
    out_dir = ROOT / "output" / name / "figure"
    out_dir.mkdir(parents=True, exist_ok=True)
    ns = nonstim(df)
    pd.DataFrame([{"area_non_stimabile_%": ns, **quantiles(df)}]) \
        .to_csv(out_dir / "numeri.csv", index=False, float_format="%.6g")

    est = df[df["estimable"] == 1]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (e, label) in zip(axes, ERRORS.items()):
        order = np.argsort(est[e].values)
        x = est[e].values[order]
        y = 100 * np.cumsum(est["area"].values[order]) / est["area"].sum()
        ax.plot(x, y, color="#2a78d6", lw=1.8)
        for q in (50, 90):
            ax.axhline(q, color="#888", lw=1, ls="--")
        if e == "angolo":
            ax.set_xscale("log")
        else:
            ax.axvline(0, color="#444", lw=1)
        ax.set_xlim(*XLIM[e])
        ax.set_ylim(0, 100)
        ax.set_xlabel(label)
        ax.set_ylabel("area stimabile con errore <= x [%]")
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"{name}: area non stimabile {ns:.2f}%")
    plt.tight_layout()
    plt.savefig(out_dir / "errori.png", dpi=140)
    plt.close()


def subsets():
    runs = {}
    for d in SUBSETS.glob("sub_*"):
        runs[tuple(int(i) for i in d.name[4:].split("-"))] = load(d / "errors.csv")
    ref = runs[tuple(range(K))]

    rows = []
    for idx in sorted(runs, key=lambda i: (len(i), i)):
        df = runs[idx]
        paired = df[(df["estimable"] == 1) & (ref["estimable"] == 1)]
        rows.append({"sottoinsieme": "+".join(map(str, idx)), "k": len(idx),
                     "area_non_stimabile_%": nonstim(df), **quantiles(paired)})
    s = pd.DataFrame(rows)
    out_dir = SUBSETS / "figure"
    out_dir.mkdir(parents=True, exist_ok=True)
    s.to_csv(out_dir / "numeri.csv", index=False, float_format="%.6g")

    # sfalsamento orizzontale dentro ogni k
    i, n = s.groupby("k").cumcount(), s.groupby("k")["k"].transform("size")
    s["x"] = s["k"] + 0.6 * (i - (n - 1) / 2) / n

    def figure(file, cols, ylabel, title):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for col, marker, color, label in cols:
            ax.scatter(s["x"], s[col], marker=marker, s=32, color=color, label=label, zorder=3)
        ax.set_xticks(sorted(s["k"].unique()))
        ax.set_xlabel("mappe usate (k)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        if len(cols) > 1:
            ax.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(out_dir / file, dpi=140)
        plt.close()

    figure("area_non_stimabile.png", [("area_non_stimabile_%", "o", "#2a78d6", None)],
           "area non stimabile [%]", "area non stimabile, un punto per sottoinsieme")
    for e, label in ERRORS.items():
        figure(f"{e}.png", [(f"{e}_q50", "o", "#2a78d6", "q50"), (f"{e}_q90", "^", "#eb6834", "q90")],
               label + (", valore assoluto" if e != "angolo" else ""),
               "facce stimabili nel sottoinsieme e con tutte le mappe")


def zone():
    rows = []
    for d in ZONE.glob("aniso_zone_*"):
        df = load(d / "errors.csv")
        n = int(d.name.removeprefix("aniso_zone_"))
        rows.append({"N": n, "d": np.sqrt(4 * df["area"].sum() / (np.pi * n)),
                     "area_non_stimabile_%": nonstim(df), **quantiles(df)})
    z = pd.DataFrame(rows).sort_values("d")
    out_dir = ZONE / "figure"
    out_dir.mkdir(parents=True, exist_ok=True)
    z.to_csv(out_dir / "numeri.csv", index=False, float_format="%.6g")

    c = load(ROOT / "output" / "aniso_const" / "errors.csv")
    ref = {"area_non_stimabile_%": nonstim(c), **quantiles(c)}

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    panels = [("area non stimabile [%]", [("area_non_stimabile_%", "o", "#2a78d6", None)])]
    for e, label in ERRORS.items():
        panels.append((label + (", valore assoluto" if e != "angolo" else ""),
                       [(f"{e}_q50", "o", "#2a78d6", "q50"), (f"{e}_q90", "^", "#eb6834", "q90")]))
    for ax, (ylabel, cols) in zip(axes, panels):
        for col, marker, color, label in cols:
            ax.plot(z["d"], z[col], marker=marker, color=color, lw=1, label=label)
            ax.axhline(ref[col], color=color, lw=1, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("diametro medio di zona d [mm]")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3, which="both")
        if len(cols) > 1:
            ax.legend(frameon=False)
    fig.suptitle("aniso_zone_N: tratteggio = aniso_const")
    plt.tight_layout()
    plt.savefig(out_dir / "zone.png", dpi=140)
    plt.close()


def roney():
    rows, per_sub = [], []
    for sc_dir in sorted(RONEY.glob("aniso_*")):
        sc = sc_dir.name
        runs = {a: {d.name: load(d / "errors.csv") for d in (sc_dir / a).glob("sub_*")} for a in ARMS}
        d = None
        if sc != "aniso_const":
            area = next(iter(runs["faccia"].values()))["area"].sum()
            d = np.sqrt(4 * area / (np.pi * int(sc.removeprefix("aniso_zone_"))))
        for arm in ARMS:
            paired, nonstim_area, both_area, tot = [], 0.0, 0.0, 0.0
            for sub, df in runs[arm].items():
                both = (runs["roney"][sub]["estimable"] == 1) & (runs["faccia"][sub]["estimable"] == 1)
                paired.append(df[both])
                nonstim_area += nonstim(df) * df["area"].sum()
                both_area += df["area"][both].sum()
                tot += df["area"].sum()
                per_sub.append({"scenario": sc, "braccio": arm, "sub": sub,
                                "area_non_stimabile_%": nonstim(df), **quantiles(df[both])})
            rows.append({"scenario": sc, "braccio": arm, "d": d, "area_non_stimabile_%": nonstim_area / tot,
                         "area_appaiata_%": 100 * both_area / tot, **quantiles(pd.concat(paired))})
    r, s = pd.DataFrame(rows), pd.DataFrame(per_sub)
    out_dir = RONEY / "figure"
    out_dir.mkdir(parents=True, exist_ok=True)
    r.to_csv(out_dir / "numeri.csv", index=False, float_format="%.6g")
    s.to_csv(out_dir / "terne.csv", index=False, float_format="%.6g")

    panels = [("area_non_stimabile_%", "area non stimabile [%]")] + \
             [(e, label + (", valore assoluto" if e != "angolo" else "")) for e, label in ERRORS.items()]

    # aniso_const: un punto per terna, sfalsato attorno al braccio
    c = s[s["scenario"] == "aniso_const"].copy()
    c["x"] = c["braccio"].map({a: i for i, a in enumerate(ARMS)}) \
        + 0.6 * (c.groupby("braccio").cumcount() / 19 - 0.5)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    for ax, (col, ylabel) in zip(axes, panels):
        if col == "area_non_stimabile_%":
            ax.scatter(c["x"], c[col], s=24, color="#2a78d6", zorder=3)
        else:
            ax.scatter(c["x"], c[f"{col}_q50"], marker="o", s=24, color="#2a78d6", label="q50", zorder=3)
            ax.scatter(c["x"], c[f"{col}_q90"], marker="^", s=24, color="#eb6834", label="q90", zorder=3)
            ax.legend(frameon=False)
        ax.set_xticks(range(len(ARMS)), list(ARMS))
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("aniso_const, 3 mappe: un punto per terna, facce stimabili nei due bracci")
    plt.tight_layout()
    plt.savefig(out_dir / "const.png", dpi=140)
    plt.close()

    z = r[r["scenario"] != "aniso_const"].sort_values("d")
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    for ax, (col, ylabel) in zip(axes, panels):
        for arm, color in ARMS.items():
            a = z[z["braccio"] == arm]
            if col == "area_non_stimabile_%":
                ax.plot(a["d"], a[col], marker="o", color=color, lw=1, label=arm)
            else:
                ax.plot(a["d"], a[f"{col}_q50"], marker="o", color=color, lw=1, label=f"{arm} q50")
                ax.plot(a["d"], a[f"{col}_q90"], marker="^", color=color, lw=1, ls="--", label=f"{arm} q90")
        ax.set_xscale("log")
        ax.set_xlabel("diametro medio di zona d [mm]")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3, which="both")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("aniso_zone_N, 3 mappe: terne messe insieme, facce stimabili nei due bracci")
    plt.tight_layout(rect=(0, 0.08, 1, 1))
    plt.savefig(out_dir / "zone.png", dpi=140)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment", nargs="?", help="nome della cartella in output/")
    ap.add_argument("--sottoinsiemi", action="store_true", help="confronto fra output/subsets/sub_*")
    ap.add_argument("--zone", action="store_true", help="confronto fra output/zone/aniso_zone_*")
    ap.add_argument("--roney", action="store_true", help="confronto roney / faccia in output/roney/")
    args = ap.parse_args()
    if args.sottoinsiemi:
        subsets()
    elif args.zone:
        zone()
    elif args.roney:
        roney()
    elif args.experiment:
        experiment(args.experiment)
    else:
        ap.error("serve ESPERIMENTO, --sottoinsiemi, --zone oppure --roney")


if __name__ == "__main__":
    main()
