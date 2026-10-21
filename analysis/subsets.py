#!/usr/bin/env python3
"""
Passo 2: confronto fra i sottoinsiemi di mappe (k = 3..6).

Legge output/subsets/sub_*/errors.csv (generati da experiments/subsets.py)
e risponde a: quante mappe servono, e quali siti contano.

Confronto APPAIATO contro k = 6: per ogni sottoinsieme S, con F_S le sue
facce stimabili, l'errore di S su F_S si confronta con l'errore di k = 6
sulle stesse facce (dq50 = q50 - q50_ref). Con meno mappe le facce
difficili diventano non stimabili e l'errore sulle sole stimabili
sembrerebbe migliore di quanto e': il confronto appaiato toglie l'effetto
del cambio di facce, e ogni numero di errore sta accanto all'area non
stimabile, separata per causa.

Tutte le statistiche sono PESATE PER AREA, con wquantile di analisi.py.

Uso, dalla radice del repository:
    python analysis/subsets.py
Tabella in output/subsets/summary.csv, figure in output/subsets/figures/.
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analisi import MIN_GRADIENTS, ROOT, rule, wquantile

OUT_DIR = ROOT / "output" / "subsets"
FIG_DIR = OUT_DIR / "figures"
K_COLORS = {3: "#86b6ef", 4: "#3987e5", 5: "#1c5cab", 6: "#0d366b"}   # rampa ordinale, un solo tono
POINT_COLOR = "#2a78d6"
CAUSE_COLORS = {"nonstim_mask_%": "#2a78d6", "nonstim_sys_%": "#eb6834"}
CAUSE_LABELS = {"nonstim_mask_%": "k_used < 3 (maschere)",
                "nonstim_sys_%": "k_used >= 3 (sigma_min o lambda_2)"}


def load():
    """Per sottoinsieme: (indici, nomi dei siti, errors.csv). Ordine (k, indici)."""
    runs = []
    for d in OUT_DIR.glob("sub_*"):
        idx = tuple(int(i) for i in d.name[4:].split("-"))
        names = list(dict.fromkeys(pd.read_csv(d / "mask_pairs.csv", usecols=["map_name"],
                                               dtype=str)["map_name"]))
        df = pd.read_csv(d / "errors.csv", usecols=["face_id", "area", "estimable", "k_used",
                                                    "cond", "err_angle"]).set_index("face_id")
        df["cond"] = df["cond"].replace([np.inf, -np.inf], np.nan)
        runs.append((idx, names, df))
    return sorted(runs, key=lambda r: (len(r[0]), r[0]))


def strip(ax, s, col, ylabel, title):
    """Un punto per sottoinsieme per k, jitter deterministico; etichetta su migliore e peggiore."""
    for k, g in s.groupby("k"):
        g = g.sort_values(col)
        x = k + np.linspace(-0.18, 0.18, len(g)) if len(g) > 1 else np.array([float(k)])
        ax.scatter(x, g[col], s=36, color=POINT_COLOR, edgecolor="white", lw=1, zorder=3)
        if len(g) > 1:
            for xi, r, va in ((x[0], g.iloc[0], "top"), (x[-1], g.iloc[-1], "bottom")):
                ax.annotate(r["siti"], (xi, r[col]), xytext=(6, -4 if va == "top" else 4),
                            textcoords="offset points", fontsize=7.5, color="#444", va=va)
    ax.axhline(0, color="#888", lw=1, ls="--", zorder=1)
    ax.set_xticks(sorted(s["k"].unique()))
    ax.set_xlabel("mappe usate (k)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)


def main():
    runs = load()
    ref = runs[-1][2]                                     # sub_0-1-2-3-4-5
    area = ref["area"]
    area_tot = area.sum()
    pct = lambda m: 100 * area[m].sum() / area_tot

    rows, n_out = [], {}
    for idx, names, df in runs:
        F = df["estimable"] == 1
        F_ref = F & (ref["estimable"] == 1)               # F_S senza le non stimabili a k = 6: S e k = 6 su queste facce
        exp_id = "sub_" + "-".join(map(str, idx))
        n_out[exp_id] = (F & ~F_ref).sum(), pct(F & ~F_ref)
        est, s_ref, r = df[F], df[F_ref], ref[F_ref]
        q50, q90 = (wquantile(s_ref["err_angle"], s_ref["area"], q) for q in (0.50, 0.90))
        q50_ref, q90_ref = (wquantile(r["err_angle"], r["area"], q) for q in (0.50, 0.90))
        rows.append({
            "experiment_id": exp_id,
            "k": len(idx),
            "siti": "+".join(names),
            "nonstim_%": pct(~F),
            "nonstim_mask_%": pct(df["k_used"] < MIN_GRADIENTS),
            "nonstim_sys_%": pct(~F & (df["k_used"] >= MIN_GRADIENTS)),
            "q50": q50, "q90": q90,
            "q50_ref": q50_ref, "q90_ref": q90_ref,
            "dq50": q50 - q50_ref, "dq90": q90 - q90_ref,
            "cond_q50": wquantile(est["cond"], est["area"], 0.50),
        })
    s = pd.DataFrame(rows)
    s.to_csv(OUT_DIR / "summary.csv", index=False, float_format="%.6g")
    fmt = lambda x: f"{x:8.3f}"
    full = s.iloc[-1]

    rule(f"CONFRONTO APPAIATO CONTRO k = 6  -  {len(s)} sottoinsiemi")
    print(f"k = 6: q50 {full['q50']:.3f}   q90 {full['q90']:.3f}   "
          f"non stim. {full['nonstim_%']:.2f}%   dq50 {full['dq50']:g}   dq90 {full['dq90']:g}")
    n = np.array([v[0] for v in n_out.values()])
    a = np.array([v[1] for v in n_out.values()])
    print(f"facce di F_S non stimabili a k = 6 (escluse dal confronto): "
          f"min {n.min()}  mediana {np.median(n):.0f}  max {n.max()} facce, "
          f"max {a.max():.3f}% dell'area")

    rule("a) PER k: dq50 = q50(S) - q50(k = 6) sulle facce stimabili in S e a k = 6, e area non stimabile")
    for k, g in s.groupby("k"):
        print(f"\nk = {k}  ({len(g)} sottoinsiemi)")
        for col in ("dq50", "nonstim_%"):
            print(f"   {col:10s} min {g[col].min():8.3f}   mediana {g[col].median():8.3f}   "
                  f"max {g[col].max():8.3f}")
        if k == 3:
            print(f"   mediana nonstim_mask_% {g['nonstim_mask_%'].median():.3f}   "
                  f"nonstim_sys_% {g['nonstim_sys_%'].median():.3f}")
        for lab, r in (("migliore", g.loc[g["dq50"].idxmin()]),
                       ("peggiore", g.loc[g["dq50"].idxmax()])):
            print(f"   {lab}: {r['siti']:28s} dq50 {r['dq50']:+7.3f}   dq90 {r['dq90']:+7.3f}   "
                  f"non stim. {r['nonstim_%']:6.2f}%")

    rule("b) LEAVE-ONE-OUT: k = 5, ordinati per dq50")
    all_names = full["siti"].split("+")
    loo = s[s["k"] == len(all_names) - 1].sort_values("dq50").copy()
    loo.insert(0, "tolto", [next(n for n in all_names if n not in r.split("+")) for r in loo["siti"]])
    print(loo[["tolto", "dq50", "dq90", "nonstim_%", "nonstim_mask_%", "nonstim_sys_%"]]
          .to_string(index=False, float_format=fmt))

    rule("c) SPEARMAN tra cond_q50 e dq50, per k")
    for k, g in s.groupby("k"):
        if len(g) < 3:
            continue
        rho, p = spearmanr(g["cond_q50"], g["dq50"])
        print(f"k = {k}: rho = {rho:+.3f}   p = {p:.2g}   (n = {len(g)})")

    rule("d) AREA NON STIMABILE IN ALMENO UN SOTTOINSIEME, per k")
    for k in sorted(s["k"].unique()):
        bad = np.logical_or.reduce([df["estimable"].values == 0 for idx, _, df in runs if len(idx) == k])
        print(f"k = {k}: {pct(bad):6.2f}%")

    # ------------------------------------------------------------- figure
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    strip(ax, s, "dq50", "dq50 = q50(S) - q50(k = 6) [gradi]",
          "errore angolare appaiato contro k = 6, un punto per sottoinsieme")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "dq50_per_k.png", dpi=140)
    plt.close()

    # area non stimabile: una barra impilata per sottoinsieme, raggruppate per k
    fig, ax = plt.subplots(figsize=(13, 5.5))
    x = np.arange(len(s)) + s["k"].values - s["k"].min()   # uno spazio fra i gruppi di k
    bottom = np.zeros(len(s))
    for col in CAUSE_COLORS:
        ax.bar(x, s[col], bottom=bottom, width=0.8, color=CAUSE_COLORS[col],
               edgecolor="white", lw=1, label=CAUSE_LABELS[col])
        bottom += s[col].values
    ax.set_xticks(x)
    ax.set_xticklabels(s["siti"], rotation=90, fontsize=6.5)
    for k, g in s.groupby("k"):
        xs = x[g.index]
        ax.text(xs.mean(), 1.0, f"k = {k}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("area non stimabile [%]")
    ax.set_title("area non stimabile per sottoinsieme, per causa", pad=18)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "nonstim_per_k.png", dpi=140)
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for k, g in s.groupby("k"):
        ax.scatter(g["cond_q50"], g["dq50"], s=40, color=K_COLORS[k],
                   edgecolor="white", lw=1, label=f"k = {k}", zorder=3)
    ax.axhline(0, color="#888", lw=1, ls="--", zorder=1)
    ax.set_xscale("log")
    ax.set_xlabel("cond mediano sulle facce stimabili di S")
    ax.set_ylabel("dq50 = q50(S) - q50(k = 6) [gradi]")
    ax.set_title("errore appaiato contro condizionamento, un punto per sottoinsieme")
    ax.grid(alpha=0.3, which="both")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "dq50_vs_cond.png", dpi=140)
    plt.close()

    rule("FIGURE")
    for f in sorted(FIG_DIR.glob("*.png")):
        print(f"  {f}")


if __name__ == "__main__":
    main()
