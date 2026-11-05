#!/usr/bin/env python3
"""
Analisi della griglia di esperimenti: per ogni scenario aniso_zone_N i 42
sottoinsiemi di almeno 3 mappe con lo stimatore per faccia
(output/SCENARIO/faccia/sub_*/errors.csv) e le 20 terne con il metodo 2 di
Roney et al. 2019 (output/SCENARIO/roney/sub_i-j-k/errors.csv).

Uso, dalla radice del repository:
    python analysis/analisi.py
Scrive tutto in output/figure/.

Errori per faccia stimabile, statistiche PESATE PER AREA:
  angolo    err_angle, gradi
  cv_l      100 (CV_L_est - CV_L_true) / CV_L_true, %
  rapporto  100 (r_est - r_true) / r_true, %, r = CV_L / CV_T
q50, q90 e q99 si calcolano sul valore assoluto; per cv_l e rapporto c'e'
anche la mediana con segno. d = sqrt(4 A / (pi N)) e' il diametro medio di
zona in mm, A l'area totale.

numeri.csv, una riga per scenario e sottoinsieme del metodo per faccia: area
non stimabile ed errori sulle facce stimabili del sottoinsieme.

confronto.csv, una riga per scenario, terna e metodo (faccia, roney): area
non stimabile del metodo da solo, area appaiata (stimabile da entrambi i
metodi con quella terna) ed errori sulle facce appaiate. Le righe con sub =
tutte mettono insieme le facce appaiate delle 20 terne dello scenario; le
loro aree sono le medie sulle terne.

Figure, ognuna in quattro PNG: _area (area non stimabile, un pannello) e
_angolo, _cv_l, _rapporto. Nei grafici 1 e 2 questi ultimi hanno un solo
pannello con q50, q90 e q99 (cerchio, quadrato, triangolo, colori Okabe-Ito),
solo punti senza linee; nei grafici 3 e 4 tre pannelli q50, q90, q99, ognuno
con la sua scala. Scale lineari da 0.
  grafico1_mappe       aniso_zone_1, in funzione di k: un punto per
                       sottoinsieme e quantile, sfalsati dentro ogni k
  grafico2_zone        tutte le mappe, in funzione di d: un punto per
                       scenario e quantile
  grafico3_zone_mappe  in funzione di d, per ogni k la mediana dei
                       sottoinsiemi e la fascia dal minimo al massimo
  grafico4_roney       un punto per terna e scenario, Roney in orizzontale e
                       per faccia in verticale, stessi limiti sui due assi,
                       un colore e un simbolo per scenario
Mediana, minimo e massimo fra sottoinsiemi contano ogni sottoinsieme uno.
Nei grafici 2 e 3 la linea tratteggiata e' il diametro del disco di Roney,
2 patch.radius_mm di experiments/roney.conf.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
OUT = ROOT / "output"
FIG = OUT / "figure"
SCENARIO_MAPPE = "aniso_zone_1"                           # grafico 1
ERRORS = ("angolo", "cv_l", "rapporto")
QUANTILES = (50, 90, 99)
YLABEL = {"area": "area non stimabile [%]",
          "angolo": "errore angolare [gradi]",
          "cv_l": "|errore relativo su CV_L| [%]",
          "rapporto": "|errore relativo su CV_L / CV_T| [%]"}
METODI = {"faccia": "per faccia", "roney": "Roney"}
BLU = "#2a78d6"
QSTILE = {f"q{q}": (c, m) for q, c, m in zip(QUANTILES, ("#0072B2", "#E69F00", "#009E73"), "os^")}  # Okabe-Ito
SIMBOLI = "osD^v"                                         # scenari del grafico 4
RONEY_DISCO = 2 * float(next(l.split("=")[1] for l in open(ROOT / "experiments" / "roney.conf")
                             if l.startswith("patch.radius_mm")))


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


def quantiles(df):
    """q50, q90 e q99 del valore assoluto di ogni errore sulle facce stimabili, poi la mediana con
    segno di cv_l e rapporto."""
    est = df[df["estimable"] == 1]
    w = est["area"]
    out = {f"{e}_q{q}": wquantile(est[e].abs(), w, q / 100) for e in ERRORS for q in QUANTILES}
    for e in ("cv_l", "rapporto"):
        out[f"{e}_mediana"] = wquantile(est[e], w, 0.50)
    return out


def tables():
    """numeri.csv e confronto.csv, come DataFrame."""
    numeri, confronto = [], []
    for sc_dir in sorted(OUT.glob("aniso_zone_*"), key=lambda p: int(p.name.removeprefix("aniso_zone_"))):
        sc = sc_dir.name
        n = int(sc.removeprefix("aniso_zone_"))
        faccia = {d.name: load(d / "errors.csv") for d in (sc_dir / "faccia").glob("sub_*")}
        roney = {d.name: load(d / "errors.csv") for d in (sc_dir / "roney").glob("sub_*")}
        area = next(iter(faccia.values()))["area"]
        head = {"scenario": sc, "N": n, "d": np.sqrt(4 * area.sum() / (np.pi * n))}

        for sub in sorted(faccia, key=lambda s: (s.count("-"), s)):
            df = faccia[sub]
            numeri.append({**head, "sub": sub, "k": sub.count("-") + 1,
                           "area_non_stimabile_%": nonstim(df), **quantiles(df)})

        for metodo, runs in (("faccia", faccia), ("roney", roney)):
            rows, paired = [], []
            for sub in sorted(roney):
                df = runs[sub]
                both = (roney[sub]["estimable"] == 1) & (faccia[sub]["estimable"] == 1)
                rows.append({**head, "sub": sub, "metodo": metodo, "area_non_stimabile_%": nonstim(df),
                             "area_appaiata_%": 100 * area[both].sum() / area.sum(), **quantiles(df[both])})
                paired.append(df[both])
            confronto += rows
            confronto.append({**head, "sub": "tutte", "metodo": metodo,
                              "area_non_stimabile_%": np.mean([r["area_non_stimabile_%"] for r in rows]),
                              "area_appaiata_%": np.mean([r["area_appaiata_%"] for r in rows]),
                              **quantiles(pd.concat(paired))})
    return pd.DataFrame(numeri), pd.DataFrame(confronto)


def figure(name, title, draw, insieme=False):
    """Quattro PNG, output/figure/NAME_{area,angolo,cv_l,rapporto}.png. draw(ax, cols, ylabel) disegna
    in un pannello le colonne di cols, {etichetta: colonna}. Con insieme q50, q90 e q99 stanno in un
    solo pannello, altrimenti un pannello per quantile."""
    for q, ylabel in YLABEL.items():
        if q == "area":
            panels = {"area non stimabile": {"area non stimabile": "area_non_stimabile_%"}}
        elif insieme:
            panels = {"q50, q90 e q99": {f"q{p}": f"{q}_q{p}" for p in QUANTILES}}
        else:
            panels = {f"q{p}": {f"q{p}": f"{q}_q{p}"} for p in QUANTILES}
        fig, axes = plt.subplots(1, len(panels), figsize=(max(4.6 * len(panels), 7.5), 4.6), squeeze=False,
                                 layout="constrained")
        for ax, (panel, cols) in zip(axes[0], panels.items()):
            draw(ax, cols, ylabel)
            ax.set_title(panel)
            ax.grid(alpha=0.3)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=min(len(labels), 5 if len(panels) > 1 else 4),
                   frameon=False)
        fig.suptitle(title)
        fig.savefig(FIG / f"{name}_{q}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)


def asse_y(ax, ylabel):
    """Asse verticale lineare, da 0."""
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)


def punti(ax, x, df, cols):
    """Un punto per riga e colonna di cols, colore e simbolo del quantile."""
    for label, col in cols.items():
        color, marker = QSTILE.get(label, (BLU, "o"))
        ax.scatter(x, df[col], s=24, color=color, marker=marker, label=label, zorder=3)


def asse_d(ax):
    """Asse orizzontale d in scala log, con la linea del diametro del disco di Roney."""
    ax.axvline(RONEY_DISCO, color="#555", lw=1, ls="--", label=f"diametro del disco di Roney ({RONEY_DISCO:g} mm)")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.set_xlabel("diametro medio di zona d [mm]")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    numeri, confronto = tables()
    numeri.to_csv(FIG / "numeri.csv", index=False, float_format="%.6g")
    confronto.to_csv(FIG / "confronto.csv", index=False, float_format="%.6g")

    # grafico 1: aniso_zone_1, in funzione di k, punti sfalsati dentro ogni k
    g1 = numeri[numeri["scenario"] == SCENARIO_MAPPE].copy()
    i, n = g1.groupby("k").cumcount(), g1.groupby("k")["k"].transform("size")
    g1["x"] = g1["k"] + 0.6 * (i - (n - 1) / 2) / n

    def draw1(ax, cols, ylabel):
        punti(ax, g1["x"], g1, cols)
        ax.set_xticks(sorted(g1["k"].unique()))
        ax.set_xlabel("mappe usate k")
        asse_y(ax, ylabel)

    figure("grafico1_mappe", f"{SCENARIO_MAPPE}: un punto per sottoinsieme di mappe", draw1, insieme=True)

    # grafico 2: tutte le mappe, in funzione di d
    g2 = numeri[numeri["k"] == numeri["k"].max()].sort_values("d")

    def draw2(ax, cols, ylabel):
        punti(ax, g2["d"], g2, cols)
        asse_d(ax)
        asse_y(ax, ylabel)

    figure("grafico2_zone", f"aniso_zone_N, per faccia, {g2['k'].iloc[0]} mappe: un punto per scenario", draw2,
           insieme=True)

    # grafico 3: in funzione di d, una linea e una fascia per k
    ks = sorted(numeri["k"].unique())
    colors = dict(zip(ks, plt.cm.viridis(np.linspace(0, 0.85, len(ks)))))

    def draw3(ax, cols, ylabel):
        (col,) = cols.values()
        for k in ks:
            s = numeri[numeri["k"] == k].groupby("d")[col]
            ax.fill_between(s.min().index, s.min(), s.max(), color=colors[k], alpha=0.15, lw=0)
            ax.plot(s.median().index, s.median(), color=colors[k], lw=2, marker="o", ms=4, label=f"k = {k}")
        asse_d(ax)
        asse_y(ax, ylabel)

    figure("grafico3_zone_mappe", "aniso_zone_N: mediana dei sottoinsiemi di k mappe, fascia dal minimo al massimo",
           draw3)

    # grafico 4: per faccia contro Roney, un punto per terna e scenario
    terne = confronto[confronto["sub"] != "tutte"]
    g4 = terne[terne["metodo"] == "roney"].merge(terne[terne["metodo"] == "faccia"],
                                                 on=["scenario", "N", "d", "sub"], suffixes=("_roney", "_faccia"))
    scenari = g4.drop_duplicates("scenario").sort_values("d")
    colors4 = dict(zip(scenari["scenario"], plt.cm.viridis(np.linspace(0, 0.85, len(scenari)))))

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for (sc, n, d), marker in zip(scenari[["scenario", "N", "d"]].itertuples(index=False), SIMBOLI):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_faccia"], s=22, color=colors4[sc], marker=marker,
                       label=f"N = {n} (d = {d:.0f} mm)" if d >= 10 else f"N = {n} (d = {d:.2g} mm)", zorder=3)
        v = pd.concat([g4[f"{col}_roney"], g4[f"{col}_faccia"]])
        lim = (0, 1.05 * v.max())
        ax.plot(lim, lim, color="#555", lw=1, label="y = x")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{METODI['roney']}: {ylabel}")
        ax.set_ylabel(f"{METODI['faccia']}: {ylabel}")

    figure("grafico4_roney", "terne di mappe: per faccia contro Roney, facce appaiate (area: ogni metodo da solo)",
           draw4)


if __name__ == "__main__":
    main()
