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
anche la mediana con segno. cond = sigma_max / sigma_min del sistema di una
faccia; i suoi quantili, pesati per area, sono sulle facce con cond calcolato
(stimabili e scartate per lambda2 <= 0). d = sqrt(4 A / (pi N)) e' il diametro medio di
zona in mm, A l'area totale.

numeri.csv, una riga per scenario e sottoinsieme del metodo per faccia: area
non stimabile, errori sulle facce stimabili del sottoinsieme e q50, q90 e q99
del cond.

confronto.csv, una riga per scenario, terna e metodo (faccia, roney): area
non stimabile del metodo da solo, area appaiata (stimabile da entrambi i
metodi con quella terna) ed errori sulle facce appaiate. Le righe con sub =
tutte mettono insieme le facce appaiate delle 20 terne dello scenario; le
loro aree sono le medie sulle terne.

Scatole dei grafici 0a e 0b: dal q25 al q75, linea al q50, baffi dal q5 al q95,
pesati per area.
grafico0_cond_errori.png, aniso_zone_1, facce con cond calcolato dei 42
sottoinsiemi messe insieme, in 10 fasce di cond di area uguale, ognuna al suo
cond mediano. Tre pannelli angolo, cv_l, rapporto, una scatola per fascia sulle
sue facce stimabili, scala log sui due assi, nel titolo la correlazione di
Spearman fra cond ed errore, non pesata; in un quadrato 2 x 2 con il pannello
dell'area: la % di area della fascia dove la stima non e' un'ellisse, cioe' le
facce con cond ma con l'autovalore minore lambda2 <= 0 (la conica stimata e'
un'iperbole). Le facce scartate per valore singolare minimo sotto soglia non
hanno cond e non compaiono.
grafico0_mappe_cond.png, aniso_zone_1, una scatola del cond per k, sulle facce
con cond calcolato di tutti i sottoinsiemi con k mappe, scala log.
grafico0_cond_area.png, aniso_zone_1, un punto per sottoinsieme da numeri.csv:
cond_q90 (log) in orizzontale, area non stimabile in verticale, colore e
simbolo per k, nel titolo il rho di Spearman sui 42 punti.

Le altre figure, ognuna in quattro PNG: _area (area non stimabile, un pannello) e
_angolo, _cv_l, _rapporto. Nei grafici 1 e 2 questi ultimi hanno un solo
pannello con q50, q90 e q99 (cerchio, quadrato, triangolo, colori Okabe-Ito),
solo punti senza linee, con l'asse verticale in scala log ed etichette in
numeri semplici; nei grafici 3 e 4 tre pannelli q50, q90, q99, ognuno con la
sua scala. Gli altri assi verticali, compresi i file _area, sono lineari da 0;
nel grafico 4 anche l'asse orizzontale, con gli stessi limiti.
  grafico1_mappe       aniso_zone_1, in funzione di k: un punto per
                       sottoinsieme e quantile, sfalsati dentro ogni k
  grafico2_zone        tutte le mappe, in funzione di d: un punto per
                       scenario e quantile
  grafico3_zone_mappe  in funzione di d, per ogni k la mediana dei
                       sottoinsiemi e la fascia dal minimo al massimo
  grafico4_roney       un punto per terna e scenario, Roney in orizzontale e
                       per faccia in verticale, stessi limiti sui due assi,
                       un colore e un simbolo per scenario; sugli assi solo 0,
                       minimo, baricentro (media degli 80 punti, croce) e
                       massimo della nuvola, con proiezioni punteggiate e
                       senza griglia
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
from matplotlib.colors import to_rgba
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter
from scipy.stats import spearmanr

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
QSTILE = {f"q{q}": (c, m) for q, c, m in zip(QUANTILES, ("#0072B2", "#D55E00", "#009E73"), "os^")}  # Okabe-Ito
SIMBOLI = "osD^v"                                         # scenari del grafico 4
COLORI_N = {1: "#000000", 8: "#0072B2", 70: "#D55E00", 600: "#009E73"}   # grafico 4, per N, Okabe-Ito
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


def cond_quantiles(df):
    """q50, q90 e q99 del cond sulle facce con cond calcolato."""
    c = df[df["cond"].notna()]
    return {f"cond_q{q}": wquantile(c["cond"], c["area"], q / 100) for q in QUANTILES}


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
                           "area_non_stimabile_%": nonstim(df), **quantiles(df), **cond_quantiles(df)})

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


def figure(name, title, draw, insieme=False, dopo=None):
    """Quattro PNG, output/figure/NAME_{area,angolo,cv_l,rapporto}.png. draw(ax, cols, ylabel) disegna
    in un pannello le colonne di cols, {etichetta: colonna}. Con insieme q50, q90 e q99 stanno in un
    solo pannello, altrimenti un pannello per quantile. dopo(fig), se c'e', si chiama prima di salvare."""
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
            ax.grid(alpha=0.3)
            draw(ax, cols, ylabel)
            ax.set_title(panel)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=min(len(labels), 5 if len(panels) > 1 else 4),
                   frameon=False)
        fig.suptitle(title)
        if dopo:
            dopo(fig)
        fig.savefig(FIG / f"{name}_{q}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)


def scala_log(ax, assi):
    """Scala log sugli assi dati ("x", "y"), con etichette in numeri semplici (1, 2, 5, 10, 20) invece
    che in potenze di 10: tacche 1-2-5 per decade, oppure 1-9 se cosi' le etichette sono meno di quattro.
    Va chiamata dopo aver disegnato i dati."""
    for a in assi:
        getattr(ax, f"set_{a}scale")("log")
        axis = getattr(ax, f"{a}axis")
        lo, hi = getattr(ax, f"get_{a}lim")()
        for subs in ((1, 2, 5), (1, 2, 3, 4, 5, 6, 7, 8, 9)):
            ticks = [t for t in LogLocator(subs=subs).tick_values(lo, hi) if lo <= t <= hi]
            if len(ticks) >= 4:
                break
        axis.set_major_locator(LogLocator(subs=subs))
        axis.set_major_formatter(FormatStrFormatter("%g"))
        axis.set_minor_formatter(NullFormatter())
        ax.grid(axis=a, which="minor", alpha=0.15)


def separa_etichette(fig):
    """Sposta su piu' righe le etichette delle tacche che si sovrappongono: in basso sull'asse x, a
    sinistra sull'asse y."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_pt = 72 / fig.dpi
    for ax in fig.axes:
        for axis, lato in ((ax.xaxis, 0), (ax.yaxis, 1)):
            ticks = sorted(axis.get_major_ticks(), key=lambda t: t.get_loc())
            boxes = [t.label1.get_window_extent(renderer) for t in ticks]
            passo = max((b.height if lato == 0 else b.width) for b in boxes) + 3
            fine = []                                     # per riga, dove finisce l'ultima etichetta
            for t, b in zip(ticks, boxes):
                a, z = (b.x0, b.x1) if lato == 0 else (b.y0, b.y1)
                r = next((i for i, f in enumerate(fine) if a > f + 3), len(fine))
                fine[r:r + 1] = [z]
                t.set_pad(t.get_pad() + r * passo * px_pt)


def asse_y(ax, ylabel, log=False):
    """Asse verticale: in scala log, oppure lineare da 0."""
    if log:
        ax.set_ylabel(f"{ylabel} (scala log)")
        scala_log(ax, "y")
    else:
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


def scatole(ax, x, x0, x1, gruppi, col):
    """Una scatola per gruppo, in x da x0 a x1: dal q25 al q75 del valore assoluto di col, linea al q50 e baffi
    dal q5 al q95, pesati per area."""
    q5, q25, q50, q75, q95 = (gruppi.apply(lambda s: wquantile(s[col].abs(), s["area"], q)).to_numpy()
                              for q in (0.05, 0.25, 0.50, 0.75, 0.95))
    blu = QSTILE["q50"][0]
    ax.vlines(x, q5, q25, color=blu, lw=1)
    ax.vlines(x, q75, q95, color=blu, lw=1)
    ax.bar(x0, q75 - q25, width=x1 - x0, bottom=q25, align="edge", facecolor=to_rgba(blu, 0.25), edgecolor=blu,
           lw=1, zorder=3)
    ax.hlines(q50, x0, x1, color=blu, lw=2, zorder=4)


def nota_scatole(fig):
    """La legenda delle scatole, sotto la figura."""
    fig.supxlabel("scatola: q25–q75 · linea: q50 · baffi: q5–q95 (pesati per area)", fontsize="medium")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    numeri, confronto = tables()
    numeri.to_csv(FIG / "numeri.csv", index=False, float_format="%.6g")
    confronto.to_csv(FIG / "confronto.csv", index=False, float_format="%.6g")

    # grafico 0a: aniso_zone_1, tutti i sottoinsiemi, errore in funzione del cond
    df = pd.concat([load(d / "errors.csv").assign(k=d.name.count("-") + 1)
                    for d in (OUT / SCENARIO_MAPPE / "faccia").glob("sub_*")])
    c = df[df["cond"].notna()].sort_values("cond")
    w = c["area"]
    c["fascia"] = np.minimum((10 * (w.cumsum() - 0.5 * w) / w.sum()).astype(int), 9)   # 10 fasce di area uguale
    est = c[c["estimable"] == 1]
    xc = c.groupby("fascia").apply(lambda s: wquantile(s["cond"], s["area"], 0.50)).to_numpy()
    h = 10 ** (0.4 * np.diff(np.log10(xc)).min())          # mezza larghezza delle scatole, costante in log
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    axes = axes.flat
    for ax, e in zip(axes, ERRORS):
        scatole(ax, xc, xc / h, xc * h, est.groupby("fascia"), e)
        ax.grid(alpha=0.3)
        ax.set_xlabel("cond (scala log)")
        ax.set_ylabel(f"{YLABEL[e]} (scala log)")
        scala_log(ax, "xy")
        ax.set_title(f"{e}, ρ di Spearman = {spearmanr(est['cond'], est[e].abs()).statistic:.2f}")
    ax = axes[3]
    ax.plot(xc, c.groupby("fascia").apply(nonstim), color=QSTILE["q50"][0], lw=1.5, marker="o", ms=5)
    ax.grid(alpha=0.3)
    ax.set_xlabel("cond (scala log)")
    scala_log(ax, "x")
    asse_y(ax, "area dove la stima non è un'ellisse [%]")
    ax.set_title("stima che non è un'ellisse")
    nota_scatole(fig)
    fig.suptitle(f"{SCENARIO_MAPPE}, tutti i sottoinsiemi: errore in funzione del cond,\nfasce di area uguale")
    fig.savefig(FIG / "grafico0_cond_errori.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # grafico 0b: aniso_zone_1, una scatola del cond per k
    ks = np.array(sorted(c["k"].unique()))
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    scatole(ax, ks, ks - 0.25, ks + 0.25, c.groupby("k"), "cond")
    ax.grid(alpha=0.3)
    ax.set_xticks(ks)
    ax.set_xlabel("mappe usate k")
    asse_y(ax, "cond", log=True)
    nota_scatole(fig)
    fig.suptitle(f"{SCENARIO_MAPPE}: distribuzione del cond per numero di mappe")
    fig.savefig(FIG / "grafico0_mappe_cond.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # grafico 0c: aniso_zone_1, area non stimabile contro cond_q90, un punto per sottoinsieme
    g0 = numeri[numeri["scenario"] == SCENARIO_MAPPE]
    rho = spearmanr(g0["cond_q90"], g0["area_non_stimabile_%"]).statistic
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    for (k, s), (color, marker) in zip(g0.groupby("k"), (("#000000", "o"), ("#0072B2", "s"), ("#D55E00", "^"),
                                                         ("#009E73", "D"))):
        ax.scatter(s["cond_q90"], s["area_non_stimabile_%"], s=28, color=color, marker=marker, label=f"k = {k}",
                   zorder=3)
    ax.grid(alpha=0.3)
    ax.set_xlabel("cond, q90 del sottoinsieme (scala log)")
    scala_log(ax, "x")
    asse_y(ax, YLABEL["area"])
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle(f"{SCENARIO_MAPPE}: area non stimabile e coda del cond, un punto per sottoinsieme "
                 f"(ρ di Spearman = {rho:.2f})")
    fig.savefig(FIG / "grafico0_cond_area.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # grafico 1: aniso_zone_1, in funzione di k, punti sfalsati dentro ogni k
    g1 = numeri[numeri["scenario"] == SCENARIO_MAPPE].copy()
    i, n = g1.groupby("k").cumcount(), g1.groupby("k")["k"].transform("size")
    g1["x"] = g1["k"] + 0.6 * (i - (n - 1) / 2) / n


    def draw1(ax, cols, ylabel):
        punti(ax, g1["x"], g1, cols)
        ax.set_xticks(sorted(g1["k"].unique()))
        ax.set_xlabel("mappe usate k")
        asse_y(ax, ylabel, log=len(cols) > 1)

    figure("grafico1_mappe", f"{SCENARIO_MAPPE}: un punto per sottoinsieme di mappe", draw1, insieme=True)

    # grafico 2: tutte le mappe, in funzione di d
    g2 = numeri[numeri["k"] == numeri["k"].max()].sort_values("d")

    def draw2(ax, cols, ylabel):
        punti(ax, g2["d"], g2, cols)
        asse_d(ax)
        asse_y(ax, ylabel, log=len(cols) > 1)

    figure("grafico2_zone", f"aniso_zone_N, per faccia, {g2['k'].iloc[0]} mappe: un punto per scenario", draw2,
           insieme=True)

    # grafico 3: in funzione di d, una linea e una fascia per k
    ks = sorted(numeri["k"].unique())
    colors = dict(zip(ks, plt.cm.viridis(np.linspace(0, 0.55, len(ks)))))

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

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for (sc, n, d), marker in zip(scenari[["scenario", "N", "d"]].itertuples(index=False), SIMBOLI):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_faccia"], s=22, color=COLORI_N[n], marker=marker,
                       label=f"N = {n} (d = {d:.0f} mm)" if d >= 10 else f"N = {n} (d = {d:.2g} mm)", zorder=3)
        x, y = g4[f"{col}_roney"], g4[f"{col}_faccia"]
        lim = (0, 1.05 * max(x.max(), y.max()))
        ax.plot(lim, lim, color="#555", lw=1, label="y = x")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        # sugli assi solo 0, minimo, baricentro e massimo della nuvola, con le proiezioni
        mx, my = x.mean(), y.mean()
        proiezione = dict(color="#444", lw=0.8, ls=":", zorder=2)
        for i in (x.idxmin(), x.idxmax()):
            ax.plot([x[i], x[i]], [0, y[i]], **proiezione)
        for i in (y.idxmin(), y.idxmax()):
            ax.plot([0, x[i]], [y[i], y[i]], **proiezione)
        ax.plot([mx, mx, 0], [0, my, my], **proiezione)
        ax.scatter(mx, my, s=170, marker="+", linewidths=2, color="#444", label="baricentro", zorder=5)
        for setter, v in ((ax.set_xticks, x), (ax.set_yticks, y)):
            vals = [0, v.min(), v.mean(), v.max()]
            setter(vals, ["0"] + [f"{t:#.3g}".rstrip(".") for t in vals[1:]])   # 19.0, non 19
        ax.grid(False)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{METODI['roney']}: {ylabel}")
        ax.set_ylabel(f"{METODI['faccia']}: {ylabel}")

    figure("grafico4_roney", "terne di mappe: per faccia contro Roney, facce appaiate (area: ogni metodo da solo)\n"
           "sugli assi: min, baricentro (media delle terne) e max", draw4, dopo=separa_etichette)


if __name__ == "__main__":
    main()
