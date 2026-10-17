#!/usr/bin/env python3
"""
Analisi del CSV di errori prodotto dallo stadio 6.

Quattro viste, tutte dal solo CSV di un esperimento (output/<experiment_id>/errors.csv):
  1. errore per regione anatomica
  2. bias su CV_L: sistematico o regionale?
  3. le facce non stimabili: dove stanno
  4. errore in funzione del condizionamento

Tutte le statistiche sono PESATE PER AREA, non per conteggio di facce:
facce di dimensione diversa non contano uguale.

I percorsi si configurano nelle costanti in cima al file.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- percorsi
ROOT     = Path(__file__).resolve().parents[1]            # radice del repository
CSV_PATH = ROOT / "output" / "aniso_const" / "errors.csv"
FIG_DIR  = ROOT / "output" / "aniso_const" / "figures"
# --------------------------------------------------------------------------

REGION_NAMES = {
    1: "corpo LA",
    2: "anello mitralico",
    3: "LIPV",
    4: "LSPV",
    5: "RIPV",
    6: "RSPV",
    7: "LAA",
}


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


def main():
    df = pd.read_csv(CSV_PATH, dtype={"experiment_id": str})
    ids = df["experiment_id"].unique()
    if len(ids) != 1:
        # face_id non e' univoco fra esperimenti: mescolarli darebbe statistiche senza senso
        raise SystemExit(f"{CSV_PATH}: attesi dati di un solo esperimento, trovati {list(ids)}")
    exp_id = ids[0]
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df["region"] = df["region_tag"].map(REGION_NAMES).fillna(
        df["region_tag"].astype(str))
    est = df[df["estimable"] == 1].copy()

    area_tot = df["area"].sum()

    rule(f"PANORAMICA  -  esperimento {exp_id}")
    print(f"facce totali       : {len(df)}")
    print(f"stimabili          : {len(est)}  "
          f"({100 * est['area'].sum() / area_tot:.2f}% dell'area)")
    print(f"non stimabili      : {len(df) - len(est)}")
    print(f"area totale        : {area_tot:.1f} mm^2")

    a = est["area"].values
    print(f"\nerrore angolare (gradi, pesato per area)")
    for q in (0.25, 0.50, 0.75, 0.90, 0.95, 0.99):
        print(f"   q{int(q * 100):02d} : {wquantile(est['err_angle'], a, q):8.3f}")
    frac20 = est.loc[est["err_angle"] < 20, "area"].sum() / est["area"].sum()
    print(f"   area entro 20 gradi : {100 * frac20:.2f}%")
    print("   (Roney 2019, Tab.3, 3 mappe: mediana 10.7-12.1 deg, ~70% entro 20)")

    # ------------------------------------------------------- 1. per regione
    rule("1. ERRORE PER REGIONE ANATOMICA")
    rows = []
    for tag, g in est.groupby("region_tag"):
        w = g["area"].values
        n_tot = (df["region_tag"] == tag).sum()
        rows.append({
            "tag": tag,
            "regione": REGION_NAMES.get(tag, str(tag)),
            "facce": n_tot,
            "non_stim": n_tot - len(g),
            "ang_q50": wquantile(g["err_angle"], w, 0.50),
            "ang_q90": wquantile(g["err_angle"], w, 0.90),
            "CVL_q50": wquantile(g["err_CV_L_rel"], w, 0.50),
            "CVT_q50": wquantile(g["err_CV_T_rel"], w, 0.50),
            "cond_q50": wquantile(g["cond"], w, 0.50),
        })
    reg = pd.DataFrame(rows).sort_values("ang_q50", ascending=False)
    print(reg.to_string(index=False, float_format=lambda x: f"{x:9.4f}"))
    print("\nAttesa: LAA e vene peggio del corpo (curvatura alta, fibre "
          "disomogenee).\nSe l'errore non segue l'anatomia, c'e' altro da capire.")

    fig, ax = plt.subplots(figsize=(9, 5))
    order = reg["tag"].tolist()
    ax.boxplot([est.loc[est["region_tag"] == t, "err_angle"].dropna()
                for t in order],
               showfliers=False)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([REGION_NAMES.get(t, str(t)) for t in order])
    ax.set_ylabel("errore angolare di fibra [gradi]")
    ax.set_title(f"{exp_id}: errore angolare per regione (box non pesati, outlier nascosti)")
    ax.grid(axis="y", alpha=0.3)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "err_per_regione.png", dpi=140)
    plt.close()

    # ----------------------------------------------------------- 2. il bias
    rule("2. BIAS SU CV_L  -  sistematico o regionale?")
    print("Un bias costante in tutte le regioni indica errore di "
          "discretizzazione\ndel gradiente (l'interpolante lineare taglia la "
          "curvatura del fronte e\nsottostima la pendenza). Uno regionale o "
          "bimodale indica altro.\n")
    for tag, g in est.groupby("region_tag"):
        w = g["area"].values
        q25 = wquantile(g["err_CV_L_rel"], w, 0.25)
        q50 = wquantile(g["err_CV_L_rel"], w, 0.50)
        q75 = wquantile(g["err_CV_L_rel"], w, 0.75)
        print(f"  {REGION_NAMES.get(tag, tag):18s} "
              f"q25 {q25:+8.4f}   q50 {q50:+8.4f}   q75 {q75:+8.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, col, name in zip(axes,
                             ["err_CV_L_rel", "err_CV_T_rel"],
                             ["CV_L", "CV_T"]):
        v = est[col].dropna()
        lo, hi = np.percentile(v, [0.5, 99.5])
        ax.hist(v.clip(lo, hi), bins=120, weights=est.loc[v.index, "area"])
        ax.axvline(0, color="k", lw=1)
        ax.axvline(wquantile(v, est.loc[v.index, "area"], 0.5),
                   color="C3", lw=1.5, ls="--", label="mediana")
        ax.set_xlabel(f"errore relativo {name}")
        ax.set_ylabel("area [mm^2]")
        ax.legend()
        ax.grid(alpha=0.3)
    fig.suptitle(f"{exp_id}: distribuzione degli errori relativi sulle velocita'")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "bias_cv.png", dpi=140)
    plt.close()

    # ------------------------------------------------ 3. le non stimabili
    rule("3. FACCE NON STIMABILI")
    bad = df[df["estimable"] == 0]
    if len(bad) == 0:
        print("nessuna")
    else:
        print(f"{len(bad)} facce, {100 * bad['area'].sum() / area_tot:.3f}% "
              f"dell'area\n")
        print("distribuzione per regione:")
        for tag, g in bad.groupby("region_tag"):
            n_tot = (df["region_tag"] == tag).sum()
            print(f"  {REGION_NAMES.get(tag, tag):18s} {len(g):5d} "
                  f"su {n_tot:6d}  ({100 * len(g) / n_tot:.2f}%)")
        print("\nSe si concentrano vicino ai siti di stimolo e sulle linee di\n"
              "collisione, i criteri degeneri catturano gia' le zone giuste.\n"
              "Se sono sparse a caso, c'e' altro da capire.")
        bad[["face_id", "region_tag", "area"]].to_csv(
            FIG_DIR / "facce_non_stimabili.csv", index=False)
        print(f"\nelenco scritto in {FIG_DIR / 'facce_non_stimabili.csv'}\n"
              "(caricalo in ParaView come selezione per vedere dove stanno)")

    # ------------------------------------------- 4. errore vs condizionamento
    rule("4. ERRORE IN FUNZIONE DEL CONDIZIONAMENTO")
    print("Questa e' la relazione fondamentale: il numero di mappe e' solo un\n"
          "proxy grezzo della copertura direzionale, che cond misura "
          "direttamente.\n")
    c = est["cond"].replace([np.inf, -np.inf], np.nan).dropna()
    print("condizionamento (pesato per area)")
    for q in (0.50, 0.90, 0.99, 1.00):
        print(f"   q{int(q * 100):03d} : "
              f"{wquantile(c, est.loc[c.index, 'area'], q):12.2f}")

    bins = np.logspace(np.log10(max(c.min(), 1.0)), np.log10(c.max()), 13)
    idx = np.digitize(est["cond"], bins)
    print("\n  bin di cond            facce    ang_q50    ang_q90   CVL_q50")
    for b in range(1, len(bins)):
        g = est[idx == b]
        if len(g) < 20:
            continue
        w = g["area"].values
        print(f"  [{bins[b-1]:8.1f},{bins[b]:9.1f})  {len(g):7d}  "
              f"{wquantile(g['err_angle'], w, 0.50):9.3f}  "
              f"{wquantile(g['err_angle'], w, 0.90):9.3f}  "
              f"{wquantile(g['err_CV_L_rel'], w, 0.50):+9.4f}")

    fig, ax = plt.subplots(figsize=(9, 5))
    sub = est.sample(min(20000, len(est)), random_state=0)
    ax.scatter(sub["cond"], sub["err_angle"], s=2, alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("numero di condizionamento")
    ax.set_ylabel("errore angolare di fibra [gradi]")
    ax.set_title(f"{exp_id}: errore vs copertura direzionale")
    ax.grid(alpha=0.3, which="both")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "err_vs_cond.png", dpi=140)
    plt.close()

    rule("FIGURE")
    for f in sorted(FIG_DIR.glob("*.png")):
        print(f"  {f}")


if __name__ == "__main__":
    main()