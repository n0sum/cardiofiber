#!/usr/bin/env python3
"""
Analisi del CSV di errori prodotto dallo stadio 6.

Cinque viste, dai CSV di un esperimento (output/<esperimento>/errors.csv e mask_pairs.csv):
  1. errore per regione anatomica
  2. bias su CV_L: sistematico o regionale?
  3. le facce non stimabili: dove stanno
  4. errore in funzione del condizionamento
  5. r = g^T M_true g - 1 sulle coppie (faccia, mappa) tenute contro escluse

Accanto a ogni numero di errore: frazione di area delle coppie escluse dalle
maschere (somma delle aree delle coppie escluse / (K * area del
sottoinsieme)), e frazione di area resa non stimabile (k_used < 3).

Tutte le statistiche sono PESATE PER AREA, non per conteggio di facce:
facce di dimensione diversa non contano uguale.

Uso, dalla radice del repository:
    python analysis/analisi.py aniso_const
Le figure vanno in output/<esperimento>/figures/.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]                # radice del repository

REGION_NAMES = {
    1: "corpo LA",
    2: "anello mitralico",
    3: "LIPV",
    4: "LSPV",
    5: "RIPV",
    6: "RSPV",
    7: "LAA",
}

# valori di masked (experiments/README.md)
MASK_LABELS = {0: "tenute", 1: "escluse"}
MASK_COLORS = {0: "#2a78d6", 1: "#eb6834"}
MIN_GRADIENTS = 3                                         # kMinGradients in Mesh.cpp


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


def mask_fracs(sub, K):
    """Frazioni di area del sottoinsieme di facce: coppie escluse e area non stimabile per maschere."""
    a = sub["area"]
    tot = a.sum()
    return {"excl": (a * sub["n_masked"]).sum() / (K * tot),
            "nonstim": a[sub["k_used"] < MIN_GRADIENTS].sum() / tot}


def fmt_fracs(f):
    return f"[escluse {100 * f['excl']:.2f}% | non stim. per maschere {100 * f['nonstim']:.2f}%]"


def main():
    ap = argparse.ArgumentParser(description="Analisi degli errori di un esperimento.")
    ap.add_argument("esperimento", help="nome della cartella in output/, es. aniso_const")
    out_dir = ROOT / "output" / ap.parse_args().esperimento
    csv_path = out_dir / "errors.csv"
    mask_csv_path = out_dir / "mask_pairs.csv"
    fig_dir = out_dir / "figures"

    df = pd.read_csv(csv_path, dtype={"experiment_id": str})
    ids = df["experiment_id"].unique()
    if len(ids) != 1:
        # face_id non e' univoco fra esperimenti: mescolarli darebbe statistiche senza senso
        raise SystemExit(f"{csv_path}: attesi dati di un solo esperimento, trovati {list(ids)}")
    exp_id = ids[0]
    fig_dir.mkdir(parents=True, exist_ok=True)
    # errori relativi con segno, (stimato - vero) / vero; NaN sulle facce non stimabili
    for cv in ("CV_L", "CV_T"):
        df[f"err_{cv}_rel"] = (df[f"{cv}_est"] - df[f"{cv}_true"]) / df[f"{cv}_true"]

    pairs = pd.read_csv(mask_csv_path, dtype={"experiment_id": str, "map_name": str})
    pids = pairs["experiment_id"].unique()
    if len(pids) != 1 or pids[0] != exp_id:
        raise SystemExit(f"{mask_csv_path}: attesi dati del solo esperimento {exp_id}, trovati {list(pids)}")
    map_names = list(dict.fromkeys(pairs["map_name"]))     # ordine di comparsa = ordine delle mappe
    K = len(map_names)
    per_face = pairs.groupby("face_id").size()
    if set(per_face.index) != set(df["face_id"]) or (per_face != K).any():
        raise SystemExit(f"{mask_csv_path}: ogni faccia di {csv_path} deve avere esattamente {K} righe")
    if not set(pairs["masked"].unique()) <= set(MASK_LABELS):
        raise SystemExit(f"{mask_csv_path}: masked fuori da {sorted(MASK_LABELS)}")
    pairs = pairs.merge(df[["face_id", "area"]], on="face_id", how="left")
    df["n_masked"] = df["face_id"].map(pairs.groupby("face_id")["masked"].sum())

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
    print(f"mappe              : {K} ({', '.join(map_names)})")
    print(f"k_used             : " + "  ".join(
        f"{k}={n}" for k, n in df["k_used"].value_counts().sort_index().items()))

    a = est["area"].values
    print(f"\nerrore angolare (gradi, pesato per area)  {fmt_fracs(mask_fracs(df, K))}")
    for q in (0.25, 0.50, 0.75, 0.90, 0.95, 0.99):
        print(f"   q{int(q * 100):02d} : {wquantile(est['err_angle'], a, q):8.3f}")
    frac20 = est.loc[est["err_angle"] < 20, "area"].sum() / est["area"].sum()
    print(f"   area entro 20 gradi : {100 * frac20:.2f}%")

    # ------------------------------------------------------- 1. per regione
    rule("1. ERRORE PER REGIONE ANATOMICA")
    rows = []
    for tag, g in est.groupby("region_tag"):
        w = g["area"].values
        n_tot = (df["region_tag"] == tag).sum()
        mf = mask_fracs(df[df["region_tag"] == tag], K)
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
            "escluse%": 100 * mf["excl"],
            "nonstim%": 100 * mf["nonstim"],
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
    plt.savefig(fig_dir / "err_per_regione.png", dpi=140)
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
              f"q25 {q25:+8.4f}   q50 {q50:+8.4f}   q75 {q75:+8.4f}   "
              f"{fmt_fracs(mask_fracs(df[df['region_tag'] == tag], K))}")

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
    plt.savefig(fig_dir / "bias_cv.png", dpi=140)
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
        print("\nUna faccia e' non stimabile se dopo la maschera di collisione le\n"
              "restano meno di 3 mappe, se sigma_min e' sotto soglia (fronti quasi\n"
              "paralleli) o se lambda_2 <= 0. Ci si aspetta che stiano sulle linee\n"
              "di collisione dei fronti. Se sono sparse a caso, c'e' altro da capire.")
        bad[["face_id", "region_tag", "area"]].to_csv(
            fig_dir / "facce_non_stimabili.csv", index=False)
        print(f"\nelenco scritto in {fig_dir / 'facce_non_stimabili.csv'}\n"
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
    print("\n  bin di cond            facce    ang_q50    ang_q90   CVL_q50   (facce stimabili del bin)")
    for b in range(1, len(bins)):
        g = est[idx == b]
        if len(g) < 20:
            continue
        w = g["area"].values
        print(f"  [{bins[b-1]:8.1f},{bins[b]:9.1f})  {len(g):7d}  "
              f"{wquantile(g['err_angle'], w, 0.50):9.3f}  "
              f"{wquantile(g['err_angle'], w, 0.90):9.3f}  "
              f"{wquantile(g['err_CV_L_rel'], w, 0.50):+9.4f}   "
              f"{fmt_fracs(mask_fracs(g, K))}")

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
    plt.savefig(fig_dir / "err_vs_cond.png", dpi=140)
    plt.close()

    # ------------------------------------------- 5. r su tenute vs escluse
    rule("5. r = g^T M_true g - 1: COPPIE TENUTE CONTRO ESCLUSE")
    print("r = 0: la coppia e' coerente con un fronte piano nel tessuto vero.\n"
          "Metrica di valutazione: non entra nella stima ne' nelle maschere.\n")
    print("coppie per mappa:")
    tab = pairs.pivot_table(index="map_name", columns="masked", values="r",
                            aggfunc="size", fill_value=0).reindex(index=map_names, columns=list(MASK_LABELS), fill_value=0)
    tab.columns = [MASK_LABELS[c] for c in tab.columns]
    tab.loc["totale"] = tab.sum()
    print(tab.to_string())
    print("\nfrazione di area delle coppie (area / (K * area totale)):")
    for c, name in MASK_LABELS.items():
        print(f"  {name:10s} {100 * pairs.loc[pairs['masked'] == c, 'area'].sum() / (K * area_tot):8.3f}%")

    pairs["abs_r"] = pairs["r"].abs()
    qs = (0.05, 0.25, 0.50, 0.75, 0.95)
    for col in ("r", "abs_r"):
        print(f"\n{'|r|' if col == 'abs_r' else 'r'} (pesato per area)")
        print("  " + f"{'':10s} {'coppie':>8s} " + " ".join(f"{'q' + format(int(q * 100), '02d'):>10s}" for q in qs))
        for c, name in MASK_LABELS.items():
            g = pairs[pairs["masked"] == c]
            vals = " ".join(f"{wquantile(g[col], g['area'], q):10.4f}" for q in qs) if len(g) else "  -"
            print(f"  {name:10s} {len(g):8d} {vals}")

    lo, hi = np.percentile(pairs["r"], [0.5, 99.5])
    bins = np.linspace(lo, hi, 121)
    fig, ax = plt.subplots(figsize=(9, 5))
    for c, name in MASK_LABELS.items():
        g = pairs[pairs["masked"] == c]
        if len(g) == 0:
            continue
        ax.hist(g["r"].clip(lo, hi), bins=bins, weights=g["area"], density=True,
                histtype="step", lw=2, color=MASK_COLORS[c], label=f"{name} ({len(g)})")
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("r = g^T M_true g - 1  (clip ai percentili 0.5 / 99.5)")
    ax.set_ylabel("densita' pesata per area")
    ax.set_title(f"{exp_id}: r sulle coppie (faccia, mappa) tenute ed escluse")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_dir / "r_mask.png", dpi=140)
    plt.close()

    rule("FIGURE")
    for f in sorted(fig_dir.glob("*.png")):
        print(f"  {f}")


if __name__ == "__main__":
    main()