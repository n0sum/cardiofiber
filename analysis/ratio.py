#!/usr/bin/env python3
"""
Scansione in r = CV_L / CV_T a CV uniforme, contro aniso_const (r = 2).

Per ogni scenario, sulle facce stimabili sia nello scenario sia in
aniso_const, pesato per area: q50 e q90 di err_angle, rapporto con
aniso_const sulle stesse facce, previsioni f(r)/f(2) con
f = (r^2 + 1) / (r^2 - 1) e g(r)/g(2) con g = r^2 / (r^2 - 1), cond e cond
direzionale mediani (dircond.py), area non stimabile dello scenario.
cdir/ref: cond direzionale mediano dello scenario / quello di aniso_const
sulle stesse facce; lontano da 1, la copertura direzionale e' cambiata.

aniso_const_half (CV dimezzate, r = 2) e' il controllo di scala: deve dare
rapporti 1.

Tutta la mesh e, a parte, il solo corpo (tag 1).

Uso, dalla radice del repository:
    python analysis/ratio.py
"""

import pandas as pd

from analisi import ROOT, rule, wquantile
from dircond import conf_value, face_conds

REF = "aniso_const"
SCENARIOS = ["aniso_r110", "aniso_r125", "aniso_r150", "aniso_const",
             "aniso_r250", "aniso_r300", "aniso_const_half"]

f = lambda r: (r ** 2 + 1) / (r ** 2 - 1)
g = lambda r: r ** 2 / (r ** 2 - 1)


def load(exp):
    df = pd.read_csv(ROOT / conf_value(exp, "output_csv"))
    df["cond_dir"] = face_conds(exp)[1]
    cv = pd.read_csv(ROOT / conf_value(exp, "region_cv_csv"), comment="#")
    return df, cv["CV_L"].iloc[0] / cv["CV_T"].iloc[0]


def main():
    ref, r_ref = load(REF)
    data = {s: load(s) for s in SCENARIOS}
    for title, sel in (("TUTTA LA MESH", lambda d: d.index >= 0), ("CORPO (tag 1)", lambda d: d["region_tag"] == 1)):
        rule(f"SCANSIONE IN r CONTRO {REF}  -  {title}, facce stimabili in entrambi")
        rows = []
        for s in SCENARIOS:
            df, r = data[s]
            sub = df[sel(df)]
            both = sel(df) & (df["estimable"] == 1) & (ref["estimable"] == 1)
            a, b = df[both], ref[both]
            q = lambda x, p: wquantile(x[p[0]].to_numpy(), x["area"].to_numpy(), p[1])
            rows.append({"scenario": s, "r": r,
                         "q50": q(a, ("err_angle", .5)), "q90": q(a, ("err_angle", .9)),
                         "oss_q50": q(a, ("err_angle", .5)) / q(b, ("err_angle", .5)),
                         "oss_q90": q(a, ("err_angle", .9)) / q(b, ("err_angle", .9)),
                         "f/f(2)": f(r) / f(r_ref), "g/g(2)": g(r) / g(r_ref),
                         "cond": q(a, ("cond", .5)), "cond_dir": q(a, ("cond_dir", .5)),
                         "cdir/ref": q(a, ("cond_dir", .5)) / q(b, ("cond_dir", .5)),
                         "nonstim_%": 100 * sub.loc[sub["estimable"] == 0, "area"].sum() / sub["area"].sum()})
        print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\nangoli in gradi; oss = scenario / aniso_const sulle stesse facce; mediane pesate per area.")


if __name__ == "__main__":
    main()
