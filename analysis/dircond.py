#!/usr/bin/env python3
"""
Cond direzionale per faccia: il numero di condizionamento del sistema dello
stimatore costruito dalle sole DIREZIONI dei gradienti LAT.

Riga per mappa: (c^2, sqrt(2) c s, s^2), con (c, s) = g/|g| in una base
ortonormale del piano della faccia (embedding di Frobenius: il risultato non
dipende dalla base). Le coppie (faccia, mappa) escluse dalle maschere si
prendono da mask_pairs.csv (masked = 1). Con normalize=False le righe sono
quelle dello stimatore C++ e cond riproduce la sua colonna cond.

Gradienti con face_gradient di cv_from_lat.py, LAT e mesh dal lat_ply del
.conf dell'esperimento.

Uso come modulo:
    from dircond import face_conds
    cond, cond_dir = face_conds("aniso_const")
"""

import numpy as np
import pandas as pd

from analisi import MIN_GRADIENTS, ROOT
from cv_from_lat import face_gradient


def conf_value(exp, key):
    for line in open(ROOT / "experiments" / f"{exp}.conf"):
        k, _, v = line.partition("=")
        if not line.startswith("#") and k.strip() == key:
            return v.strip()
    raise SystemExit(f"{exp}.conf: manca {key}")


def read_lat_ply(path):
    """Vertici, facce e mappe {nome: lat} dal PLY di multiLatToPly.py."""
    with open(path, "rb") as f:
        names, props = {}, []
        while (line := f.readline().decode("ascii").strip()) != "end_header":
            p = line.split()
            if p[:2] == ["element", "vertex"]:
                nv = int(p[2])
            elif p[:2] == ["element", "face"]:
                nf = int(p[2])
            elif p[0] == "property" and p[1] == "float":
                props.append(p[2])
            elif p[0] == "comment" and p[1].startswith("lat_"):
                names[p[1]] = p[2]
        v = np.fromfile(f, dtype=[(q, "<f4") for q in props], count=nv)
        fc = np.fromfile(f, dtype=[("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")], count=nf)
    P = np.stack([v["x"], v["y"], v["z"]], 1).astype(np.float64)
    E = np.stack([fc["a"], fc["b"], fc["c"]], 1).astype(np.int64)
    return P, E, {names[k]: v[k].astype(np.float64) for k in names}


def system_cond(P, E, lats, masked, normalize):
    """cond = smax/smin del sistema K x 3 per faccia; NaN con meno di MIN_GRADIENTS righe.
    masked: {nome mappa: bool per faccia}. Le righe escluse diventano zero: non cambiano i valori singolari."""
    p0, p1, p2 = P[E[:, 0]], P[E[:, 1]], P[E[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    n /= np.linalg.norm(n, axis=1)[:, None]
    u = (p1 - p0) / np.linalg.norm(p1 - p0, axis=1)[:, None]
    w = np.cross(n, u)
    rows = []
    for name, lat in lats.items():
        g = face_gradient(P, E, lat)[0]
        c, s = (g * u).sum(1), (g * w).sum(1)
        if normalize:                                              # gradiente nullo: nessuna direzione, riga zero
            norm = np.hypot(c, s)
            c = np.divide(c, norm, out=np.zeros_like(c), where=norm > 0)
            s = np.divide(s, norm, out=np.zeros_like(s), where=norm > 0)
        R = np.stack([c * c, np.sqrt(2.0) * c * s, s * s], 1)
        rows.append(np.where(masked[name][:, None], 0.0, R))
    sv = np.linalg.svd(np.stack(rows, 1), compute_uv=False)        # (F, 3)
    with np.errstate(divide="ignore", invalid="ignore"):           # rango < 3: cond infinito
        cond = sv[:, 0] / sv[:, -1]
    used = sum(~m for m in masked.values())
    return np.where(used >= MIN_GRADIENTS, cond, np.nan)


def face_conds(exp):
    """(cond con le righe dello stimatore, cond direzionale) per faccia dell'esperimento."""
    P, E, lats = read_lat_ply(ROOT / conf_value(exp, "lat_ply"))
    mp = pd.read_csv(ROOT / conf_value(exp, "output_mask_csv"))
    masked = {name: (g.sort_values("face_id")["masked"].to_numpy() == 1) for name, g in mp.groupby("map_name")}
    lats = {name: lats[name] for name in masked}                   # le sole mappe dell'esperimento
    return system_cond(P, E, lats, masked, False), system_cond(P, E, lats, masked, True)
