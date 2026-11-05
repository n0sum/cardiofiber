#!/usr/bin/env python3
"""
Metodo 2 di Roney et al. 2019 su terne di mappe LAT, per il confronto con lo
stimatore per faccia (estimator/src): patch geodetica, fit planare per
mappa, ellisse centrata sui vettori velocita' g / |g|^2.

Con tre punti la conica passa esattamente per tutti e tre: e' anche la
soluzione del fit a residui radiali di Roney, quindi l'algoritmo genetico non
serve.

Uso, dalla radice del repository:
    python estimator/roney.py experiments/roney.conf aniso_zone_1 aniso_zone_50 ...
Legge data/reference/LA_LAT_SCENARIO.ply e, solo dopo la stima,
data/reference/LA_GT_SCENARIO.ply. Le terne sono quelle delle mappe lat_k
presenti nel LAT PLY, cioe' che convergono; k e' l'indice del sito. Scrive
output/roney/SCENARIO/roney/sub_i-j-k/errors.csv, con le colonne di
exportErrorsAsCsv in Mesh.cpp.

Unita': gradienti in ms/mm, velocita' in mm/ms.
"""

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]                # radice del repository
KEYS = ["patch.radius_mm", "patch.n_points", "patch.seed", "fit.max_residual", "mask.min_singular_value"]
COLUMNS = ["experiment_id", "face_id", "area", "region_tag", "estimable", "k_used", "cond",
           "CV_L_est", "CV_T_est", "CV_L_true", "CV_T_true", "err_angle"]


def read_conf(path):
    kv = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if line:
            k, v = line.split("=")
            kv[k.strip()] = v.strip()
    return {k: float(kv[k]) for k in KEYS}               # KeyError se manca una chiave


def read_ply(path):
    """PLY binario di simulator/convert/toPly.py: vertici tutti float, facce
    'list uchar int' seguita da proprieta' float/int."""
    f = open(path, "rb")
    elems, cur = {}, None
    while (line := f.readline().decode("ascii").strip()) != "end_header":
        w = line.split()
        if w[0] == "element":
            cur = w[1]
            elems[cur] = [int(w[2]), []]
        elif w[0] == "property" and w[1] == "list":
            elems[cur][1] += [("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")]
        elif w[0] == "property":
            elems[cur][1].append((w[2], {"float": "<f4", "int": "<i4"}[w[1]]))
    out = {}
    for name, (n, dt) in elems.items():
        out[name] = np.fromfile(f, dtype=np.dtype(dt), count=n)
    return out


def geometry(V, F, radius, n_points, seed):
    """Per faccia: indici dei punti scelti (riempiti con -1 fino a n_points), coordinate 2D
    (X, Y) dei punti, mappa A (3x2) dal piano MDS a R^3, normale della faccia."""
    nv, nf = len(V), len(F)
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    w = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
    G = coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(nv, nv)).tocsr()

    # vicini entro il raggio di ogni vertice, a blocchi per non tenere la matrice densa
    near = [None] * nv
    for s in range(0, nv, 2000):
        idx = np.arange(s, min(s + 2000, nv))
        D = dijkstra(G, indices=idx, limit=radius)
        for r, v in enumerate(idx):
            j = np.flatnonzero(np.isfinite(D[r]))
            near[v] = (j, D[r, j])

    rng = np.random.default_rng(seed)
    pts = -np.ones((nf, n_points), dtype=np.int64)
    XY = np.zeros((nf, n_points, 2))
    A = np.zeros((nf, 3, 2))
    c_all = V[F].mean(axis=1)
    for i in range(nf):
        c = c_all[i]
        dist = np.full(nv, np.inf)
        for v in F[i]:
            j, d = near[v]
            dist[j] = np.minimum(dist[j], np.linalg.norm(c - V[v]) + d)
        patch = np.flatnonzero(dist <= radius)

        Dp = dijkstra(G[patch][:, patch])
        n = len(patch)
        J = np.eye(n) - 1.0 / n
        lam, U = np.linalg.eigh(-0.5 * J @ (Dp ** 2) @ J)
        xy = U[:, [-1, -2]] * np.sqrt(np.maximum(lam[[-1, -2]], 0))   # MDS classica, centrata

        P3 = V[patch] - V[patch].mean(axis=0)
        A[i] = np.linalg.lstsq(xy, P3, rcond=None)[0].T

        sel = rng.choice(n, n_points, replace=False) if n > n_points else np.arange(n)
        pts[i, :len(sel)] = patch[sel]
        XY[i, :len(sel)] = xy[sel]

    nrm = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    area = np.linalg.norm(nrm, axis=1) / 2
    return pts, XY, A, nrm / np.linalg.norm(nrm, axis=1)[:, None], area


def estimate(X, ok, A, nrm, min_sv):
    """X: (nf, 3, 2) vettori velocita', ok: (nf,) tutte e tre le mappe valide.
    Conica x^T P x = 1 per i tre punti; un'ellisse di semiassi CV_L, CV_T ha conica M^-1, M = P^-1."""
    R = np.stack([X[..., 0] ** 2, np.sqrt(2) * X[..., 0] * X[..., 1], X[..., 1] ** 2], axis=-1)
    R[~ok] = np.eye(3)                                    # righe fittizie: la faccia resta non stimabile
    ok = ok & (np.linalg.svd(R, compute_uv=False)[:, -1] > min_sv)
    R[~ok] = np.eye(3)
    p = np.linalg.solve(R, np.ones(R.shape[:2])[..., None])[..., 0]
    P = np.stack([np.stack([p[:, 0], p[:, 1] / np.sqrt(2)], -1),
                  np.stack([p[:, 1] / np.sqrt(2), p[:, 2]], -1)], -2)
    lp, U = np.linalg.eigh(P)                             # crescenti; M ha gli stessi autovettori
    ok = ok & (lp[:, 0] > 0)
    lp = np.where(ok[:, None], lp, 1.0)
    lam, lam2, fib2 = 1 / lp[:, 0], 1 / lp[:, 1], U[:, :, 0]   # lambda1 di M = 1 / autovalore minore di P
    f = np.einsum("fij,fj->fi", A, fib2)
    f -= np.sum(f * nrm, axis=1)[:, None] * nrm
    with np.errstate(divide="ignore", invalid="ignore"):
        f /= np.linalg.norm(f, axis=1)[:, None]
    return ok, np.sqrt(lam), np.sqrt(lam2), f


def main():
    conf = read_conf(sys.argv[1])
    radius, n_points = conf["patch.radius_mm"], int(conf["patch.n_points"])

    geo = None
    for scenario in sys.argv[2:]:
        lat = read_ply(ROOT / "data" / "reference" / f"LA_LAT_{scenario}.ply")
        F = np.stack([lat["face"][k] for k in "abc"], axis=1).astype(np.int64)
        if geo is None:                                   # la mesh e' la stessa in tutti gli scenari
            V = np.stack([lat["vertex"][k] for k in "xyz"], axis=1).astype(np.float64)
            geo = geometry(V, F, radius, n_points, int(conf["patch.seed"]))
            pts, XY, A, nrm, area = geo
            real = pts >= 0
            D = np.concatenate([real[..., None] * 1.0, XY], axis=-1)   # righe [1, X, Y], nulle se riempite
            pinv = np.linalg.pinv(D)                      # (nf, 3, n_points)
            n_real = real.sum(axis=1)

        maps = [int(p.removeprefix("lat_")) for p in lat["vertex"].dtype.names if p.startswith("lat_")]
        g, valid = [], []
        for k in maps:
            t = np.where(real, lat["vertex"][f"lat_{k}"].astype(np.float64)[np.maximum(pts, 0)], 0.0)
            coef = np.einsum("fij,fj->fi", pinv, t)
            res = np.einsum("fjk,fk->fj", D, coef) - t   # nulli sulle righe riempite
            rms = np.sqrt((res ** 2).sum(axis=1) / n_real)
            gk = coef[:, 1:]
            gn = np.linalg.norm(gk, axis=1)
            valid.append((gn > 0) & (rms <= conf["fit.max_residual"] * gn * radius))
            g.append(gk)
        g, valid = np.stack(g, axis=1), np.stack(valid, axis=1)        # (nf, mappe, 2), (nf, mappe)

        results = {}                                      # terna di siti -> stima
        for pos in itertools.combinations(range(len(maps)), 3):         # posizioni in maps
            s = list(pos)
            sub = tuple(maps[i] for i in pos)
            k_used, ok = valid[:, s].sum(axis=1), valid[:, s].all(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                x = g[:, s] / (g[:, s] ** 2).sum(axis=-1, keepdims=True)
            x = np.where(ok[:, None, None], x, 0.0)
            results[sub] = (k_used, *estimate(x, ok, A, nrm, conf["mask.min_singular_value"]))

        gt = read_ply(ROOT / "data" / "reference" / f"LA_GT_{scenario}.ply")["face"]
        fgt = np.stack([gt[k] for k in ("fx", "fy", "fz")], axis=1).astype(np.float64)
        for sub, (k_used, ok, cvl, cvt, f) in results.items():
            exp_id = "roney_sub_" + "-".join(map(str, sub))
            ang = np.degrees(np.arccos(np.minimum(1, np.abs(np.sum(f * fgt, axis=1)))))
            df = pd.DataFrame({"experiment_id": exp_id, "face_id": np.arange(len(F)), "area": area,
                               "region_tag": gt["elem_tag"], "estimable": ok.astype(int), "k_used": k_used,
                               "cond": np.nan, "CV_L_est": np.where(ok, cvl, np.nan),
                               "CV_T_est": np.where(ok, cvt, np.nan), "CV_L_true": gt["CV_L"],
                               "CV_T_true": gt["CV_T"], "err_angle": np.where(ok, ang, np.nan)})[COLUMNS]
            out = ROOT / "output" / "roney" / scenario / "roney" / ("sub_" + "-".join(map(str, sub)))
            out.mkdir(parents=True, exist_ok=True)
            df.to_csv(out / "errors.csv", index=False, float_format="%.9e")


if __name__ == "__main__":
    main()
