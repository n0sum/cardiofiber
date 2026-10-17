#!/usr/bin/env python3
"""
Diagnosi per-stadio del refinement (ICP+ZoomOut) sulla functional map
tra le due mesh cardiache, e test del termine orientation-preserving
nativo di geomfum per stabilizzare l'init contro il flip di simmetria.

Convenzione (identica a lat_transfer.py, verificata su main.cpp):
    mesh1elaborata.ply = 2-LA.mesh    = NORMALE = TARGET = b
    mesh2elaborata.ply = 2-LA-FA.mesh = FA      = SOURCE = a
    mesh_src = TriangleMesh(V2, F2)   (a)
    mesh_tgt = TriangleMesh(V1, F1)   (b)
    p2p = converter(fmap, basis_src, basis_tgt)  shape=[n_tgt], values in [0, n_src)

Nessun transfer LAT, nessun ancoraggio spaziale: solo qualità della fmap.
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from landmarks import get_matched_landmarks

try:
    from geomfum.shape import TriangleMesh
    from geomfum.descriptor.spectral import (
        HeatKernelSignature, WaveKernelSignature,
        LandmarkHeatKernelSignature, LandmarkWaveKernelSignature,
    )
    from geomfum.descriptor.pipeline import DescriptorPipeline, L2InnerNormalizer
    from geomfum.functional_map import (
        FunctionalMap, SDPFactorBuilder, LBCFactorBuilder, MultFactorBuilder,
        SpectralDescriptorPreservation, LBCommutativityEnforcing,
        OperatorCommutativityEnforcing, FactorSum, ScipyMinimize,
    )
    from geomfum.refine import IcpRefiner, ZoomOut
    from geomfum.convert import P2pFromFmConverter, SinkhornP2pFromFmConverter
except ImportError as exc:
    sys.exit(f"Import error: {exc}\nInstall geomfum senza torch (vedi memory).")


# ── PLY I/O ───────────────────────────────────────────────────────────────────

def read_ply(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    end = next(i for i, l in enumerate(lines) if l.strip() == 'end_header')
    header = [l.strip() for l in lines[:end]]
    body = lines[end + 1:]
    n_verts, n_faces = 0, 0
    vert_props = []
    in_vert, in_face = False, False
    for line in header:
        if line.startswith('element vertex'):
            n_verts = int(line.split()[-1]); in_vert, in_face = True, False
        elif line.startswith('element face'):
            n_faces = int(line.split()[-1]); in_vert, in_face = False, True
        elif line.startswith('property') and in_vert and 'list' not in line:
            vert_props.append(line.split()[-1])
    V = np.array([[float(v) for v in body[i].split()] for i in range(n_verts)])
    F = np.array(
        [list(map(int, body[n_verts + i].split()[1:4])) for i in range(n_faces)],
        dtype=np.int32,
    )
    return {name: V[:, i] for i, name in enumerate(vert_props)}, F


# ── Integrità mesh ────────────────────────────────────────────────────────────

def _count_boundary_loops(boundary_edges):
    if not boundary_edges:
        return 0
    adj = defaultdict(set)
    for a, b in boundary_edges:
        adj[a].add(b); adj[b].add(a)
    visited, n_loops = set(), 0
    for start in adj:
        if start in visited:
            continue
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            for nb in adj[v]:
                if nb not in visited:
                    stack.append(nb)
        n_loops += 1
    return n_loops


def mesh_integrity(V, F):
    n_v = len(V)
    row = np.concatenate([F[:, 0], F[:, 1], F[:, 2]])
    col = np.concatenate([F[:, 1], F[:, 2], F[:, 0]])
    adj = csr_matrix((np.ones(len(row)), (row, col)), shape=(n_v, n_v))
    n_comp, _ = connected_components(adj, directed=False)

    edge_count = defaultdict(int)
    for fac in F:
        for a, b in [(fac[0], fac[1]), (fac[1], fac[2]), (fac[2], fac[0])]:
            edge_count[(min(a, b), max(a, b))] += 1
    boundary_edges = [(a, b) for (a, b), c in edge_count.items() if c == 1]
    nonmanifold_edges = [(a, b) for (a, b), c in edge_count.items() if c > 2]
    n_loops = _count_boundary_loops(boundary_edges)
    return dict(n_components=n_comp, n_boundary_loops=n_loops,
                n_nonmanifold_edges=len(nonmanifold_edges))


# ── Metriche fmap / p2p ───────────────────────────────────────────────────────

def diagonality(fmap, band=3):
    total = np.sum(fmap ** 2)
    if total == 0:
        return 0.0
    rows, cols = np.indices(fmap.shape)
    mask = np.abs(rows - cols) <= band
    return float(np.sum(fmap[mask] ** 2) / total)


def trace_metric(fmap):
    k = min(fmap.shape)
    return float(np.trace(fmap[:k, :k]))


def p2p_from_fmap(fmap, basis_src, basis_tgt, use_sinkhorn=True):
    if use_sinkhorn:
        try:
            return np.asarray(SinkhornP2pFromFmConverter()(fmap, basis_src, basis_tgt), dtype=np.int64)
        except Exception:
            pass
    return np.asarray(P2pFromFmConverter()(fmap, basis_src, basis_tgt), dtype=np.int64)


def corr_metrics(V_src, V_tgt, p2p, basis_src, basis_tgt):
    """Correlazioni coordinate (x,y,z) e autovettori LBO 1-3, trasferiti src->tgt via p2p."""
    out = {}
    for i, name in enumerate(['x', 'y', 'z']):
        fn_xfer = V_src[p2p, i]
        fn_native = V_tgt[:, i]
        std_x, std_n = fn_xfer.std(), fn_native.std()
        out[f'corr_{name}'] = (float(np.corrcoef(fn_xfer, fn_native)[0, 1])
                                if std_x > 1e-10 and std_n > 1e-10 else float('nan'))
    for ev_idx in [1, 2, 3]:
        ev_src = np.asarray(basis_src.full_vecs[:, ev_idx], dtype=float)
        ev_tgt = np.asarray(basis_tgt.full_vecs[:, ev_idx], dtype=float)
        ev_xfer = ev_src[p2p]
        if np.dot(ev_xfer, ev_tgt) < 0:
            ev_xfer = -ev_xfer
        std_x, std_n = ev_xfer.std(), ev_tgt.std()
        out[f'corr_ev{ev_idx}'] = (float(np.corrcoef(ev_xfer, ev_tgt)[0, 1])
                                    if std_x > 1e-10 and std_n > 1e-10 else float('nan'))
    return out


def aggregate_score(corr):
    """Media delle correlazioni (autovettori LBO + coordinate), con penalità se
    una correlazione coordinata è negativa. Non dipende da tr(fmap): tra atri di
    forma diversa la mappa corretta non è vicina all'identità, quindi il segno
    della traccia non è un discriminante affidabile (vedi prompt/memoria)."""
    vals = [corr.get(f'corr_{k}') for k in ('ev1', 'ev2', 'ev3', 'x', 'y', 'z')]
    vals = [v for v in vals if v is not None and not np.isnan(v)]
    if not vals:
        return -1e9
    s = float(np.mean(vals))
    for k in ('x', 'y', 'z'):
        c = corr.get(f'corr_{k}')
        if c is not None and not np.isnan(c) and c < 0:
            s -= 0.5
    return s


def is_true_flip(corr, agg_threshold=0.15):
    """Rileva un flip geometrico vero: correlazione coordinata dominante negativa
    (riflessione sinistra/destra) oppure correlazione aggregata (coordinate +
    autovettori LBO) troppo bassa. Il segno di tr(fmap) NON è usato: è solo un
    proxy valido quando la mappa è vicina all'identità, il che non è il caso qui."""
    coord = [corr[f'corr_{k}'] for k in ('x', 'y', 'z')]
    coord = [c for c in coord if not np.isnan(c)]
    dominant_neg_coord = any(c < 0 for c in coord)
    vals = [corr[f'corr_{k}'] for k in ('ev1', 'ev2', 'ev3', 'x', 'y', 'z')]
    vals = [v for v in vals if not np.isnan(v)]
    low_agg = (np.mean(vals) < agg_threshold) if vals else True
    return bool(dominant_neg_coord or low_agg)


def evaluate(label, fmap, V_src, V_tgt, basis_src, basis_tgt, log):
    tr = trace_metric(fmap)
    diag = diagonality(fmap)
    p2p = p2p_from_fmap(fmap, basis_src, basis_tgt)
    corr = corr_metrics(V_src, V_tgt, p2p, basis_src, basis_tgt)
    sc = aggregate_score(corr)
    flip = is_true_flip(corr)
    flip_str = "*** FLIP (geom) ***" if flip else "ok"
    log(f"  [{label:22s}] shape={str(fmap.shape):12s} tr={tr:+7.3f}  diag={diag:.3f}  score={sc:+.3f}  "
        f"corr(x,y,z)=({corr['corr_x']:+.3f},{corr['corr_y']:+.3f},{corr['corr_z']:+.3f})  "
        f"corr(ev1,ev2,ev3)=({corr['corr_ev1']:+.3f},{corr['corr_ev2']:+.3f},{corr['corr_ev3']:+.3f})  [{flip_str}]")
    return dict(label=label, tr=tr, diag=diag, score=sc, flip=flip, **corr, fmap=fmap, p2p=p2p)


def compute_fmap_orient(mesh_src, mesh_tgt, descr_src, descr_tgt,
                         w_descr, w_lbc, w_mult, w_orient):
    proj_src = mesh_src.basis.project(descr_src)
    proj_tgt = mesh_tgt.basis.project(descr_tgt)
    k_src = mesh_src.basis.spectrum_size
    k_tgt = mesh_tgt.basis.spectrum_size

    factors = [
        SpectralDescriptorPreservation(proj_src, proj_tgt, weight=w_descr),
        LBCommutativityEnforcing.from_bases(mesh_src.basis, mesh_tgt.basis, weight=w_lbc),
        OperatorCommutativityEnforcing.from_multiplication(
            mesh_src.basis, descr_src, mesh_tgt.basis, descr_tgt, weight=w_mult
        ),
        OperatorCommutativityEnforcing.from_orientation(
            mesh_src, descr_src, mesh_tgt, descr_tgt, weight=w_orient
        ),
    ]
    objective = FactorSum(factors)
    x0 = np.zeros((k_tgt, k_src))
    res = ScipyMinimize(method="L-BFGS-B").minimize(objective, x0, fun_jac=objective.gradient)
    return res.x.reshape(k_tgt, k_src)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mesh1', help='target/normale/b (mesh1elaborata.ply)')
    ap.add_argument('mesh2', help='source/FA/a (mesh2elaborata.ply)')
    ap.add_argument('--outdir', default='diag_out')
    ap.add_argument('--k', type=int, default=120)
    ap.add_argument('--use-k', type=int, default=30)
    ap.add_argument('--n-times', type=int, default=100)
    ap.add_argument('--w-descr', type=float, default=1.0)
    ap.add_argument('--w-lbc', type=float, default=1e-2)
    ap.add_argument('--w-mult', type=float, default=1e-1)
    ap.add_argument('--icp-nit', type=int, default=10)
    ap.add_argument('--zo-nit', type=int, default=6)
    ap.add_argument('--zo-step', type=int, default=5)
    ap.add_argument('--orient-ws', type=str, default='0.5,1,5')
    ap.add_argument('--n-times-landmarks', type=int, default=3,
                     help='HKS/WKS time steps per landmark (default 3)')
    ap.add_argument('--repick', action='store_true',
                     help='Forza il re-picking dei landmark anche se la cache esiste')
    ap.add_argument('--no-landmarks', action='store_true',
                     help='Disabilita i descrittori landmark (baseline per confronto)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    report = []

    def log(line):
        print(line)
        report.append(line)

    log(f"Leggo mesh1 (target/normale/b): {args.mesh1}")
    props1, F1 = read_ply(args.mesh1)
    V1 = np.column_stack([props1['x'], props1['y'], props1['z']])
    log(f"Leggo mesh2 (source/FA/a): {args.mesh2}")
    props2, F2 = read_ply(args.mesh2)
    V2 = np.column_stack([props2['x'], props2['y'], props2['z']])

    log("\n" + "=" * 70)
    log("PASSO 0 — Integrità (assunta post fix_holes.py)")
    log("=" * 70)
    rep1 = mesh_integrity(V1, F1)
    rep2 = mesh_integrity(V2, F2)
    for name, r in [('mesh1/target/normale', rep1), ('mesh2/source/FA', rep2)]:
        log(f"  {name}: comp={r['n_components']} loop={r['n_boundary_loops']} nonmanifold={r['n_nonmanifold_edges']}")
        assert r['n_components'] == 1, f"{name}: componenti connesse = {r['n_components']} (atteso 1)"
        assert r['n_boundary_loops'] == 0, f"{name}: loop di bordo = {r['n_boundary_loops']} (atteso 0)"
        assert r['n_nonmanifold_edges'] == 0, f"{name}: spigoli non-manifold = {r['n_nonmanifold_edges']} (atteso 0)"
    log("  Assert integrità: OK")

    if not args.no_landmarks:
        log("\n" + "=" * 70)
        log("PASSO 0.5 — Landmark anatomici (pre-fmap)")
        log("=" * 70)
        lm1, lm2 = get_matched_landmarks(
            args.mesh1, V1, F1, args.mesh2, V2, F2, repick=args.repick,
        )
        log(f"  landmark: {len(lm1)} punti per mesh, ordine corrispondente garantito")

    mesh_src = TriangleMesh(V2, F2)   # a = source = FA
    mesh_tgt = TriangleMesh(V1, F1)   # b = target = normale
    if not args.no_landmarks:
        mesh_src.landmark_indices = np.asarray(lm2, dtype=np.int64)
        mesh_tgt.landmark_indices = np.asarray(lm1, dtype=np.int64)
    mesh_src.vertices -= mesh_src.vertices.mean(axis=0)
    mesh_tgt.vertices -= mesh_tgt.vertices.mean(axis=0)
    V_src = np.array(mesh_src.vertices)
    V_tgt = np.array(mesh_tgt.vertices)

    K = args.k
    log(f"\nLBO spectrum K={K} (src n={len(V_src)}, tgt n={len(V_tgt)}) ...")
    mesh_src.laplacian.find_spectrum(spectrum_size=K, set_as_basis=True)
    mesh_tgt.laplacian.find_spectrum(spectrum_size=K, set_as_basis=True)

    log(f"Descrittori HKS+WKS (n_times={args.n_times})"
        + ("" if args.no_landmarks else
           f" + LandmarkHKS+LandmarkWKS (n_times_landmarks={args.n_times_landmarks})")
        + " ...")
    hks = HeatKernelSignature.from_registry(n_domain=args.n_times)
    wks = WaveKernelSignature.from_registry(n_domain=args.n_times)
    steps = [hks, wks]
    if not args.no_landmarks:
        lm_hks = LandmarkHeatKernelSignature.from_registry(n_domain=args.n_times_landmarks)
        lm_wks = LandmarkWaveKernelSignature.from_registry(n_domain=args.n_times_landmarks)
        steps += [lm_hks, lm_wks]
    pipeline = DescriptorPipeline(steps + [L2InnerNormalizer()])
    descr_src = pipeline.apply(mesh_src)
    descr_tgt = pipeline.apply(mesh_tgt)

    use_k = args.use_k
    mesh_src.basis.use_k = use_k
    mesh_tgt.basis.use_k = use_k

    log(f"\n{'=' * 70}\nPASSO 1 — Valutazione per-stadio (pipeline attuale, use_k={use_k})\n{'=' * 70}")
    solver = FunctionalMap(factor_builders=[
        SDPFactorBuilder(weight=args.w_descr),
        LBCFactorBuilder(weight=args.w_lbc),
        MultFactorBuilder(weight=args.w_mult),
    ])
    fmap_a = solver(mesh_src.basis, mesh_tgt.basis, descr_src, descr_tgt)
    res_a = evaluate('(a) init grezzo', fmap_a, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    fmap_b = IcpRefiner(nit=args.icp_nit)(fmap_a.copy(), mesh_src.basis, mesh_tgt.basis)
    res_b = evaluate('(b) solo ICP', fmap_b, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    fmap_c = ZoomOut(nit=args.zo_nit, step=args.zo_step)(fmap_a.copy(), mesh_src.basis, mesh_tgt.basis)
    res_c = evaluate('(c) solo ZoomOut', fmap_c, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    fmap_d = IcpRefiner(nit=args.icp_nit)(fmap_a.copy(), mesh_src.basis, mesh_tgt.basis)
    fmap_d = ZoomOut(nit=args.zo_nit, step=args.zo_step)(fmap_d, mesh_src.basis, mesh_tgt.basis)
    res_d = evaluate('(d) ICP+ZoomOut', fmap_d, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    fmap_g = ZoomOut(nit=30, step=1)(fmap_a.copy(), mesh_src.basis, mesh_tgt.basis)
    res_g = evaluate('(g) landmark+ZO(step=1)', fmap_g, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    log(f"\n{'=' * 70}\nPASSO 2 — Init rafforzato con orientation-preserving\n{'=' * 70}")
    orient_results = []
    for w in [float(x) for x in args.orient_ws.split(',')]:
        fmap_o = compute_fmap_orient(mesh_src, mesh_tgt, descr_src, descr_tgt,
                                      args.w_descr, args.w_lbc, args.w_mult, w)
        res_o = evaluate(f'(a\') W={w:g}', fmap_o, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)
        orient_results.append((w, res_o))

    non_flip_orient = [t for t in orient_results if not t[1]['flip']]
    pool = non_flip_orient if non_flip_orient else orient_results
    best_w, best_res_o = max(pool, key=lambda t: t[1]['score'])
    log(f"\n  Migliore: W={best_w:g}  (score={best_res_o['score']:+.3f}, tr={best_res_o['tr']:+.3f}, "
        f"flip={best_res_o['flip']})")

    log(f"\n{'=' * 70}\nPASSO 3 — Refinement stabile dall'init rafforzato (W={best_w:g})\n{'=' * 70}")
    fmap_e = ZoomOut(nit=30, step=1)(best_res_o['fmap'].copy(), mesh_src.basis, mesh_tgt.basis)
    res_e = evaluate('(e) reinforced+ZO(step=1)', fmap_e, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    fmap_f = IcpRefiner(nit=args.icp_nit)(best_res_o['fmap'].copy(), mesh_src.basis, mesh_tgt.basis)
    fmap_f = ZoomOut(nit=args.zo_nit, step=args.zo_step)(fmap_f, mesh_src.basis, mesh_tgt.basis)
    res_f = evaluate('(f) reinforced+ICP+ZO', fmap_f, V_src, V_tgt, mesh_src.basis, mesh_tgt.basis, log)

    all_results = [res_a, res_b, res_c, res_d, res_g, best_res_o, res_e, res_f]

    def _complexity(res):
        label = res['label']
        if label.startswith("(a')"):
            return 1  # init rafforzato con orientation-preserving
        return {
            '(a) init grezzo': 0,
            '(b) solo ICP': 1,
            '(c) solo ZoomOut': 1,
            '(g) landmark+ZO(step=1)': 1,
            '(d) ICP+ZoomOut': 2,
            '(e) reinforced+ZO(step=1)': 2,
            '(f) reinforced+ICP+ZO': 3,
        }.get(label, 99)

    # Selezione: massimizza la correlazione aggregata (score), scartando solo le
    # varianti con flip geometrico vero (corr_x<0 o correlazione aggregata bassa).
    # Il segno di tr(fmap) NON è usato come discriminante (vedi aggregate_score/is_true_flip).
    valid = [r for r in all_results if not r['flip']]
    if valid:
        best_score = max(r['score'] for r in valid)
        near_best = [r for r in valid if best_score - r['score'] <= 0.02]
        best_final = min(near_best, key=_complexity)
    else:
        best_final = max(all_results, key=lambda r: r['score'])

    log(f"\n{'=' * 70}\nPASSO 4 — Tabella accettazione (prima/dopo)\n{'=' * 70}")
    header = (f"  {'stadio':28s} {'score':>7s} {'tr':>8s} {'diag':>6s} {'corr_x':>7s} {'corr_y':>7s} "
              f"{'corr_z':>7s} {'corr_ev1':>8s} {'corr_ev2':>8s} {'corr_ev3':>8s} {'flip':>6s}")
    log(header)
    for r in all_results:
        log(f"  {r['label']:28s} {r['score']:+7.3f} {r['tr']:+8.3f} {r['diag']:6.3f} {r['corr_x']:+7.3f} {r['corr_y']:+7.3f} "
            f"{r['corr_z']:+7.3f} {r['corr_ev1']:+8.3f} {r['corr_ev2']:+8.3f} {r['corr_ev3']:+8.3f} {str(r['flip']):>6s}")
    log(f"\n  Scelta finale: {best_final['label']}  (score={best_final['score']:+.3f}, "
        f"flip geometrico: {best_final['flip']}, tr={best_final['tr']:+.3f})")
    log(f"  Verifica anti-mascheramento: corr_x della scelta finale = {best_final['corr_x']:+.3f} "
        f"({'NEGATIVA -> possibile flip mascherato!' if best_final['corr_x'] < 0 else 'positiva, nessun flip mascherato'})")

    report_path = os.path.join(args.outdir, 'report_stages.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    log(f"\nReferto: {report_path}")


if __name__ == '__main__':
    main()
