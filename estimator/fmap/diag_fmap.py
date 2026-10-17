#!/usr/bin/env python3
"""
Diagnosi qualità del functional map tra due mesh cardiache.

NON modifica la pipeline. Intercetta fmap iniziale, fmap finale, p2p e
produce: integrità mesh, PNG della fmap, PLY di funzioni geometriche
trasferite, errore di round-trip, referto testuale.

Usage:
    python diag_fmap.py mesh1.ply mesh2.ply [--outdir diag_out] [--k 120] ...

Convenzione: mesh1 = target/FA (b), mesh2 = source/normale (a).
             p2p shape=[n_b], valori in [0, n_a).
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from geomfum.shape import TriangleMesh
    from geomfum.descriptor.spectral import HeatKernelSignature, WaveKernelSignature
    from geomfum.descriptor.pipeline import DescriptorPipeline, L2InnerNormalizer
    from geomfum.functional_map import (
        FunctionalMap, SDPFactorBuilder, LBCFactorBuilder, MultFactorBuilder,
    )
    from geomfum.refine import IcpRefiner, ZoomOut
    from geomfum.convert import P2pFromFmConverter, SinkhornP2pFromFmConverter
except ImportError as exc:
    sys.exit(f"Import error: {exc}\nInstall geomfum senza torch (vedi memory).")


# ── PLY I/O (identico a lat_transfer.py) ─────────────────────────────────────

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


def write_ply(path, props, faces):
    n_verts = len(next(iter(props.values())))
    with open(path, 'w') as out:
        out.write("ply\nformat ascii 1.0\ncomment Written by diag_fmap.py\n")
        out.write(f"element vertex {n_verts}\n")
        for name in props:
            out.write(f"property float {name}\n")
        out.write(f"element face {len(faces)}\n")
        out.write("property list uchar int vertex_indices\nend_header\n")
        data = np.column_stack(list(props.values()))
        for row in data:
            out.write(' '.join(f'{v:.7g}' for v in row) + '\n')
        for fac in faces:
            out.write(f"3 {fac[0]} {fac[1]} {fac[2]}\n")


def write_q(path, V, F, quality):
    write_ply(path, {'x': V[:, 0], 'y': V[:, 1], 'z': V[:, 2],
                     'quality': np.asarray(quality, dtype=np.float64)}, F)


# ── Test 1: Integrità mesh ────────────────────────────────────────────────────

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
    n_v, n_f = len(V), len(F)

    # Componenti connesse
    row = np.concatenate([F[:, 0], F[:, 1], F[:, 2]])
    col = np.concatenate([F[:, 1], F[:, 2], F[:, 0]])
    adj = csr_matrix((np.ones(len(row)), (row, col)), shape=(n_v, n_v))
    n_comp, comp_labels = connected_components(adj, directed=False)
    comp_sizes = sorted(np.bincount(comp_labels).tolist(), reverse=True) if n_comp > 1 else []

    # Spigoli
    edge_count = defaultdict(int)
    for fac in F:
        for a, b in [(fac[0], fac[1]), (fac[1], fac[2]), (fac[2], fac[0])]:
            edge_count[(min(a, b), max(a, b))] += 1
    boundary_edges = [(a, b) for (a, b), c in edge_count.items() if c == 1]
    nonmanifold_edges = [(a, b) for (a, b), c in edge_count.items() if c > 2]
    n_loops = _count_boundary_loops(boundary_edges)

    # Facce degeneri
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    n_degen = int(np.sum(areas < 1e-10))

    # Vertici duplicati
    unique_v = np.unique(np.round(V, 6), axis=0)
    n_dupl = n_v - len(unique_v)

    return dict(
        n_verts=n_v, n_faces=n_f, n_components=n_comp, comp_sizes=comp_sizes,
        n_boundary_edges=len(boundary_edges), n_nonmanifold_edges=len(nonmanifold_edges),
        n_boundary_loops=n_loops, n_degenerate_faces=n_degen,
        min_face_area=float(areas.min()), n_duplicate_verts=n_dupl,
    )


def format_integrity(label, r):
    lines = [
        f"\n── Integrità: {label} ──",
        f"  Vertici         : {r['n_verts']}",
        f"  Facce           : {r['n_faces']}",
        f"  Componenti conn.: {r['n_components']}"
            + ("  *** PROBLEMA — spettro LBO si rompe su componenti multiple ***"
               if r['n_components'] > 1 else "  [ok]"),
    ]
    if r['n_components'] > 1:
        lines.append(f"  Dimensioni comp : {r['comp_sizes']} vertici")
    lines += [
        f"  Loop di bordo   : {r['n_boundary_loops']}  ({r['n_boundary_edges']} spigoli liberi)",
        f"  Non-manifold    : {r['n_nonmanifold_edges']}"
            + ("  [PROBLEMA]" if r['n_nonmanifold_edges'] > 0 else "  [ok]"),
        f"  Facce degeneri  : {r['n_degenerate_faces']}  (area min {r['min_face_area']:.2e})"
            + ("  [PROBLEMA]" if r['n_degenerate_faces'] > 0 else "  [ok]"),
        f"  Vertici duplicati: {r['n_duplicate_verts']}"
            + ("  [PROBLEMA]" if r['n_duplicate_verts'] > 0 else "  [ok]"),
    ]
    return lines


# ── Test 2: fmap PNG + diagonalità ────────────────────────────────────────────

def diagonality(fmap, band=3):
    """Frazione di ||fmap||^2_F contenuta entro |i-j|<=band dalla diagonale."""
    total = np.sum(fmap ** 2)
    if total == 0:
        return 0.0
    rows, cols = np.indices(fmap.shape)
    mask = np.abs(rows - cols) <= band
    return float(np.sum(fmap[mask] ** 2) / total)


def save_fmap_png(fmap, path, title):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    vmax = max(np.abs(fmap).max(), 1e-6)
    im = axes[0].imshow(fmap, cmap='bwr', vmin=-vmax, vmax=vmax,
                        origin='upper', aspect='auto')
    axes[0].set_title(title)
    axes[0].set_xlabel('src eigenindex'); axes[0].set_ylabel('tgt eigenindex')
    plt.colorbar(im, ax=axes[0])

    # Log-scale energy per diagonale
    k = min(fmap.shape)
    diag_energy = [np.sum(np.diagonal(fmap, d) ** 2)
                   for d in range(-(k - 1), k)]
    ds = np.arange(-(k - 1), k)
    axes[1].semilogy(ds, np.array(diag_energy) + 1e-30, linewidth=1)
    axes[1].axvline(0, color='r', linewidth=0.8, linestyle='--', label='diag')
    axes[1].set_title('Energia per diagonale'); axes[1].set_xlabel('offset dalla diag')
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ── Test 3: funzioni geometriche trasferite ───────────────────────────────────

def geometric_transfer_tests(V_src, F_src, V_tgt, F_tgt, p2p, basis_src, basis_tgt, outdir):
    lines = ["\n── Test 3: trasferimento funzioni geometriche ──"]
    n_src, n_tgt = len(V_src), len(V_tgt)

    # 3a) Coordinate (centrate)
    for i, name in enumerate(['x', 'y', 'z']):
        fn_src = V_src[:, i]
        fn_xfer = fn_src[p2p]
        fn_native = V_tgt[:, i]

        write_q(os.path.join(outdir, f'src_{name}.ply'), V_src, F_src, fn_src)
        write_q(os.path.join(outdir, f'tgt_{name}_xfer.ply'), V_tgt, F_tgt, fn_xfer)

        # Accordo di segno (emisfero): se la mappa è buona, sign(xfer) ~ sign(native)
        sign_match = 100 * np.mean(np.sign(fn_xfer) == np.sign(fn_native))
        # Correlazione lineare
        std_x = fn_xfer.std(); std_n = fn_native.std()
        corr = (float(np.corrcoef(fn_xfer, fn_native)[0, 1])
                if std_x > 1e-10 and std_n > 1e-10 else float('nan'))
        verdict = "ok" if sign_match > 75 else "WARN — possibile flip o mappa sbagliata"
        lines.append(f"  coord {name}: accordo segno {sign_match:.1f}%,  corr {corr:.3f}  [{verdict}]")

    # 3b) Distanza euclidea dal vertice più vicino al centroide (ref=centroid-nearest)
    ref_src = int(np.argmin(np.linalg.norm(V_src, axis=1)))
    dist_src = np.linalg.norm(V_src - V_src[ref_src], axis=1)
    dist_xfer = dist_src[p2p]
    ref_tgt = int(np.argmin(np.linalg.norm(V_tgt, axis=1)))
    dist_native = np.linalg.norm(V_tgt - V_tgt[ref_tgt], axis=1)

    write_q(os.path.join(outdir, 'src_dist_ref.ply'), V_src, F_src, dist_src)
    write_q(os.path.join(outdir, 'tgt_dist_ref_xfer.ply'), V_tgt, F_tgt, dist_xfer)
    write_q(os.path.join(outdir, 'tgt_dist_ref_native.ply'), V_tgt, F_tgt, dist_native)
    lines.append("  dist da ref: PLY salvati (src_dist_ref, tgt_dist_ref_xfer, tgt_dist_ref_native)")

    # 3c) Autovettori LBO 1, 2, 3 (funzioni globali più informative)
    for ev_idx in [1, 2, 3]:
        ev_src = np.asarray(basis_src.full_vecs[:, ev_idx], dtype=float)
        ev_native = np.asarray(basis_tgt.full_vecs[:, ev_idx], dtype=float)
        ev_xfer = ev_src[p2p]

        # Risolvi ambiguità di segno: allinea il trasferito al nativo
        if np.dot(ev_xfer, ev_native) < 0:
            ev_xfer = -ev_xfer

        write_q(os.path.join(outdir, f'src_lbo{ev_idx}.ply'), V_src, F_src, ev_src)
        write_q(os.path.join(outdir, f'tgt_lbo{ev_idx}_xfer.ply'), V_tgt, F_tgt, ev_xfer)
        write_q(os.path.join(outdir, f'tgt_lbo{ev_idx}_native.ply'), V_tgt, F_tgt, ev_native)

        # Correlazione xfer vs native
        std_x = ev_xfer.std(); std_n = ev_native.std()
        corr = (float(np.corrcoef(ev_xfer, ev_native)[0, 1])
                if std_x > 1e-10 and std_n > 1e-10 else float('nan'))
        verdict = "buono" if abs(corr) > 0.6 else ("parziale" if abs(corr) > 0.3 else "SCARSO")
        lines.append(f"  LBO eigvec {ev_idx}: corr(xfer, native)={corr:.3f}  [{verdict}]")

    # 3d) Indicatore regione: segno di x (split emisfero)
    hemi_src = np.sign(V_src[:, 0]).astype(float)
    hemi_xfer = hemi_src[p2p]
    write_q(os.path.join(outdir, 'src_hemi_x.ply'), V_src, F_src, hemi_src)
    write_q(os.path.join(outdir, 'tgt_hemi_x_xfer.ply'), V_tgt, F_tgt, hemi_xfer)
    lines.append("  indicatore emisfero x: salvato in src_hemi_x, tgt_hemi_x_xfer")

    return lines


# ── Test 5: direzione p2p ─────────────────────────────────────────────────────

def direction_test(p2p, n_src, n_tgt):
    lines = ["\n── Test 5: direzione a/b ──"]
    ok_len = len(p2p) == n_tgt
    ok_max = p2p.max() < n_src
    lines += [
        f"  len(p2p)={len(p2p)}  n_target={n_tgt}  → {'OK' if ok_len else 'ERRORE'}",
        f"  p2p.max()={p2p.max()}  n_source={n_src}  → {'OK' if ok_max else 'ERRORE'}",
        "  Direzione: " + (
            "CORRETTA (a=src, b=tgt)" if (ok_len and ok_max)
            else "*** INVERTITA — a/b scambiati ***"
        ),
    ]
    return lines, ok_len and ok_max


# ── Conclusione ───────────────────────────────────────────────────────────────

def make_conclusion(rep1, rep2, diag_init, diag_final, tr_init, dir_ok):
    suspects = []
    multi = rep1['n_components'] > 1 or rep2['n_components'] > 1
    if multi:
        suspects.append(
            "(a) MESH NON INTEGRA / MULTI-COMPONENTE — causa radice più probabile;\n"
            "     lo spettro LBO mescola componenti disconnesse → descrittori inconsistenti"
        )
    if not dir_ok:
        suspects.append("(d) DIREZIONE a/b INVERTITA — scambiare gli argomenti mesh1/mesh2")
    if tr_init < 0:
        suspects.append(
            "(b) FLIP DI SIMMETRIA — tr(C_init)<0; la mappa specchia l'atrio;\n"
            "     aumentare --w-orient o aggiungere un termine orientation-preserving"
        )
    if diag_final < 0.5:
        suspects.append(
            "(e) DESCRITTORI DEBOLI — diagonalità fmap finale < 50%;\n"
            "     possibili cause: mesh troppo diverse, K troppo piccolo, pesi sbagliati"
        )
    loops1, loops2 = rep1['n_boundary_loops'], rep2['n_boundary_loops']
    if loops1 != loops2:
        suspects.append(
            f"(c) ANATOMIE NON COINCIDENTI — loop bordo diversi (mesh1={loops1}, mesh2={loops2});\n"
            "     le aperture non matchano → regioni senza corrispondenza"
        )
    if not suspects:
        suspects.append(
            "nessun segnale critico rilevato automaticamente;\n"
            "verificare visivamente i PLY delle coordinate e degli autovettori LBO in MeshLab"
        )

    sep = "═" * 55
    lines = [f"\n{sep}", "CONCLUSIONE — SOSPETTATI:"]
    for i, s in enumerate(suspects, 1):
        lines.append(f"  {i}. {s}")
    lines.append(sep)
    return lines


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Diagnosi qualità functional map tra mesh cardiache (solo lettura, niente correzioni).'
    )
    parser.add_argument('mesh1', help='Target (.ply) – FA (b)')
    parser.add_argument('mesh2', help='Source (.ply) – normale (a)')
    parser.add_argument('--outdir',   default='diag_out', help='Cartella output (default: diag_out)')
    parser.add_argument('--k',        type=int,   default=120,  help='Spettro LBO (default 120)')
    parser.add_argument('--use-k',    type=int,   default=30,   help='use_k iniziale fmap (default 30)')
    parser.add_argument('--n-times',  type=int,   default=100,  help='Passi temporali HKS/WKS (default 100)')
    parser.add_argument('--w-descr',  type=float, default=1.0)
    parser.add_argument('--w-lbc',    type=float, default=1e-2)
    parser.add_argument('--w-mult',   type=float, default=1e-1)
    parser.add_argument('--icp-nit',  type=int,   default=10)
    parser.add_argument('--zo-nit',   type=int,   default=6)
    parser.add_argument('--zo-step',  type=int,   default=5)
    parser.add_argument('--no-sinkhorn', action='store_true')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    report = []

    def log(lines):
        if isinstance(lines, str):
            lines = [lines]
        for l in lines:
            print(l)
        report.extend(lines)

    # ── Lettura ────────────────────────────────────────────────────────────────
    print(f"Leggo mesh1 (target/FA): {args.mesh1}")
    props1, F1 = read_ply(args.mesh1)
    V1 = np.column_stack([props1['x'], props1['y'], props1['z']])
    print(f"Leggo mesh2 (source/normale): {args.mesh2}")
    props2, F2 = read_ply(args.mesh2)
    V2 = np.column_stack([props2['x'], props2['y'], props2['z']])

    # ── Test 1: Integrità ──────────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("TEST 1 — INTEGRITÀ MESH")
    log("=" * 55)
    rep1 = mesh_integrity(V1, F1)
    rep2 = mesh_integrity(V2, F2)
    log(format_integrity("mesh1/target/FA", rep1))
    log(format_integrity("mesh2/source/normale", rep2))
    loops_match = rep1['n_boundary_loops'] == rep2['n_boundary_loops']
    log(f"\n  Loop bordo: mesh1={rep1['n_boundary_loops']}, mesh2={rep2['n_boundary_loops']}"
        + ("  [ok — stesse aperture anatomiche]" if loops_match
           else "  [WARN — anatomie non coincidenti]"))

    # ── Pipeline geomfum ───────────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("PIPELINE GEOMFUM")
    log("=" * 55)

    mesh_src = TriangleMesh(V2, F2)   # a = source = normale
    mesh_tgt = TriangleMesh(V1, F1)   # b = target = FA

    mesh_src.vertices -= mesh_src.vertices.mean(axis=0)
    mesh_tgt.vertices -= mesh_tgt.vertices.mean(axis=0)
    V_src = np.array(mesh_src.vertices)   # centrate
    V_tgt = np.array(mesh_tgt.vertices)

    K = args.k
    log(f"LBO spectrum K={K} ...")
    mesh_src.laplacian.find_spectrum(spectrum_size=K, set_as_basis=True)
    mesh_tgt.laplacian.find_spectrum(spectrum_size=K, set_as_basis=True)

    log(f"Descrittori HKS+WKS (n_times={args.n_times}) ...")
    hks = HeatKernelSignature.from_registry(n_domain=args.n_times)
    wks = WaveKernelSignature.from_registry(n_domain=args.n_times)
    pipeline = DescriptorPipeline([hks, wks, L2InnerNormalizer()])
    descr_src = pipeline.apply(mesh_src)
    descr_tgt = pipeline.apply(mesh_tgt)

    use_k = args.use_k
    mesh_src.basis.use_k = use_k
    mesh_tgt.basis.use_k = use_k

    log(f"FunctionalMap (use_k={use_k}) ...")
    solver = FunctionalMap(factor_builders=[
        SDPFactorBuilder(weight=args.w_descr),
        LBCFactorBuilder(weight=args.w_lbc),
        MultFactorBuilder(weight=args.w_mult),
    ])
    fmap_init = solver(mesh_src.basis, mesh_tgt.basis, descr_src, descr_tgt)

    # ── Test 2: fmap iniziale ──────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("TEST 2 — ISPEZIONE FMAP")
    log("=" * 55)
    k_sq = min(fmap_init.shape)
    tr_init = float(np.trace(fmap_init[:k_sq, :k_sq]))
    diag_init = diagonality(fmap_init)
    flip_warn = "  *** FLIP RILEVATO — tr<0 ***" if tr_init < 0 else "  [ok]"
    log(f"  fmap iniziale: shape={fmap_init.shape},  tr={tr_init:.3f},  diag={diag_init:.3f}{flip_warn}")
    png_init = os.path.join(args.outdir, 'fmap_init.png')
    save_fmap_png(fmap_init, png_init, f'fmap iniziale  diag={diag_init:.2f}  tr={tr_init:.2f}')
    if HAS_MPL:
        log(f"  → {png_init}")

    # Refinement
    log(f"ICP (nit={args.icp_nit}) + ZoomOut (nit={args.zo_nit}, step={args.zo_step}) ...")
    fmap = fmap_init.copy()
    if args.icp_nit > 0:
        fmap = IcpRefiner(nit=args.icp_nit)(fmap, mesh_src.basis, mesh_tgt.basis)
    fmap = ZoomOut(nit=args.zo_nit, step=args.zo_step)(fmap, mesh_src.basis, mesh_tgt.basis)

    k_sq2 = min(fmap.shape)
    tr_final = float(np.trace(fmap[:k_sq2, :k_sq2]))
    diag_final = diagonality(fmap)
    log(f"  fmap finale  : shape={fmap.shape},  tr={tr_final:.3f},  diag={diag_final:.3f}")
    png_final = os.path.join(args.outdir, 'fmap_final.png')
    save_fmap_png(fmap, png_final, f'fmap finale  diag={diag_final:.2f}  tr={tr_final:.2f}')
    if HAS_MPL:
        log(f"  → {png_final}")

    # ── P2P ───────────────────────────────────────────────────────────────────
    log("Calcolo p2p ...")
    if args.no_sinkhorn:
        p2p = P2pFromFmConverter()(fmap, mesh_src.basis, mesh_tgt.basis)
    else:
        try:
            p2p = SinkhornP2pFromFmConverter()(fmap, mesh_src.basis, mesh_tgt.basis)
        except Exception as exc:
            log(f"  [warn] Sinkhorn fallito ({exc}), uso NN")
            p2p = P2pFromFmConverter()(fmap, mesh_src.basis, mesh_tgt.basis)
    p2p = np.asarray(p2p, dtype=np.int64)

    # ── Test 5: direzione ─────────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("TEST 5 — DIREZIONE A/B")
    log("=" * 55)
    dir_lines, dir_ok = direction_test(p2p, n_src=len(V_src), n_tgt=len(V_tgt))
    log(dir_lines)

    # ── Test 3: funzioni geometriche ──────────────────────────────────────────
    log("\n" + "=" * 55)
    log("TEST 3 — FUNZIONI GEOMETRICHE TRASFERITE")
    log("=" * 55)
    geom_lines = geometric_transfer_tests(
        V_src, F2, V_tgt, F1, p2p,
        mesh_src.basis, mesh_tgt.basis,
        args.outdir,
    )
    log(geom_lines)
    log(f"  PLY salvati in {args.outdir}/")

    # ── Test 4: round-trip ────────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("TEST 4 — ROUND-TRIP BIJETTIVITÀ")
    log("=" * 55)
    try:
        inv_p2p = np.asarray(
            P2pFromFmConverter()(fmap.T, mesh_tgt.basis, mesh_src.basis),
            dtype=np.int64,
        )
        k_rt = inv_p2p[p2p]
        n_tgt_n = len(V_tgt)
        rt_err = np.linalg.norm(V_tgt[k_rt] - V_tgt[np.arange(n_tgt_n)], axis=1)
        diag_len = float(np.linalg.norm(V_tgt.max(axis=0) - V_tgt.min(axis=0)))
        thresh = 0.05 * diag_len
        pct_good = 100.0 * float(np.mean(rt_err < thresh))
        rt_verdict = (
            "BUONO (>70%)" if pct_good > 70
            else ("ACCETTABILE" if pct_good > 50 else "SCARSO — mappa inaffidabile")
        )
        log([
            f"  errore RT medio  : {rt_err.mean():.3f} mm,  mediana: {np.median(rt_err):.3f} mm",
            f"  RT 90° percentile: {np.percentile(rt_err, 90):.3f} mm",
            f"  diagonale bbox   : {diag_len:.1f} mm,  soglia 5%: {thresh:.1f} mm",
            f"  dentro soglia    : {pct_good:.1f}%  [{rt_verdict}]",
        ])
    except Exception as exc:
        log(f"  FALLITO: {exc}")

    # ── Conclusione ────────────────────────────────────────────────────────────
    log("\n" + "=" * 55)
    log("CONCLUSIONE")
    log("=" * 55)
    log(make_conclusion(rep1, rep2, diag_init, diag_final, tr_init, dir_ok))

    # ── Scrivi referto ─────────────────────────────────────────────────────────
    report_path = os.path.join(args.outdir, 'report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    print(f"\nReferto completo: {report_path}")
    print(f"PNG fmap:         {args.outdir}/fmap_init.png  {args.outdir}/fmap_final.png")
    print(f"PLY funzioni:     {args.outdir}/src_*.ply  {args.outdir}/tgt_*_xfer.ply")


if __name__ == '__main__':
    main()
