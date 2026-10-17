#!/usr/bin/env python3
"""
Transfer LAT (quality) from mesh2 (source, FA) to mesh1 (target, normale)
using Functional Maps via geomfum.

Convenzione (identica a diag_stages.py, verificata su main.cpp):
    mesh1elaborata.ply = 2-LA.mesh    = NORMALE = TARGET = b
    mesh2elaborata.ply = 2-LA-FA.mesh = FA      = SOURCE = a
    mesh_src = TriangleMesh(V2, F2)   (a, dona il LAT)
    mesh_tgt = TriangleMesh(V1, F1)   (b, mantiene geometria/topologia)
    p2p = converter(fmap, basis_src, basis_tgt)  shape=[n_tgt], values in [0, n_src)
    lat_transferred = lat2[p2p]  (pullback per indice, nessuna ricostruzione spettrale)

Config di produzione (hardcoded), validata in diag_stages.py come variante
"(g) landmark+ZO(step=1)": corr_x≈0.68, corr_ev1≈0.64-0.74, nessun flip vero.
NON è più un esperimento: niente sweep di pesi/varianti ad ogni run, niente
ICP, niente termine orientation-preserving (tr(fmap)<0 è un proxy di flip
inaffidabile tra atri di forma diversa - vedi diag_stages.py).

Usage:
    python lat_transfer.py mesh1.ply mesh2.ply output.ply

Output PLY ha la stessa geometria/topologia di mesh1 con:
  quality             - LAT nativo di mesh1 (normale/target)
  quality_transferred - LAT di mesh2 (FA/source) trasferito via pullback
  oor_flag            - richiesto dal contratto PLY di Mesh::importPlyDualLat;
                         sempre 0 sotto pullback puro (non può esserci overshoot
                         perché lat_transferred è un sottoinsieme di lat2 per
                         costruzione), mantenuto per compatibilità con il C++.
"""

import sys
import numpy as np

from landmarks import get_matched_landmarks

try:
    from geomfum.shape import TriangleMesh
    from geomfum.descriptor.spectral import (
        HeatKernelSignature, WaveKernelSignature,
        LandmarkHeatKernelSignature, LandmarkWaveKernelSignature,
    )
    from geomfum.descriptor.pipeline import DescriptorPipeline, L2InnerNormalizer
    from geomfum.functional_map import FunctionalMap, SDPFactorBuilder, LBCFactorBuilder, MultFactorBuilder
    from geomfum.refine import ZoomOut
    from geomfum.convert import SinkhornP2pFromFmConverter, P2pFromFmConverter
except ImportError as exc:
    sys.exit(
        f"Import error: {exc}\n"
        "Install:\n"
        "  pip install 'geomfum@git+https://github.com/3diglab/geomfum.git@main'\n"
        "  pip install pot robust_laplacian geomstats"
    )


# ─── Config di produzione (validata, vedi diag_stages.py variante "(g)") ──────

LBO_K = 120
FMAP_USE_K = 30
N_TIMES = 100
N_TIMES_LANDMARKS = 3
W_DESCR = 1.0
W_LBC = 1e-2
W_MULT = 1e-1
ZOOMOUT_NIT = 30
ZOOMOUT_STEP = 1


# ─── PLY I/O ─────────────────────────────────────────────────────────────────

def read_ply(path):
    """Read ASCII PLY. Returns (props dict of float arrays, faces int32 array)."""
    with open(path, 'r') as f:
        lines = f.readlines()

    end = next(i for i, l in enumerate(lines) if l.strip() == 'end_header')
    header = [l.strip() for l in lines[:end]]
    body = lines[end + 1:]

    n_verts = 0
    n_faces = 0
    vert_props = []
    in_vert = False
    in_face = False

    for line in header:
        if line.startswith('element vertex'):
            n_verts = int(line.split()[-1])
            in_vert, in_face = True, False
        elif line.startswith('element face'):
            n_faces = int(line.split()[-1])
            in_vert, in_face = False, True
        elif line.startswith('property') and in_vert and 'list' not in line:
            vert_props.append(line.split()[-1])

    V = np.array([[float(v) for v in body[i].split()] for i in range(n_verts)])
    F = np.array(
        [list(map(int, body[n_verts + i].split()[1:4])) for i in range(n_faces)],
        dtype=np.int32,
    )
    props = {name: V[:, i] for i, name in enumerate(vert_props)}
    return props, F


def write_ply(path, props, faces):
    """Write ASCII PLY with given vertex property dict and face array."""
    n_verts = len(next(iter(props.values())))
    with open(path, 'w') as out:
        out.write("ply\nformat ascii 1.0\ncomment Written by lat_transfer.py\n")
        out.write(f"element vertex {n_verts}\n")
        for name in props:
            out.write(f"property float {name}\n")
        out.write(f"element face {len(faces)}\n")
        out.write("property list uchar int vertex_indices\nend_header\n")
        data = np.column_stack(list(props.values()))
        for row in data:
            out.write(' '.join(f'{v:.7g}' for v in row) + '\n')
        for f in faces:
            out.write(f"3 {f[0]} {f[1]} {f[2]}\n")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    argv = sys.argv[1:]
    repick = '--repick' in argv
    argv = [a for a in argv if a != '--repick']
    if len(argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} mesh1.ply mesh2.ply output.ply [--repick]")
    mesh1_path, mesh2_path, output_path = argv

    # ── Read ──────────────────────────────────────────────────────────────────
    print(f"Reading target (mesh1/normale): {mesh1_path}")
    props1, F1 = read_ply(mesh1_path)
    V1 = np.column_stack([props1['x'], props1['y'], props1['z']])
    lat1 = props1['quality']
    print(f"  {len(V1)} vertices, {len(F1)} faces")

    print(f"Reading source (mesh2/FA): {mesh2_path}")
    props2, F2 = read_ply(mesh2_path)
    V2 = np.column_stack([props2['x'], props2['y'], props2['z']])
    lat2 = props2['quality']
    print(f"  {len(V2)} vertices, {len(F2)} faces")

    # ── Landmark (da cache, o picking interattivo pyvista se mancante) ───────
    lm1, lm2 = get_matched_landmarks(mesh1_path, V1, F1, mesh2_path, V2, F2, repick=repick)
    print(f"  landmark: {len(lm1)} punti per mesh")

    # ── Build geomfum meshes ──────────────────────────────────────────────────
    print("Building geomfum meshes and centering...")
    mesh_src = TriangleMesh(V2, F2)   # a = source = FA
    mesh_tgt = TriangleMesh(V1, F1)   # b = target = normale

    mesh_src.landmark_indices = np.asarray(lm2, dtype=np.int64)
    mesh_tgt.landmark_indices = np.asarray(lm1, dtype=np.int64)

    mesh_src.vertices -= mesh_src.vertices.mean(axis=0)
    mesh_tgt.vertices -= mesh_tgt.vertices.mean(axis=0)

    # ── Spectral basis ────────────────────────────────────────────────────────
    print(f"Computing LBO spectrum (K={LBO_K}) ...")
    mesh_src.laplacian.find_spectrum(spectrum_size=LBO_K, set_as_basis=True)
    mesh_tgt.laplacian.find_spectrum(spectrum_size=LBO_K, set_as_basis=True)

    # ── Descriptors (geometry + landmarks: HKS + WKS + LandmarkHKS + LandmarkWKS) ─
    print(f"Computing descriptors (n_times={N_TIMES}, n_times_landmarks={N_TIMES_LANDMARKS}) ...")
    hks = HeatKernelSignature.from_registry(n_domain=N_TIMES)
    wks = WaveKernelSignature.from_registry(n_domain=N_TIMES)
    lm_hks = LandmarkHeatKernelSignature.from_registry(n_domain=N_TIMES_LANDMARKS)
    lm_wks = LandmarkWaveKernelSignature.from_registry(n_domain=N_TIMES_LANDMARKS)
    pipeline = DescriptorPipeline([hks, wks, lm_hks, lm_wks, L2InnerNormalizer()])
    descr_src = pipeline.apply(mesh_src)
    descr_tgt = pipeline.apply(mesh_tgt)
    print(f"  descriptor shapes: src {descr_src.shape}, tgt {descr_tgt.shape}")

    # ── Functional map: solver grezzo + ZoomOut(step=1) ───────────────────────
    # Niente ICP, niente termine orientation-preserving: config validata "(g)".
    mesh_src.basis.use_k = FMAP_USE_K
    mesh_tgt.basis.use_k = FMAP_USE_K
    print(f"Computing functional map (initial use_k={FMAP_USE_K}) ...")

    solver = FunctionalMap(factor_builders=[
        SDPFactorBuilder(weight=W_DESCR),
        LBCFactorBuilder(weight=W_LBC),
        MultFactorBuilder(weight=W_MULT),
    ])
    fmap = solver(mesh_src.basis, mesh_tgt.basis, descr_src, descr_tgt)
    print(f"  fmap shape: {fmap.shape}")

    print(f"ZoomOut (nit={ZOOMOUT_NIT}, step={ZOOMOUT_STEP}) ...")
    fmap = ZoomOut(nit=ZOOMOUT_NIT, step=ZOOMOUT_STEP)(fmap, mesh_src.basis, mesh_tgt.basis)
    print(f"  fmap final shape: {fmap.shape}")

    # ── Point-to-point correspondence (Sinkhorn, fallback NN) ─────────────────
    print("Computing p2p via Sinkhorn ...")
    try:
        p2p = SinkhornP2pFromFmConverter()(fmap, mesh_src.basis, mesh_tgt.basis)
    except Exception as exc:
        print(f"  [warn] Sinkhorn failed ({exc}), falling back to NN")
        p2p = P2pFromFmConverter()(fmap, mesh_src.basis, mesh_tgt.basis)

    p2p = np.asarray(p2p, dtype=np.int64)
    print(f"  p2p shape: {p2p.shape},  index range: [{p2p.min()}, {p2p.max()}]")

    # ── Pullback transfer: copy lat, no spectral reconstruction ──────────────
    # Preserva la piena risoluzione della sorgente e non può fare overshoot:
    # lat_transferred è per costruzione un sottoinsieme di lat2.
    lat_transferred = lat2[p2p]
    print(f"  lat2 range: [{lat2.min():.2f}, {lat2.max():.2f}] ms  "
          f"| transferred range: [{lat_transferred.min():.2f}, {lat_transferred.max():.2f}] ms")
    assert lat_transferred.min() >= lat2.min() and lat_transferred.max() <= lat2.max(), (
        "pullback-by-index invariant violated - p2p deve indicizzare solo dentro lat2"
    )

    # ── Write main output ─────────────────────────────────────────────────────
    out_props = dict(props1)
    out_props['quality_transferred'] = lat_transferred.astype(np.float32)
    # oor_flag richiesto dal contratto PLY di Mesh::importPlyDualLat; sempre 0
    # sotto pullback puro (mantenuto solo per compatibilità con il C++).
    out_props['oor_flag'] = np.zeros(len(lat_transferred), dtype=np.float32)

    print(f"Writing {output_path} ...")
    write_ply(output_path, out_props, F1)

    print("\nDone.")
    print(f"  mesh1 LAT (normale)   : [{lat1.min():.2f}, {lat1.max():.2f}] ms")
    print(f"  mesh2 LAT (FA)        : [{lat2.min():.2f}, {lat2.max():.2f}] ms")
    print(f"  transferred LAT       : [{lat_transferred.min():.2f}, {lat_transferred.max():.2f}] ms")


if __name__ == '__main__':
    main()
