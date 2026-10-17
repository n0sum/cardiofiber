#!/usr/bin/env python3
"""
Esporta mesh + fibre ground truth in un PLY binario leggibile con happly.

Le fibre diventano proprieta' per-FACCIA (fx, fy, fz) insieme al tag di
regione (elem_tag).

Legge le fonti CANONICHE del bundle (file testuali openCARP), non i .vtk:
i .vtk passati per pyvista possono attraversare filtri che riordinano le
celle, rompendo l'allineamento faccia <-> fibra. A differenza di un
disallineamento sui vertici, questo NON produrrebbe nessun errore visibile:
solo un campo di fibre sbagliato.

Lato C++:
    face.getProperty<float>("fx")  ->  ground truth per faccia
"""

from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
ROOT      = Path(__file__).resolve().parents[2]           # radice del repository
MESH_BASE = ROOT / "data" / "augmenta"  / "LA_endo_with_fiber"   # senza estensione
OUT_PATH  = ROOT / "data" / "reference" / "LA_fibers_gt.ply"
# --------------------------------------------------------------------------


def read_pts(path):
    n = int(open(path).readline())
    pts = np.loadtxt(path, skiprows=1, dtype=np.float64)
    if pts.shape != (n, 3):
        raise SystemExit(f"{path}: attesi {n} punti, trovati {pts.shape}")
    return pts


def read_elem(path):
    """Righe: 'Tr v0 v1 v2 tag'. Qui serve anche il tag."""
    n = int(open(path).readline())
    cols = np.loadtxt(path, skiprows=1, usecols=(1, 2, 3, 4), dtype=np.int64)
    if cols.shape != (n, 4):
        raise SystemExit(f"{path}: attesi {n} triangoli, trovati {cols.shape}")
    return cols[:, :3], cols[:, 3].astype(np.int32)


def read_lon(path, n_faces):
    """Prima riga: 1 (solo f) o 2 (ortotropo, f e t). Qui serve solo f."""
    with open(path) as f:
        kind = int(f.readline())
    if kind not in (1, 2):
        raise SystemExit(f"{path}: prima riga '{kind}', attesi 1 o 2")

    vals = np.loadtxt(path, skiprows=1, dtype=np.float64)
    if vals.ndim == 1:
        vals = vals.reshape(1, -1)
    if vals.shape[0] != n_faces:
        raise SystemExit(f"{path}: {vals.shape[0]} righe, attese {n_faces}")
    if vals.shape[1] != 3 * kind:
        raise SystemExit(f"{path}: {vals.shape[1]} colonne, attese {3 * kind}")

    return vals[:, :3]          # la trasversale t si ricava dalla normale


def main():
    mesh_base = str(MESH_BASE)

    pts            = read_pts(mesh_base + ".pts")
    tri, elem_tag  = read_elem(mesh_base + ".elem")
    fib            = read_lon(mesh_base + ".lon", tri.shape[0])
    n_verts, n_faces = pts.shape[0], tri.shape[0]

    if tri.min() < 0 or tri.max() >= n_verts:
        raise SystemExit("indici di vertice fuori range nel .elem")

    # ---- controllo normalizzazione ---------------------------------------
    norms = np.linalg.norm(fib, axis=1)
    if norms.min() < 1e-6:
        raise SystemExit(f"{(norms < 1e-6).sum()} fibre di modulo nullo nel .lon")
    fib = (fib / norms[:, None]).astype(np.float32)

    # ---- quanto le fibre sono gia' tangenti alle facce -------------------
    e1 = pts[tri[:, 1]] - pts[tri[:, 0]]
    e2 = pts[tri[:, 2]] - pts[tri[:, 0]]
    nrm = np.cross(e1, e2)
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    dev = np.degrees(np.arcsin(np.clip(np.abs(np.einsum("ij,ij->i", fib, nrm)), 0, 1)))

    # ---- vertici: x y z ---------------------------------------------------
    vdtype = np.dtype([("x", "<f4"), ("y", "<f4"), ("z", "<f4")])
    vrec = np.empty(n_verts, dtype=vdtype)
    vrec["x"], vrec["y"], vrec["z"] = pts[:, 0], pts[:, 1], pts[:, 2]

    # ---- facce: uchar 3 + 3 int32 + fx fy fz + elem_tag -------------------
    fdtype = np.dtype([("n", "u1"),
                       ("a", "<i4"), ("b", "<i4"), ("c", "<i4"),
                       ("fx", "<f4"), ("fy", "<f4"), ("fz", "<f4"),
                       ("elem_tag", "<i4")])
    frec = np.empty(n_faces, dtype=fdtype)
    frec["n"] = 3
    frec["a"], frec["b"], frec["c"] = tri[:, 0], tri[:, 1], tri[:, 2]
    frec["fx"], frec["fy"], frec["fz"] = fib[:, 0], fib[:, 1], fib[:, 2]
    frec["elem_tag"] = elem_tag

    header = ["ply", "format binary_little_endian 1.0"]
    header.append(f"comment mesh {mesh_base}")
    header.append("comment ground truth: fibre per faccia dal .lon")
    header.append(f"element vertex {n_verts}")
    header += ["property float x", "property float y", "property float z"]
    header.append(f"element face {n_faces}")
    header.append("property list uchar int vertex_indices")
    header += ["property float fx", "property float fy", "property float fz"]
    header.append("property int elem_tag")
    header.append("end_header")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        vrec.tofile(f)
        frec.tofile(f)

    print(f"scritto {OUT_PATH}")
    print(f"  vertici {n_verts}   facce {n_faces}")
    print(f"  estensione (mm): {np.round(pts.max(0) - pts.min(0), 2)}")
    tags, counts = np.unique(elem_tag, return_counts=True)
    print("  tag di regione:")
    for t, c in zip(tags, counts):
        print(f"    {t}: {c} facce")
    print(f"  deviazione fibra dal piano della faccia: "
          f"mediana {np.median(dev):.3f}deg   max {dev.max():.3f}deg")


if __name__ == "__main__":
    main()