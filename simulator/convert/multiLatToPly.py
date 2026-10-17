#!/usr/bin/env python3
"""
Esporta mesh + N mappe LAT in un PLY binario leggibile con happly.

Le LAT diventano proprieta' per-vertice chiamate lat_0, lat_1, ... ; i nomi dei
siti di pacing finiscono nei commenti dell'header.

Legge le fonti CANONICHE del bundle (file testuali openCARP + npz), non i .vtk:
i .vtk passati per pyvista possono attraversare filtri che riordinano le celle
(extract_surface), rompendo l'allineamento faccia <-> dato.

I percorsi si configurano nelle costanti in cima al file.

Lato C++:
    v.getProperty<float>("lat_k")  ->  vQualities[k]
"""

from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
ROOT      = Path(__file__).resolve().parents[2]           # radice del repository
MESH_BASE = ROOT / "data" / "augmenta"  / "LA_endo_with_fiber"   # senza estensione
NPZ_PATH  = ROOT / "data" / "reference" / "LAT_all.npz"
OUT_PATH  = ROOT / "data" / "reference" / "LA_LAT.ply"
# --------------------------------------------------------------------------


def read_pts(path):
    n = int(open(path).readline())
    pts = np.loadtxt(path, skiprows=1, dtype=np.float64)
    if pts.shape != (n, 3):
        raise SystemExit(f"{path}: attesi {n} punti, trovati {pts.shape}")
    return pts


def read_elem(path):
    n = int(open(path).readline())
    # righe: "Tr v0 v1 v2 tag"  -> saltiamo la colonna di testo e il tag
    tri = np.loadtxt(path, skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    if tri.shape != (n, 3):
        raise SystemExit(f"{path}: attesi {n} triangoli, trovati {tri.shape}")
    return tri


def main():
    mesh_base = str(MESH_BASE)

    pts = read_pts(mesh_base + ".pts")
    tri = read_elem(mesh_base + ".elem")
    n_verts, n_faces = pts.shape[0], tri.shape[0]

    data = np.load(NPZ_PATH)
    names = list(data.files)
    maps = []
    for name in names:
        a = np.asarray(data[name]).ravel()
        if a.size != n_verts:
            raise SystemExit(f"mappa '{name}': {a.size} valori, attesi {n_verts}")
        maps.append(a.astype(np.float32))

    if tri.min() < 0 or tri.max() >= n_verts:
        raise SystemExit("indici di vertice fuori range nel .elem")

    # ---- vertici: x y z lat_0 lat_1 ... (interleaved, float32) -------------
    vprops = ["x", "y", "z"] + [f"lat_{k}" for k in range(len(maps))]
    vdtype = np.dtype([(p, "<f4") for p in vprops])
    vrec = np.empty(n_verts, dtype=vdtype)
    vrec["x"] = pts[:, 0]
    vrec["y"] = pts[:, 1]
    vrec["z"] = pts[:, 2]
    for k, m in enumerate(maps):
        vrec[f"lat_{k}"] = m

    # ---- facce: uchar 3 + 3 int32 (dtype packed, niente padding) ----------
    fdtype = np.dtype([("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")])
    frec = np.empty(n_faces, dtype=fdtype)
    frec["n"] = 3
    frec["a"] = tri[:, 0]
    frec["b"] = tri[:, 1]
    frec["c"] = tri[:, 2]

    header = ["ply", "format binary_little_endian 1.0"]
    header.append(f"comment mesh {mesh_base}")
    header.append(f"comment lat source {NPZ_PATH}")
    for k, name in enumerate(names):
        header.append(f"comment lat_{k} {name}")
    header.append(f"element vertex {n_verts}")
    header += [f"property float {p}" for p in vprops]
    header.append(f"element face {n_faces}")
    header.append("property list uchar int vertex_indices")
    header.append("end_header")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        vrec.tofile(f)
        frec.tofile(f)

    print(f"scritto {OUT_PATH}")
    print(f"  vertici {n_verts}   facce {n_faces}   mappe {len(maps)}")
    for k, name in enumerate(names):
        print(f"    lat_{k} = {name}")
    print(f"  estensione (mm): {np.round(pts.max(0) - pts.min(0), 2)}")


if __name__ == "__main__":
    main()