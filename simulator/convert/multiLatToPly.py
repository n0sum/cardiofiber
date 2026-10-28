#!/usr/bin/env python3
"""
Esporta mesh + N mappe LAT in un PLY binario leggibile con happly.

Le LAT diventano proprieta' per-vertice chiamate lat_0, lat_1, ...

Le mappe sono gli act.igb delle run openCARP di uno scenario, una per riga di
pacing_sites.csv. Lo scenario e' il nome del CSV delle CV (es. aniso_const):
legge data/reference/sims/SCENARIO/run_k/act.igb e scrive
data/reference/LA_LAT_SCENARIO.ply.

Uso, dalla radice del repository:
    python simulator/convert/multiLatToPly.py aniso_const

Legge le fonti CANONICHE (file testuali openCARP + act.igb), non i .vtk:
i .vtk passati per pyvista possono attraversare filtri che riordinano le celle
(extract_surface), rompendo l'allineamento faccia <-> dato.

I percorsi si configurano nelle costanti in cima al file.

Lato C++:
    v.getProperty<float>("lat_k")  ->  vQualities[k]
"""

import csv
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
ROOT      = Path(__file__).resolve().parents[2]           # radice del repository
MESH_BASE = ROOT / "data" / "augmenta"  / "LA_endo_with_fiber"   # senza estensione
SITES     = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
SIMS_DIR  = ROOT / "data" / "reference" / "sims"                  # SCENARIO/run_k/act.igb
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


def read_igb(path, n):
    # header di 1024 byte, poi un float32 per vertice
    a = np.fromfile(path, dtype="<f4", offset=1024)
    if a.size != n:
        raise SystemExit(f"{path}: {a.size} valori, attesi {n}")
    return a


def main():
    scenario = sys.argv[1]
    sims_dir = SIMS_DIR / scenario
    out_path = ROOT / "data" / "reference" / f"LA_LAT_{scenario}.ply"
    mesh_base = str(MESH_BASE)

    pts = read_pts(mesh_base + ".pts")
    tri = read_elem(mesh_base + ".elem")
    n_verts, n_faces = pts.shape[0], tri.shape[0]

    maps = []
    for k in range(len(list(csv.DictReader(open(SITES))))):
        a = read_igb(sims_dir / f"run_{k}" / "act.igb", n_verts)
        if not (a >= 0).all():                         # FIM non convergito: LAT = -1
            raise SystemExit(f"mappa {k}: {int((~(a >= 0)).sum())} nodi con LAT negativa o non finita (FIM non convergito)")
        maps.append(a)

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
    header.append(f"comment lat source {sims_dir}/run_*/act.igb")
    header.append(f"element vertex {n_verts}")
    header += [f"property float {p}" for p in vprops]
    header.append(f"element face {n_faces}")
    header.append("property list uchar int vertex_indices")
    header.append("end_header")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        vrec.tofile(f)
        frec.tofile(f)


if __name__ == "__main__":
    main()