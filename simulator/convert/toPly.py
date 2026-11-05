#!/usr/bin/env python3
"""
Da openCARP a PLY binario: i due ingressi di Tesi0 per uno scenario.

    data/reference/LA_LAT_NOME.ply   per vertice x y z (mm) e lat_0 ... lat_5 (ms):
                                     un act.igb per riga di pacing_sites.csv
    data/reference/LA_GT_NOME.ply    per faccia fibra fx fy fz, CV_L CV_T (mm/ms)
                                     ed elem_tag: la verita' per la valutazione

Coordinate e fibre vengono da data/augmenta/ (mm), triangoli e tag da
MESH_BASE.elem, la mesh in um su cui ha girato run_sims.py; le CV per tag
dal CSV dello scenario. Si leggono i file testuali di openCARP e non i .vtk:
un filtro che riordina le celle disallineerebbe facce e fibre senza errori.

Uso, dalla radice del repository:
    python simulator/convert/toPly.py NOME CSV MESH_BASE
    python simulator/convert/toPly.py aniso_zone_50 experiments/region_cv/aniso_random.csv data/reference/aniso_zone_50/LA_um
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT     = Path(__file__).resolve().parents[2]                    # radice del repository
AUGMENTA = ROOT / "data" / "augmenta" / "LA_endo_with_fiber"      # senza estensione
SITES    = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
SIMS_DIR = ROOT / "data" / "reference" / "sims"                   # NOME/run_k/act.igb


def write_ply(path, vertex, face):
    """vertex e face sono array strutturati; face comincia con n, a, b, c, la lista dei vertici."""
    ply_type = {np.dtype("<f4"): "float", np.dtype("<i4"): "int"}
    header = ["ply", "format binary_little_endian 1.0", f"element vertex {len(vertex)}"]
    header += [f"property {ply_type[vertex.dtype[p]]} {p}" for p in vertex.dtype.names]
    header += [f"element face {len(face)}", "property list uchar int vertex_indices"]
    header += [f"property {ply_type[face.dtype[p]]} {p}" for p in face.dtype.names[4:]]
    header += ["end_header"]
    with open(path, "wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        vertex.tofile(f)
        face.tofile(f)


def main():
    name, cv_csv, mesh_base = sys.argv[1], sys.argv[2], sys.argv[3]

    pts = np.loadtxt(f"{AUGMENTA}.pts", skiprows=1)
    fib = np.loadtxt(f"{AUGMENTA}.lon", skiprows=1)[:, :3]          # fibra; la sheet non serve
    elem = np.loadtxt(f"{mesh_base}.elem", skiprows=1, usecols=(1, 2, 3, 4), dtype=np.int64)
    tri, tag = elem[:, :3], elem[:, 3]

    maps = range(len(list(csv.DictReader(open(SITES)))))           # la mappa k e' la riga k
    lat = []
    for k in maps:
        a = np.fromfile(SIMS_DIR / name / f"run_{k}" / "act.igb", dtype="<f4", offset=1024)   # header di 1024 byte
        if not (a >= 0).all():                                       # FIM non convergito: LAT = -1
            raise SystemExit(f"mappa {k}: {int((~(a >= 0)).sum())} nodi con LAT negativa o non finita")
        lat.append(a)

    fib = (fib / np.linalg.norm(fib, axis=1)[:, None]).astype(np.float32)
    cv = {int(r["region_tag"]): (float(r["CV_L"]), float(r["CV_T"]))
          for r in csv.DictReader(l for l in open(cv_csv) if not l.startswith("#"))}
    cv_l = np.array([cv[t][0] for t in tag], dtype=np.float32)
    cv_t = np.array([cv[t][1] for t in tag], dtype=np.float32)

    xyz = [pts[:, 0], pts[:, 1], pts[:, 2]]
    xyz_names = [("x", "<f4"), ("y", "<f4"), ("z", "<f4")]
    triangles = [np.full(len(tri), 3), tri[:, 0], tri[:, 1], tri[:, 2]]
    triangle_names = [("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")]
    out = ROOT / "data" / "reference"

    write_ply(out / f"LA_LAT_{name}.ply",
              np.rec.fromarrays(xyz + lat, dtype=xyz_names + [(f"lat_{k}", "<f4") for k in maps]),
              np.rec.fromarrays(triangles, dtype=triangle_names))
    write_ply(out / f"LA_GT_{name}.ply",
              np.rec.fromarrays(xyz, dtype=xyz_names),
              np.rec.fromarrays(triangles + [fib[:, 0], fib[:, 1], fib[:, 2], cv_l, cv_t, tag],
                                dtype=triangle_names + [("fx", "<f4"), ("fy", "<f4"), ("fz", "<f4"),
                                                        ("CV_L", "<f4"), ("CV_T", "<f4"), ("elem_tag", "<i4")]))


if __name__ == "__main__":
    main()
