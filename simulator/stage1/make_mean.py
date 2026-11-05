import h5py, numpy as np, meshio
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # radice del repository
H5   = ROOT / "data" / "external" / "biatrial_model_completeFit.h5"
OUT  = ROOT / "data" / "stage1"

# forma media dell'SSM biatriale: model/mean e' [x0, y0, z0, x1, ...] (coincide con representer/points)
with h5py.File(H5, "r") as f:
    pts  = f["model/mean"][()].reshape(-1, 3)
    tris = f["representer/cells"][()].T.astype(np.int64)
print(f"vertici: {len(pts)}   triangoli: {len(tris)}")

OUT.mkdir(parents=True, exist_ok=True)
meshio.write_points_cells(OUT / "ssm_mean_surface.vtk", pts, [("triangle", tris)])
print(f"scritto: {OUT / 'ssm_mean_surface.vtk'}")
