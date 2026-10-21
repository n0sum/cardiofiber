import h5py, numpy as np, meshio
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # radice del repository
H5   = ROOT / "data" / "external" / "biatrial_model_completeFit.h5"
OUT  = ROOT / "data" / "stage1"

OUT.mkdir(parents=True, exist_ok=True)
with h5py.File(H5, "r") as f:
    mean  = f["model/mean"][()].reshape(-1)     # (3N,)
    pref  = f["representer/points"][()]          # (3, N)
    cells = f["representer/cells"][()]           # (3, M)

N = pref.shape[1]

# due possibili ordinamenti del vettore media: verifichiamo quale combacia
cand1 = mean.reshape(N, 3)      # per-punto  [x0,y0,z0, x1,...]
cand2 = mean.reshape(3, N).T    # per-componente [x0,x1,..., y0,...]
d1 = np.abs(cand1 - pref.T).max()
d2 = np.abs(cand2 - pref.T).max()
print(f"diff per-punto: {d1:.4g}   diff per-componente: {d2:.4g}")

pts = cand1 if d1 <= d2 else cand2
tris = cells.T.astype(np.int64)
print(f"vertici: {N}   triangoli: {tris.shape[0]}")

meshio.write_points_cells(OUT / "ssm_mean_surface.vtk", pts, [("triangle", tris)])
print(f"scritto: {OUT / 'ssm_mean_surface.vtk'}")
