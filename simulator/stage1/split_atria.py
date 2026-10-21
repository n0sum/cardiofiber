import pyvista as pv, numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # radice del repository
DIR  = ROOT / "data" / "stage1"

mesh = pv.read(DIR / "ssm_mean_surface.vtk")
conn = mesh.connectivity()                 # assegna RegionId
rid  = np.asarray(conn.cell_data["RegionId"])
nreg = int(rid.max()) + 1
print(f"Componenti connesse trovate: {nreg}\n")

for r in range(nreg):
    body = conn.threshold([r - 0.5, r + 0.5], scalars="RegionId")
    surf = body.extract_surface(algorithm="dataset_surface")   # default di pyvista 0.48.4, fissato
    c = np.asarray(surf.points)
    print(f"RegionId {r}:  vertici={surf.n_points:6d}   "
          f"centro_x={c[:,0].mean():8.2f}   "
          f"x∈[{c[:,0].min():.1f}, {c[:,0].max():.1f}]")
    surf.save(DIR / f"component_{r}.vtk")

print(f"\nScritti component_0.vtk, component_1.vtk, ... in {DIR}")
