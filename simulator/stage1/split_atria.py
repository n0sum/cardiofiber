import pyvista as pv, numpy as np

mesh = pv.read("ssm_mean_surface.vtk")
conn = mesh.connectivity()                 # assegna RegionId
rid  = np.asarray(conn.cell_data["RegionId"])
nreg = int(rid.max()) + 1
print(f"Componenti connesse trovate: {nreg}\n")

for r in range(nreg):
    body = conn.threshold([r - 0.5, r + 0.5], scalars="RegionId")
    surf = body.extract_surface()
    c = np.asarray(surf.points)
    print(f"RegionId {r}:  vertici={surf.n_points:6d}   "
          f"centro_x={c[:,0].mean():8.2f}   "
          f"x∈[{c[:,0].min():.1f}, {c[:,0].max():.1f}]")
    surf.save(f"component_{r}.vtk")

print("\nScritti component_0.vtk, component_1.vtk, ...")
