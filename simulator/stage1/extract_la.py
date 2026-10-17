import pyvista as pv, numpy as np

LA_REGION = 0   # <-- se il conteggio vene dice component_1, metti 1

mesh = pv.read("ssm_mean_surface.vtk")

# timbra ogni vertice col suo indice originale (0..39868) PRIMA di dividere
mesh.point_data["orig_id"] = np.arange(mesh.n_points, dtype=np.int64)

conn = mesh.connectivity()
rid  = np.asarray(conn.cell_data["RegionId"])

la_cells = np.where(rid == LA_REGION)[0]
la_surf  = conn.extract_cells(la_cells).extract_surface()

# gli indici originali sopravvivono nel campo che abbiamo creato noi
la_point_ids = np.asarray(la_surf.point_data["orig_id"]).astype(np.int64)
np.save("la_vertex_mask.npy", np.sort(la_point_ids))

la_surf.save("ssm_mean_LA.vtk")
print(f"LA: {la_surf.n_points} vertici, {la_surf.n_cells} triangoli")
print(f"maschera salvata: la_vertex_mask.npy ({la_point_ids.size} indici)")
print(f"range indici: [{la_point_ids.min()}, {la_point_ids.max()}]  max atteso < 39869")
