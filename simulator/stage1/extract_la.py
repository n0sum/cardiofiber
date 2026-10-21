import pyvista as pv, numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # radice del repository
DIR  = ROOT / "data" / "stage1"
MASK = ROOT / "simulator" / "stage1" / "la_vertex_mask.npy"

# Verso dei triangoli: AugmentA inverte il verso di tutti i triangoli rispetto
# all'SSM (elem == ssm_tri[:, ::-1]), stesso insieme e stesso ordine. La
# connettivita' di riferimento per tutte le forme e' LA_endo_with_fiber.elem,
# mai representer/cells.

LA_REGION = 0   # <-- se il conteggio vene dice component_1, metti 1

mesh = pv.read(DIR / "ssm_mean_surface.vtk")

# timbra ogni vertice col suo indice originale (0..39868) PRIMA di dividere
mesh.point_data["orig_id"] = np.arange(mesh.n_points, dtype=np.int64)

conn = mesh.connectivity()
rid  = np.asarray(conn.cell_data["RegionId"])

la_cells = np.where(rid == LA_REGION)[0]
la_surf  = conn.extract_cells(la_cells).extract_surface(algorithm="dataset_surface")   # default di pyvista 0.48.4, fissato

# gli indici originali sopravvivono nel campo che abbiamo creato noi
la_point_ids = np.asarray(la_surf.point_data["orig_id"]).astype(np.int64)
np.save(MASK, np.sort(la_point_ids))

la_surf.save(DIR / "ssm_mean_LA.vtk")
print(f"LA: {la_surf.n_points} vertici, {la_surf.n_cells} triangoli")
print(f"maschera salvata: la_vertex_mask.npy ({la_point_ids.size} indici)")
print(f"range indici: [{la_point_ids.min()}, {la_point_ids.max()}]  max atteso < 39869")
