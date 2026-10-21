import h5py, numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # radice del repository
H5   = ROOT / "data" / "external" / "biatrial_model_completeFit.h5"

def walk(name, obj):
    if isinstance(obj, h5py.Dataset):
        print(f"[dataset] {name}  shape={obj.shape}  dtype={obj.dtype}")
    else:
        print(f"[group]   {name}")

with h5py.File(H5, "r") as f:
    f.visititems(walk)
