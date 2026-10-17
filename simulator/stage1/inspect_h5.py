import h5py, numpy as np

path = "biatrial_model_completeFit.h5"   # <-- metti qui il nome esatto del tuo file .h5

def walk(name, obj):
    if isinstance(obj, h5py.Dataset):
        print(f"[dataset] {name}  shape={obj.shape}  dtype={obj.dtype}")
    else:
        print(f"[group]   {name}")

with h5py.File(path, "r") as f:
    f.visititems(walk)
