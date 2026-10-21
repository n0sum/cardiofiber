#!/usr/bin/env python3
"""
Mesh AugmentA (mm) -> mesh openCARP (um).

Solo le coordinate del .pts si moltiplicano per 1000; .elem e .lon si
copiano invariati.
"""

import shutil
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- percorsi
ROOT    = Path(__file__).resolve().parents[3]                    # radice del repository
SRC     = ROOT / "data" / "augmenta"  / "LA_endo_with_fiber"     # senza estensione
DST     = ROOT / "data" / "reference" / "LA_um"
# --------------------------------------------------------------------------


def main():
    pts = np.loadtxt(f"{SRC}.pts", skiprows=1, dtype=np.float64) * 1000.0
    DST.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(f"{DST}.pts", pts, fmt="%.4f", header=str(len(pts)), comments="")
    shutil.copyfile(f"{SRC}.elem", f"{DST}.elem")
    shutil.copyfile(f"{SRC}.lon", f"{DST}.lon")
    print(f"scritto {DST}.{{pts,elem,lon}}   punti {len(pts)}")
    print(f"  estensione (mm): {np.round((pts.max(0) - pts.min(0)) / 1000, 2)}")


if __name__ == "__main__":
    main()
