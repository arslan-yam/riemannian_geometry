"""Prepare a terrain for every method alike: repair folded faces (rgeo.terrain.repair_height_field) and
translate the chart to the origin (the competitors' harness is float32). Writes OFF.
usage: python3 scripts/competitors/prepare_mesh.py in.off out.off
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rgeo.terrain import read_off, write_off, repair_height_field

V, F = read_off(sys.argv[1])
V2, F2 = repair_height_field(V, F, verbose=True)
V2 = np.array(V2, dtype=float)
V2[:, :2] -= V2[:, :2].min(0)
moved = int((np.abs(np.asarray(V)[:, :2] - np.asarray(V)[:, :2].min(0) - V2[:, :2]) > 1e-9).any(1).sum())
write_off(sys.argv[2], V2, F2)
print(f"{os.path.basename(sys.argv[1])}: {len(V2)} vertices, {len(F2)} faces, {moved} vertices moved by the repair")
