#regenerate the benchmark meshes from the EAR-Oracle repository's datasets: repair folded faces,
#translate to the origin (their harness is float32), crop GunnisonForest to ~10k and ~1k vertices,
#and keep the full 100k GunnisonForest for the scale experiment
#
#usage: python3 scripts/make_terrains.py <path to weighted_distance_oracle clone> [--large]
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rgeo.terrain import read_off, write_off, repair_height_field

src = sys.argv[1]
large = "--large" in sys.argv
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = os.path.join(root, "data", "terrains")
os.makedirs(os.path.join(out, "large"), exist_ok=True)


def save(name, V, F, sub="."):
    V, F = repair_height_field(V, F, verbose=True)
    V = V.copy(); V[:, :2] -= V[:, :2].min(0)
    write_off(os.path.join(out, sub, f"{name}.off"), V, F)
    print(f"{name}: {len(V)} vertices, {len(F)} faces")


for name in ["HorseMount", "BigMount", "HeadlightMount", "RobinsonMountain"]:
    save(name, *read_off(os.path.join(src, "datasets", "small", f"{name}.off")))

V, F = read_off(os.path.join(src, "datasets", "median", "Simply_GunnisonForest.off"))
lo = V[:, :2].min(0); hi = V[:, :2].max(0); c = 0.5 * (lo + hi); ext = hi - lo
for frac, name in ((0.315, "Gunnison10k"), (0.10, "Gunnison1k")):
    w = frac * ext
    inside = np.all(np.abs(V[:, :2] - c) <= 0.5 * w, axis=1)
    keepF = inside[F].all(1)
    used = np.unique(F[keepF]); remap = -np.ones(len(V), dtype=int); remap[used] = np.arange(len(used))
    save(name, V[used], remap[F[keepF]])
if large:
    save("GunnisonForest100k", V, F, "large")
    V, F = read_off(os.path.join(src, "datasets", "median", "Simply_Laramie.off"))
    save("Laramie100k", V, F, "large")
