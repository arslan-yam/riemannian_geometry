"""Face weights for a weighted terrain: uniform in (1, 5), the competitors' own distribution (as in
bench/terrain_bench.py --weighted), one per line in face order (the harness's face_weight.query format).
usage: python3 scripts/competitors/make_weights.py mesh.off seed out.w
"""
import sys
import numpy as np

with open(sys.argv[1]) as f:
    assert f.readline().strip() == "OFF"
    nv, nf, _ = map(int, f.readline().split())
w = np.random.default_rng(int(sys.argv[2])).uniform(1.000001, 5.0, nf)
np.savetxt(sys.argv[3], w, fmt="%.6f")
print(f"{nf} face weights -> {sys.argv[3]}")
