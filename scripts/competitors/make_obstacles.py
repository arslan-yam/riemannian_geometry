"""Obstacles on a prepared terrain: the faces whose chart centroid lies inside an obstacle are removed, so
the obstacle is a hole in the surface and paths go around it (the same circles and rectangles as the
smooth-metric prototype, rgeo/obstacles.py). Kinds:
  blobs  random disks, radius 2-6% of the chart side, until a fraction --cover of the chart area is gone
  walls  --walls long thin rotated rectangles (4 mean edges thick, 40-70% of the chart side), each with two
         gaps of 4% of the side: barriers that force long detours
Then the mesh is cleaned for the exact solvers: at every vertex where holes touch (several fans of faces)
only the largest fan stays, and only the largest face-connected piece is kept; unused vertices go.
usage: python3 scripts/competitors/make_obstacles.py in.off out.off --kind blobs --cover 0.1 [--seed 1]
       python3 scripts/competitors/make_obstacles.py in.off out.off --kind walls --walls 6 [--seed 1]
writes out.off and out.json (the obstacles, what was removed, the cleaning)
"""
import argparse, json, os, sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rgeo.terrain import read_off, write_off

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--kind", choices=["blobs", "walls"], required=True)
ap.add_argument("--cover", type=float, default=0.1)
ap.add_argument("--walls", type=int, default=6)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

V, F = read_off(a.inp)
V = np.asarray(V, dtype=float); F = np.asarray(F, dtype=np.int64)
xy = V[:, :2]
C = xy[F].mean(axis=1)                                     # face centroids in the chart
P = xy[F]
proj = 0.5 * np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
lo, hi = xy.min(0), xy.max(0)
L = float((hi - lo).max())
E = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
edge = float(np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1).mean())
rng = np.random.default_rng(a.seed)
removed = np.zeros(len(F), dtype=bool)
obstacles = []

if a.kind == "blobs":
    while proj[removed].sum() < a.cover * proj.sum():
        c = lo + rng.random(2) * (hi - lo)
        r = L * rng.uniform(0.02, 0.06)
        removed |= ((C - c) ** 2).sum(1) < r * r
        obstacles.append({"type": "circle", "cx": float(c[0]), "cy": float(c[1]), "r": float(r)})
else:
    thick = 4 * edge
    for _ in range(a.walls):
        c = lo + (0.1 + 0.8 * rng.random(2)) * (hi - lo)
        th = rng.uniform(0, np.pi)
        half = 0.5 * L * rng.uniform(0.4, 0.7)
        u = np.array([np.cos(th), np.sin(th)]); w = np.array([-u[1], u[0]])
        s = (C - c) @ u; t = (C - c) @ w                       # along / across the wall
        inside = (np.abs(s) < half) & (np.abs(t) < thick / 2)
        gaps = [float(g) for g in rng.uniform(-half, half, 2)]
        for g in gaps:
            inside &= ~(np.abs(s - g) < 0.02 * L)
        removed |= inside
        obstacles.append({"type": "wall", "cx": float(c[0]), "cy": float(c[1]), "angle": float(th), "half_length": float(half),
                          "thickness": float(thick), "gaps": gaps, "gap_width": float(0.04 * L)})

keep = ~removed
area_removed = float(proj[removed].sum() / proj.sum())


def edge_keys(Fk):
    e = np.concatenate([Fk[:, [0, 1]], Fk[:, [1, 2]], Fk[:, [2, 0]]])
    e.sort(axis=1)
    return e[:, 0] * len(V) + e[:, 1]


def largest_piece(Fk):
    """mask of the faces in the largest piece connected through shared edges"""
    k = edge_keys(Fk)
    fid = np.tile(np.arange(len(Fk)), 3)
    o = np.argsort(k, kind="stable"); k, fid = k[o], fid[o]
    same = np.nonzero(k[1:] == k[:-1])[0]
    A = coo_matrix((np.ones(len(same)), (fid[same], fid[same + 1])), shape=(len(Fk), len(Fk)))
    nc, lab = connected_components(A, directed=False)
    return lab == np.bincount(lab).argmax(), nc


def split_fans(Fk):
    """mask of the faces to drop: at a vertex with several fans (faces around it connected through edges
    at it), all but the largest fan"""
    drop = np.zeros(len(Fk), dtype=bool)
    vf = [[] for _ in range(len(V))]
    for f, (x, y, z) in enumerate(Fk):
        vf[x].append(f); vf[y].append(f); vf[z].append(f)
    # an interior vertex has a closed ring; only vertices on the boundary can have several fans
    nb = np.zeros(len(V), dtype=np.int64)
    k = edge_keys(Fk)
    uk, cnt = np.unique(k, return_counts=True)
    bk = uk[cnt == 1]
    np.add.at(nb, bk // len(V), 1); np.add.at(nb, bk % len(V), 1)
    nfix = 0
    for v in np.nonzero(nb > 2)[0]:                       # more than one boundary pair meets at v
        fs = vf[v]
        par = list(range(len(fs)))
        def find(i):
            while par[i] != i:
                par[i] = par[par[i]]; i = par[i]
            return i
        for i in range(len(fs)):
            for j in range(i + 1, len(fs)):
                if len(set(Fk[fs[i]]) & set(Fk[fs[j]])) == 2:  # share an edge at v
                    par[find(i)] = find(j)
        roots = [find(i) for i in range(len(fs))]
        if len(set(roots)) > 1:
            big = max(set(roots), key=roots.count)
            for i, r in enumerate(roots):
                if r != big: drop[fs[i]] = True
            nfix += 1
    return drop, nfix


Fk = F[keep]
fans_fixed = 0
pieces_dropped = 0
for _ in range(20):
    drop, nfix = split_fans(Fk)
    Fk = Fk[~drop]; fans_fixed += nfix
    m, nc = largest_piece(Fk)
    pieces_dropped += int((~m).sum())
    Fk = Fk[m]
    if nfix == 0 and nc == 1:
        break
used = np.unique(Fk)
remap = -np.ones(len(V), dtype=np.int64); remap[used] = np.arange(len(used))
V2, F2 = V[used], remap[Fk]
write_off(a.out, V2, F2)
info = {"input": a.inp, "kind": a.kind, "seed": a.seed, "cover_target": a.cover if a.kind == "blobs" else None,
        "faces_in": int(len(F)), "faces_out": int(len(F2)), "vertices_out": int(len(V2)),
        "removed_area_frac": area_removed, "removed_by_obstacles": int(removed.sum()),
        "faces_dropped_by_cleaning": int(keep.sum() - len(F2)), "fans_fixed": fans_fixed,
        "mean_edge": edge, "obstacles": obstacles}
json.dump(info, open(os.path.splitext(a.out)[0] + ".json", "w"), indent=1)
print(f"{os.path.basename(a.out)}: {len(obstacles)} {a.kind}, {area_removed:.1%} of the area removed, "
      f"{len(F2)} of {len(F)} faces kept ({fans_fixed} touching vertices split)")
