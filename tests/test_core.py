"""Invariants of the core pieces. Run: python3 -m pytest tests -q  (about a minute)."""
import os, sys
import numpy as np
import pytest
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from rgeo.manifolds import PoincareDisk, GaussianHills
from rgeo.obstacles import Circle
from rgeo.graphs import build_graph, dijkstra
from rgeo.oracles import H2HIndex, SeparatorCutIndex
from rgeo.terrain import TerrainMesh, ExactTerrain, read_off

MESH = os.path.join(ROOT, "data", "terrains", "HorseMount.off")


@pytest.fixture(scope="module")
def poincare():
    m = PoincareDisk(rmax=0.92)
    pts, adj = build_graph(m, (-0.9, 0.9, -0.9, 0.9), 12, 12, [])
    rng = np.random.default_rng(0)
    qp = [(int(a), int(b)) for a, b in rng.integers(len(pts), size=(80, 2)) if a != b]
    return m, pts, adj, qp


@pytest.fixture(scope="module")
def horse():
    V, F = read_off(MESH)
    return V, F, TerrainMesh(V, F)


def test_h2h_equals_dijkstra(poincare):
    m, pts, adj, qp = poincare
    h = H2HIndex(adj)
    for s, t in qp:
        assert h.query(s, t) == pytest.approx(dijkstra(adj, s, t)[0][t], rel=1e-12)


def test_cut_index_is_an_upper_bound_and_beats_the_graph(poincare):
    #without interpolation every answer is the length of an admissible curve: >= the truth
    m, pts, adj, qp = poincare
    ix = SeparatorCutIndex(pts, adj, m, k=8, leaf_size=16, interpolate=False, rounds=1, segments=12)
    truth = np.array([PoincareDisk.exact_dist(pts[s], pts[t]) for s, t in qp])
    est = np.array([ix.query(s, t) for s, t in qp])
    graph = np.array([dijkstra(adj, s, t)[0][t] for s, t in qp])
    assert (est >= truth * (1 - 1e-9)).all()
    assert (est / truth - 1).mean() < 0.75 * (graph / truth - 1).mean() #0.75% vs 1.23% on this tiny grid


def test_exact_walk_matches_mesh_edges(horse):
    V, F, tm = horse
    E = np.array(sorted({tuple(sorted((int(t[i]), int(t[(i + 1) % 3])))) for t in F for i in range(3)}))
    L = tm.segment_lengths_batch(tm.xy[E[:, 0]], tm.xy[E[:, 1]])
    assert np.allclose(L, np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1), rtol=1e-12)


def test_walk_gradient_matches_finite_differences(horse):
    V, F, tm = horse
    rng = np.random.default_rng(0)
    A = tm.xy[F[rng.integers(len(F), size=50)]].mean(1); B = tm.xy[F[rng.integers(len(F), size=50)]].mean(1)
    L, gA, gB = tm.segment_lengths_batch(A, B, grad=True)
    ok = np.isfinite(L)
    h = 1e-3
    for c in range(2):
        E = np.zeros_like(B); E[:, c] = h
        fd = (tm.segment_lengths_batch(A, B + E) - tm.segment_lengths_batch(A, B - E)) / (2 * h)
        assert np.percentile(np.abs(fd - gB[:, c])[ok], 95) < 1e-5


def test_exact_labels_reproduce_mmp(horse):
    import pygeodesic.geodesic as pg
    V, F, tm = horse
    src = [300, 700, 400] #interior vertices; each source is a point 1 mm towards a face centroid
    X = np.array([tm.xy[v] + 1e-3 * (tm.xy[F[tm.vfaces[v][0]]].mean(0) - tm.xy[v]) / np.linalg.norm(
        tm.xy[F[tm.vfaces[v][0]]].mean(0) - tm.xy[v]) for v in src])
    D = ExactTerrain(tm).distances_from_points(X, [10, 500, 790])
    geo = pg.PyGeodesicAlgorithmExact(V, F)
    for i, s in enumerate(src):
        ref, _ = geo.geodesicDistances(np.array([s]), np.array([10, 500, 790]))
        assert np.allclose(D[i, [10, 500, 790]], ref, atol=0.01) #1 mm source shift -> < 1 cm


def test_weighted_terrain_refuses_exact_labels(horse):
    V, F, _ = horse
    tw = TerrainMesh(V, F, weights=np.full(len(F), 2.0))
    with pytest.raises(ValueError):
        ExactTerrain(tw)


def test_refinement_respects_obstacles():
    from rgeo.refine import refine_polyline
    m = GaussianHills([(3.0, 3.0, 2.2, 0.6)]); obs = [Circle(5.0, 5.0, 1.6)]
    line = np.array([[2.0, 5.0], [5.0, 7.0], [8.0, 5.0]]) #goes round the obstacle
    d, c = refine_polyline(m, line, obs, rounds=2, segments=20)
    t = np.linspace(0, 1, 50)[:, None]
    samples = np.vstack([c[i] + t * (c[i + 1] - c[i]) for i in range(len(c) - 1)])
    assert not obs[0].blocked_batch(samples).any()
