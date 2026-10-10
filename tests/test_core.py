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


# ---------------------------------------------------------------------- exact obstacle test
def test_exact_segment_test_contains_dense_sampling():
    from rgeo.obstacles import Rectangle, segments_blocked
    rng = np.random.default_rng(1)
    A = rng.uniform(0, 10, (3000, 2)); B = rng.uniform(0, 10, (3000, 2))
    t = np.linspace(0, 1, 4001)
    for ob in (Circle(5.0, 5.0, 1.7), Rectangle(3.0, 6.5, 2.0, 4.0)):
        exact = segments_blocked(A, B, [ob])
        dense = np.array([ob.blocked_batch(a + t[:, None] * (b - a)).any() for a, b in zip(A, B)])
        assert (exact | ~dense).all() #every sampled hit is an exact hit
        assert (exact != dense).mean() < 0.002 #the rest are grazing segments the samples step over


# ---------------------------------------------------------------------- sub-meshes
def test_face_grid_matches_brute_force(horse):
    from rgeo.terrain.submesh import FaceGrid
    V, F, tm = horse
    grid = FaceGrid(tm)
    P = tm.xy[F]; flo, fhi = P.min(1), P.max(1)
    rng = np.random.default_rng(2)
    lo, hi = tm.xy.min(0), tm.xy.max(0)
    for _ in range(50):
        c = rng.uniform(lo, hi); h = rng.uniform(0, 0.3) * (hi - lo)
        x0, x1, y0, y1 = c[0] - h[0], c[0] + h[0], c[1] - h[1], c[1] + h[1]
        ref = np.nonzero((flo[:, 0] <= x1) & (fhi[:, 0] >= x0) & (flo[:, 1] <= y1) & (fhi[:, 1] >= y0))[0]
        assert np.array_equal(np.sort(grid.faces_in_box(x0, x1, y0, y1)), ref)


def test_inserted_points_leave_the_surface_unchanged(horse):
    import pygeodesic.geodesic as pg
    from rgeo.terrain.submesh import extract, insert_points
    V, F, tm = horse
    faces = np.arange(len(F))
    rng = np.random.default_rng(3)
    f0 = int(rng.integers(len(F))) #40 points in one face (delaunay path) + 60 in random faces (1->3 splits)
    X = np.vstack([rng.dirichlet([1, 1, 1], 40) @ tm.xy[F[f0]],
                   np.array([rng.dirichlet([1, 1, 1]) @ tm.xy[F[f]] for f in rng.integers(len(F), size=60)])])
    V0, F0, used, remap = extract(tm, faces)
    V2, F2, idx, Xu = insert_points(tm, V0, F0, faces, X)
    area = lambda V, F: 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum()
    assert area(V2, F2) == pytest.approx(area(V0, F0), rel=1e-12)
    tg = np.array([5, 200, 650], dtype=np.int32)
    d0, _ = pg.PyGeodesicAlgorithmExact(V0, F0).geodesicDistances(np.array([400], np.int32), tg)
    d2, _ = pg.PyGeodesicAlgorithmExact(V2, F2.astype(np.int32)).geodesicDistances(np.array([400], np.int32), tg)
    assert np.allclose(d0, d2, rtol=1e-9)


# ---------------------------------------------------------------------- certified index
@pytest.fixture(scope="module")
def certified(horse):
    import pygeodesic.geodesic as pg
    from rgeo.oracles.certified import CertifiedCutIndex
    V, F, tm = horse
    ix = CertifiedCutIndex(tm, eps=0.2, gap_factor=1.0)
    rng = np.random.default_rng(4)
    qp = [(int(a), int(b)) for a, b in rng.integers(len(V), size=(160, 2)) if a != b]
    geo = pg.PyGeodesicAlgorithmExact(V, F)
    truth = np.array([geo.geodesicDistance(s, t)[0] for s, t in qp])
    return ix, qp, truth


def test_certified_bounds_hold_and_are_within_eps(certified):
    ix, qp, truth = certified
    B = np.array([ix.bounds(s, t) for s, t in qp])
    lb, ub, est = B.T
    assert (lb <= truth).all() and (ub >= truth).all()
    assert (ub <= (1 + ix.eps) * truth).all() and (lb >= (1 - ix.eps) * truth).all()
    assert ((lb <= est) & (est <= ub)).all()
    assert np.abs(est / truth - 1).mean() < 1e-3 #measured: ~0.07% at eps = 0.2, g = 1 edge


def test_compiled_query_matches_python(certified):
    from rgeo.oracles import cquery
    if not cquery.available():
        pytest.skip("no C compiler")
    ix, qp, _ = certified
    cq = cquery.CompiledQuery(ix)
    P = np.array([ix.bounds(s, t) for s, t in qp])
    C = np.array([cq.bounds(s, t) for s, t in qp])
    assert np.allclose(C, P, rtol=1e-12, atol=0)
    assert np.array_equal(cq.bounds_batch([s for s, _ in qp], [t for _, t in qp]), C)


def test_steiner_labels_are_admissible(horse):
    #on an unweighted terrain the steiner-graph labels must stay above the exact distances
    import pygeodesic.geodesic as pg
    from rgeo.oracles.certified import CertifiedCutIndex
    from rgeo.terrain.submesh import FaceGrid
    from rgeo.terrain.steiner import SteinerLabels
    V, F, tm = horse
    ix = CertifiedCutIndex(tm, eps=0.2, gap_factor=1.0, oracle=SteinerLabels(tm, FaceGrid(tm), m=2))
    rng = np.random.default_rng(5)
    qp = [(int(a), int(b)) for a, b in rng.integers(len(V), size=(60, 2)) if a != b]
    geo = pg.PyGeodesicAlgorithmExact(V, F)
    truth = np.array([geo.geodesicDistance(s, t)[0] for s, t in qp])
    B = np.array([ix.bounds(s, t) for s, t in qp])
    assert (B[:, 0] <= truth).all() and (B[:, 1] >= truth * (1 - 1e-6)).all()
    assert (B[:, 1] / truth - 1).mean() < 0.02
