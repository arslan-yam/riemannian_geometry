#a triangulated terrain surface z = f(x, y) as a 2-d riemannian manifold in the (x, y) chart.
#
#this is exactly the object the terrain distance oracle literature works on: a piecewise-linear
#height field over a triangulation, whose geodesics are the MMP / CGAL "exact" ones. the induced
#metric is g = I + grad f grad f^T with grad f constant on every face, so lengths of chart-straight
#segments are computed EXACTLY by walking the faces the segment crosses and summing the lifted
#3-d pieces. quadrature would be wrong by O(face jump) wherever a segment crosses an edge, which
#is the same order as the errors being measured, so it is not used anywhere on this manifold
import numpy as np
from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree
from ..manifolds import Manifold, build_metric


class TerrainMesh(Manifold):
    #weights: optional per-face cost multipliers (the weighted-region terrains of the oracle
    #papers: path cost = sum over faces of weight * 3-d length inside the face, Snell refraction
    #at edges). the metric becomes w_f^2 (I + grad f grad f^T); the exact face walk and its
    #gradient carry the weight per piece, so refinement works unchanged. no exact solver exists
    #for weighted terrains, so ExactTerrain refuses them and labels must be refined
    def __init__(self, V, F, shift=True, weights=None):
        V = np.asarray(V, dtype=float)
        F = np.asarray(F, dtype=int)
        self.origin = V[:, :2].min(0) if shift else np.zeros(2) #keep coordinates O(1e4), not 1e7
        self.V = V.copy()
        self.V[:, :2] -= self.origin
        self.F = F
        self.tri = Triangulation(self.V[:, 0], self.V[:, 1], F)
        self.finder = self.tri.get_trifinder()
        self.neighbors = self.tri.neighbors #[f, k] shares the edge from F[f,k] to F[f,(k+1)%3]

        #plane z = a x + b y + c on every face, so grad f = (a, b) is constant per face
        P = self.V[F] #(m, 3, 3)
        x, y, z = P[..., 0], P[..., 1], P[..., 2]
        det = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        self.a = ((z[:, 1] - z[:, 0]) * (y[:, 2] - y[:, 0]) - (z[:, 2] - z[:, 0]) * (y[:, 1] - y[:, 0])) / det
        self.b = ((x[:, 1] - x[:, 0]) * (z[:, 2] - z[:, 0]) - (x[:, 2] - x[:, 0]) * (z[:, 1] - z[:, 0])) / det
        self.c = z[:, 0] - self.a * x[:, 0] - self.b * y[:, 0]
        self.w = np.ones(len(F)) if weights is None else np.asarray(weights, dtype=float)
        self.weighted = weights is not None and not np.allclose(self.w, 1.0)
        self.lambda_min = float(self.w.min()) #|v|_g = w |Mv| >= w_min |v|
        self.xy = self.V[:, :2]
        self.bbox = (self.xy[:, 0].min(), self.xy[:, 0].max(), self.xy[:, 1].min(), self.xy[:, 1].max())
        self.vtree = cKDTree(self.xy)
        self.ev0 = self.xy[F] #(m, 3, 2) edge start points, edge k runs from vertex k to k+1
        self.ee = np.roll(self.ev0, -1, axis=1) - self.ev0
        self.vfaces = [[] for _ in range(len(V))] #faces around each vertex, for tolerant location
        for fi, tri in enumerate(F):
            for v in tri:
                self.vfaces[v].append(fi)

    # ------------------------------------------------------------------ chart helpers
    def faces_of(self, pts, tol=1e-7):
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        f = np.asarray(self.finder(pts[:, 0], pts[:, 1]))
        miss = np.nonzero(f < 0)[0]
        if len(miss) == 0:
            return f
        #the trapezoid-map finder has no tolerance: a point on the boundary edge (which is
        #where every shrunk boundary-edge segment starts) can fall a rounding error outside.
        #retry among the faces around the nearest vertex with a barycentric tolerance
        _, nv = self.vtree.query(pts[miss])
        for q, v in zip(miss, nv):
            best, bestw = -1, -np.inf
            for fc in self.vfaces[v]:
                tri = self.F[fc]
                x0, y0 = self.xy[tri[0]]; x1, y1 = self.xy[tri[1]]; x2, y2 = self.xy[tri[2]]
                det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
                if det == 0:
                    continue
                l1 = ((pts[q, 0] - x0) * (y2 - y0) - (x2 - x0) * (pts[q, 1] - y0)) / det
                l2 = ((x1 - x0) * (pts[q, 1] - y0) - (pts[q, 0] - x0) * (y1 - y0)) / det
                w = min(l1, l2, 1.0 - l1 - l2)
                if w > bestw:
                    best, bestw = fc, w
            if bestw >= -tol:
                f[q] = best
        return f

    def height(self, pts):
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        f = self.faces_of(pts)
        ok = f >= 0
        z = np.zeros(len(pts))
        z[ok] = self.a[f[ok]] * pts[ok, 0] + self.b[f[ok]] * pts[ok, 1] + self.c[f[ok]]
        return z

    def lift(self, pts):
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        return np.column_stack([pts, self.height(pts)])

    def in_domain(self, p):
        return bool(self.finder(float(p[0]), float(p[1])) >= 0)

    def in_domain_batch(self, points):
        return self.faces_of(points) >= 0

    # ------------------------------------------------------------------ metric (piecewise constant)
    def metric(self, p):
        f = int(self.finder(float(p[0]), float(p[1])))
        if f < 0:
            return np.eye(2)
        fx, fy = self.a[f], self.b[f]
        return self.w[f] ** 2 * np.array([[1.0 + fx * fx, fx * fy], [fx * fy, 1.0 + fy * fy]])

    def metric_batch(self, points):
        f = self.faces_of(points)
        fx = np.where(f >= 0, self.a[np.maximum(f, 0)], 0.0)
        fy = np.where(f >= 0, self.b[np.maximum(f, 0)], 0.0)
        w = np.where(f >= 0, self.w[np.maximum(f, 0)], 1.0)
        return build_metric(fx, fy) * (w * w)[:, None, None]

    def metric_grad_batch(self, points, h=None):
        z = np.zeros((len(points), 2, 2)) #constant on every face
        return z, z.copy()

    # ------------------------------------------------------------------ exact polyhedral lengths
    def segment_lengths_batch(self, A, B, grad=False, max_hops=4096, shrink=1e-7):
        #exact lifted length of every chart-straight segment A[i] -> B[i]. all segments walk the
        #triangulation in lock-step: one vectorised hop per iteration. with grad=True also returns
        #dL/dA and dL/dB. L = sum_k (t_{k+1} - t_k) s_k with s_k = |M_k d| the lifted speed on
        #face k and t_k the edge-crossing parameters, so the gradient has two parts: the speeds
        #change with d = B - A, and the crossing parameters t_k move when the segment rotates.
        #
        #degeneracies are the norm, not the exception, on a mesh graph: every graph edge IS a
        #mesh edge, so segments start at vertices and run along edges. the walk therefore runs
        #on the segment shrunk by a hair at both ends (A' = A + eps d, B' = B - eps d), which
        #takes the endpoints off the vertices, and adds the two cut-off pieces back exactly
        #(speed is constant inside a face). a segment that still gets stuck (passes exactly
        #through a vertex) is relocated by point location a hair further along
        A = np.atleast_2d(np.asarray(A, dtype=float)); B = np.atleast_2d(np.asarray(B, dtype=float))
        N = len(A)
        L = np.zeros(N)
        gA = np.zeros((N, 2)); gB = np.zeros((N, 2))
        d0 = B - A
        A = A + shrink * d0
        B = B - shrink * d0
        fa, fb = self.faces_of(A), self.faces_of(B)
        ok = (fa >= 0) & (fb >= 0)
        L[~ok] = np.inf
        if not ok.any():
            return (L, gA, gB) if grad else L

        d = B - A
        face = fa.copy()
        t = np.zeros(N)
        active = ok.copy()
        s_first = np.zeros(N); s_last = np.zeros(N)
        first = np.ones(N, dtype=bool)
        for _ in range(max_hops):
            idx = np.nonzero(active)[0]
            if len(idx) == 0:
                break
            f = face[idx]
            #the three edges k: from vertex k to vertex (k+1)%3, matching self.neighbors
            v0 = self.ev0[f] #(n, 3, 2)
            e = self.ee[f]
            di = d[idx]
            dd = di[:, None, :] #(n, 1, 2)
            aa = A[idx][:, None, :]
            #solve aa + tt*dd = v0 + ss*e  for tt, ss  (2x2 cramer)
            den = dd[..., 0] * (-e[..., 1]) - dd[..., 1] * (-e[..., 0])
            rhs = v0 - aa
            with np.errstate(divide="ignore", invalid="ignore"):
                tt = (rhs[..., 0] * (-e[..., 1]) - rhs[..., 1] * (-e[..., 0])) / den
                ss = (dd[..., 0] * rhs[..., 1] - dd[..., 1] * rhs[..., 0]) / den
            tcur = t[idx][:, None]
            #a segment running along an edge is not a crossing of it: 0/0 would give noise
            transversal = np.abs(den) > 1e-10 * np.hypot(dd[..., 0], dd[..., 1]) * np.hypot(e[..., 0], e[..., 1])
            #a piece that lies ON an edge costs the cheaper of the two faces (a path hugging the
            #edge on the cheap side), which is also what the mesh graph charges for that edge
            along = ~transversal & (np.abs(rhs[..., 0] * e[..., 1] - rhs[..., 1] * e[..., 0])
                                    <= 1e-9 * np.hypot(e[..., 0], e[..., 1]) * (1.0 + np.hypot(rhs[..., 0], rhs[..., 1])))
            valid = transversal & np.isfinite(tt) & (ss >= -1e-9) & (ss <= 1 + 1e-9) & (tt > tcur + 1e-12) & (tt <= 1.0 + 1e-12)
            tt = np.where(valid, tt, np.inf)
            k = tt.argmin(1)
            ar = np.arange(len(idx))
            texit = tt[ar, k]
            at_end = face[idx] == fb[idx]
            stuck = ~np.isfinite(texit) & ~at_end
            #stuck: no edge crossing ahead but not in the target face. relocate a hair ahead:
            #outside the mesh means the segment leaves it (a start exactly at a boundary vertex
            #heading out), the same face means finish here, another face means hop there.
            #the probe offset must exceed the location tolerance of faces_of, or it lands back
            #in the face it just left
            out_now = np.zeros(len(idx), dtype=bool)
            if stuck.any():
                probe = A[idx[stuck]] + np.minimum(t[idx[stuck]] + 1e-5, 1.0)[:, None] * di[stuck]
                nf = self.faces_of(probe)
                moved = (nf >= 0) & (nf != face[idx[stuck]])
                sidx = idx[stuck]
                face[sidx[moved]] = nf[moved]
                out_now[stuck] = nf < 0
                stuck[stuck] = moved
            ends_here = (~np.isfinite(texit) | at_end & (texit >= 1.0 - 1e-12)) & ~stuck & ~out_now
            texit = np.where(ends_here | (texit > 1.0), 1.0, texit)
            texit = np.where(stuck, t[idx], texit) #relocated: no progress this round

            #lifted piece from t to texit inside this face: length (texit - t) * s, s the weighted speed
            af, bf, wf = self.a[f], self.b[f], self.w[f].copy()
            if self.weighted and along.any():
                ka = along.argmax(1); onedge = along.any(1)
                nb = self.neighbors[f[onedge], ka[onedge]]
                wf[onedge] = np.minimum(wf[onedge], np.where(nb >= 0, self.w[np.maximum(nb, 0)], np.inf))
            s = wf * np.sqrt(di[:, 0] ** 2 + di[:, 1] ** 2 + (af * di[:, 0] + bf * di[:, 1]) ** 2)
            dt = texit - t[idx]
            L[idx] += dt * s
            fi = first[idx] & ~stuck
            s_first[idx[fi]] = s[fi]; first[idx[fi]] = False
            s_last[idx[~stuck]] = s[~stuck]
            done = (ends_here | (texit >= 1.0 - 1e-12)) & ~stuck
            nxt = self.neighbors[f, k].copy()
            left = ((~done) & (nxt < 0) & ~stuck) | out_now
            if (left & ~out_now).any():
                #a boundary edge was hit. at a boundary vertex the crossing can be assigned to
                #the wrong one of the two edges meeting there, so before declaring the segment
                #outside, ask point location a hair beyond the crossing
                chk = left & ~out_now
                probe = A[idx[chk]] + np.minimum(texit[chk] + 1e-5, 1.0)[:, None] * di[chk]
                nf = self.faces_of(probe)
                back = nf >= 0
                lidx = np.nonzero(chk)[0]
                nxt[lidx[back]] = nf[back]
                left[lidx[back]] = False
            if grad:
                #speed part: d s / d d = M^T u with u the unit lifted direction (constant on the piece)
                with np.errstate(divide="ignore", invalid="ignore"):
                    uz = (af * di[:, 0] + bf * di[:, 1]) / s
                    ux, uy = di[:, 0] / s, di[:, 1] / s
                gs = np.column_stack([ux + af * uz, uy + bf * uz]) * (wf * wf)[:, None] #d(w|Md|)/dd = w^2 M^T M d / s
                gs = np.where(np.isfinite(gs), gs, 0.0)
                gB[idx] += dt[:, None] * gs
                gA[idx] -= dt[:, None] * gs
                #crossing part: t_k = cross(v0 - A, e) / cross(d, e) moves when d moves
                cr = ~done & ~left & ~stuck
                if cr.any():
                    j = idx[cr]; nf = nxt[cr]; dj = d[j]
                    s_next = self.w[nf] * np.sqrt(dj[:, 0] ** 2 + dj[:, 1] ** 2 + (self.a[nf] * dj[:, 0] + self.b[nf] * dj[:, 1]) ** 2)
                    ek = e[ar[cr], k[cr]] #(c, 2)
                    cross = dj[:, 0] * ek[:, 1] - dj[:, 1] * ek[:, 0]
                    w = np.column_stack([ek[:, 1], -ek[:, 0]]) / cross[:, None]
                    ds = (s[cr] - s_next)[:, None]
                    gB[j] += ds * (-texit[cr][:, None]) * w
                    gA[j] += ds * (texit[cr] - 1.0)[:, None] * w
            #hop; a -1 neighbour means the segment leaves the mesh: mark infinite
            L[idx[left]] = np.inf
            active[idx[done | left]] = False
            mv = ~(done | left | stuck)
            face[idx[mv]] = nxt[mv]
            t[idx[mv]] = texit[mv]
        else:
            L[active] = np.inf
        #add back the two shrunk-off pieces: eps |d0| worth of the first and last face speeds
        tail = shrink / (1.0 - 2.0 * shrink)
        L[ok] += tail * (s_first[ok] + s_last[ok])
        return (L, gA, gB) if grad else L

    def segment_lenght(self, p, q, n=None):
        return float(self.segment_lengths_batch(np.asarray(p)[None], np.asarray(q)[None])[0])

    # ------------------------------------------------------------------ energy on the polyhedral surface
    # E = m sum_i L_i^2 with L_i the EXACT lifted length of the chart-straight segment i, so the
    # energy is continuous in the vertices (chart quadrature would jump whenever a sample point
    # crosses a face edge) and its gradient is 2 m L_i dL_i/dP at the two ends of every segment
    def _penalise_outside(self, A, B, L, gA, gB):
        #a segment that leaves the surface has no length. give it a large finite one that
        #shrinks as the segment shrinks, so the optimiser is pushed back instead of stalling
        #on an infinite objective (which stops L-BFGS dead for the whole stack)
        bad = ~np.isfinite(L)
        if bad.any():
            d = B[bad] - A[bad]
            nd = np.linalg.norm(d, axis=1)
            L = L.copy(); gA = gA.copy(); gB = gB.copy()
            L[bad] = 100.0 * nd
            u = 100.0 * d / np.maximum(nd, 1e-12)[:, None]
            gA[bad] = -u; gB[bad] = u
        return L, gA, gB

    def energy_and_grad_points(self, points, m):
        P = np.asarray(points, dtype=float)
        L, gA, gB = self.segment_lengths_batch(P[:-1], P[1:], grad=True)
        L, gA, gB = self._penalise_outside(P[:-1], P[1:], L, gA, gB)
        e = m * float((L * L).sum())
        grad = np.zeros_like(P)
        grad[:-1] += 2.0 * L[:, None] * gA
        grad[1:] += 2.0 * L[:, None] * gB
        return e, m * grad[1:-1].ravel()

    def energy_and_grad_stack(self, P, m):
        B = P.shape[0]
        A0, B0 = P[:, :-1].reshape(-1, 2), P[:, 1:].reshape(-1, 2)
        L, gA, gB = self.segment_lengths_batch(A0, B0, grad=True)
        L, gA, gB = self._penalise_outside(A0, B0, L, gA, gB)
        e = m * float((L * L).sum())
        Lr = L.reshape(B, m, 1)
        grad = np.zeros_like(P)
        grad[:, :-1] += 2.0 * Lr * gA.reshape(B, m, 2)
        grad[:, 1:] += 2.0 * Lr * gB.reshape(B, m, 2)
        return e, m * grad[:, 1:-1]

    # ------------------------------------------------------------------ cut line vs mesh edges
    def edge_crossings(self, axis, sv, lo, hi):
        #cut coordinates where the line {x[axis] = sv} crosses mesh edges, within [lo, hi].
        #on a weighted terrain the cheap corridors are the edges, and a weighted geodesic
        #crosses a cut where a cheap edge does; uniform portals miss those narrow dips
        if not hasattr(self, "_edges"):
            E = set()
            for tri in self.F:
                for i in range(3):
                    u, v = int(tri[i]), int(tri[(i + 1) % 3])
                    E.add((min(u, v), max(u, v)))
            self._edges = np.array(sorted(E), dtype=int)
        A = self.xy[self._edges[:, 0]]; B = self.xy[self._edges[:, 1]]
        da, db = A[:, axis] - sv, B[:, axis] - sv
        m = da * db < 0
        t = da[m] / (da[m] - db[m])
        c = A[m, 1 - axis] + t * (B[m, 1 - axis] - A[m, 1 - axis])
        c = c[(c >= lo) & (c <= hi)]
        return np.sort(c)

    # ------------------------------------------------------------------ graph on the mesh
    def mesh_graph(self):
        #vertices are the graph nodes, mesh edges the graph edges, weights the exact 3-d edge
        #lengths. this is the graph every terrain oracle in the literature starts from
        n = len(self.V)
        edges = set()
        for t in self.F:
            for i in range(3):
                u, v = int(t[i]), int(t[(i + 1) % 3])
                edges.add((min(u, v), max(u, v)))
        edge_w = {}
        for fi, t in enumerate(self.F):
            for i in range(3):
                e = (min(int(t[i]), int(t[(i + 1) % 3])), max(int(t[i]), int(t[(i + 1) % 3])))
                edge_w[e] = min(edge_w.get(e, np.inf), float(self.w[fi])) #walk the edge on its cheaper side
        adj = [dict() for _ in range(n)]
        for u, v in edges:
            w = edge_w[(u, v)] * float(np.linalg.norm(self.V[u] - self.V[v]))
            adj[u][v] = w; adj[v][u] = w
        return [tuple(p) for p in self.xy], adj
