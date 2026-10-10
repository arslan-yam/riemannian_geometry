"""Multi-seed Dijkstra with early termination (C kernel `_cdijkstra.c`, scipy fallback).

`multi_dijkstra(G, seeds, targets, limit)`:
  G        scipy CSR matrix, nonnegative weights
  seeds    list (one per source) of (nodes, initial distances): a point inside a face is the set
           of that face's boundary nodes at their straight distances, so no augmented graph is built
  targets  one array shared by every source, or a list with one array per source
  limit    optional distance bound (targets farther than it come back inf)
returns    (sources x targets) array for shared targets, else a list of arrays.

A source's search stops as soon as all of its targets are settled, so the distances are those of
a full Dijkstra while the work is only what the farthest target needs.
"""
import ctypes
import os
import numpy as np
from scipy.sparse import csr_matrix, vstack as sp_vstack, hstack as sp_hstack
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from .._cbuild import load

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cdijkstra.c")
_P = ctypes.c_void_p
_I = ctypes.c_int64
_LIB = None


def _lib():
    global _LIB
    if _LIB is None:
        lib = load(_SRC)
        if lib is None:
            _LIB = False
        else:
            lib.cd_multi.argtypes = [_I, _P, _P, _P, _I, _P, _P, _P, _I, _P, _P, ctypes.c_double, _P]
            lib.cd_multi.restype = _I
            _LIB = lib
    return _LIB


def available():
    return bool(_lib())


def multi_dijkstra(G, seeds, targets, limit=np.inf, force_python=False):
    G = G.tocsr()
    n = G.shape[0]
    shared = not isinstance(targets, (list, tuple))
    ns = len(seeds)
    seed_ptr = np.zeros(ns + 1, dtype=np.int64)
    seed_ptr[1:] = np.cumsum([len(s[0]) for s in seeds])
    seed_node = np.concatenate([np.asarray(s[0], dtype=np.int64) for s in seeds]) if ns else np.zeros(0, np.int64)
    seed_dist = np.concatenate([np.asarray(s[1], dtype=np.float64) for s in seeds]) if ns else np.zeros(0)
    lib = None if force_python else _lib()
    if lib:
        indptr = np.ascontiguousarray(G.indptr, dtype=np.int64)
        indices = np.ascontiguousarray(G.indices, dtype=np.int32)
        data = np.ascontiguousarray(G.data, dtype=np.float64)
        if shared:
            T = np.ascontiguousarray(targets, dtype=np.int64)
            out = np.empty((ns, len(T)))
            tp = None
        else:
            tp = np.zeros(ns + 1, dtype=np.int64)
            tp[1:] = np.cumsum([len(t) for t in targets])
            T = np.concatenate([np.asarray(t, dtype=np.int64) for t in targets]) if ns else np.zeros(0, np.int64)
            out = np.empty(int(tp[-1]))
        rc = lib.cd_multi(n, indptr.ctypes.data, indices.ctypes.data, data.ctypes.data, ns,
                          seed_ptr.ctypes.data, seed_node.ctypes.data, seed_dist.ctypes.data,
                          len(T), None if tp is None else tp.ctypes.data, T.ctypes.data, float(limit), out.ctypes.data)
        if rc != 0:
            raise MemoryError("cd_multi ran out of memory")
        return out if shared else [out[tp[i]:tp[i + 1]] for i in range(ns)]
    #fallback: one virtual node per source joined to its seeds, scipy dijkstra in chunks
    res = []
    for c0 in range(0, ns, 32):
        c1 = min(ns, c0 + 32)
        k = c1 - c0
        r = np.repeat(np.arange(k), np.diff(seed_ptr[c0:c1 + 1]))
        sl = slice(seed_ptr[c0], seed_ptr[c1])
        Pk = csr_matrix((seed_dist[sl], (r, seed_node[sl])), shape=(k, n))
        Gaug = sp_vstack([sp_hstack([G, Pk.T]), sp_hstack([Pk, csr_matrix((k, k))])]).tocsr()
        D = sp_dijkstra(Gaug, indices=n + np.arange(k), limit=limit)
        for i in range(k):
            T = targets if shared else targets[c0 + i]
            res.append(D[i, np.asarray(T, dtype=np.int64)])
    return np.array(res).reshape(ns, -1) if shared else res
