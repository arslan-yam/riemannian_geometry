"""Compiled query for CertifiedCutIndex (C via ctypes, built on first use with the system compiler).

`CompiledQuery(index)` packs the index into flat arrays once; `bounds(s, t)` then runs the same
arithmetic as `CertifiedCutIndex.bounds` in C (about 1 us instead of ~30 us of numpy calls), and
`bounds_batch(S, T)` answers many pairs in one call. Results agree with the Python query to
rounding (compiled with -ffp-contract=off; tests/test_core.py checks it). If no compiler is
available, `available()` is False and the Python query remains the reference implementation.
"""
import ctypes
import os
import numpy as np
from .._cbuild import load

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cquery.c")
_LIB = None
_P = ctypes.c_void_p
_I = ctypes.c_int64
_D = ctypes.c_double


class _CIndex(ctypes.Structure):
    _fields_ = [("n", _I), ("depth", _I), ("exact", _I), ("interp", _I),
                ("Rn2", _D), ("lam", _D), ("f32rel", _D),
                ("anc", _P), ("row", _P), ("depth_of", _P), ("xy", _P), ("V3", _P),
                ("near_ptr", _P), ("near_ids", _P), ("near_d", _P),
                ("L_off", _P), ("P_off", _P), ("M_off", _P), ("K", _P), ("L", _P),
                ("gap", _P), ("arc", _P), ("run", _P),
                ("foot_r", _P), ("foot_arc", _P), ("foot_run", _P),
                ("E_off", _P), ("E_n", _P), ("end_idx", _P), ("end_gap", _P), ("lb_ok", _P)]


def _lib():
    global _LIB
    if _LIB is None:
        lib = load(_SRC)
        if lib is None:
            _LIB = False
        else:
            lib.cq_bounds.argtypes = [_P, _I, _I, _P, _P]
            lib.cq_bounds.restype = None
            lib.cq_bounds_batch.argtypes = [_P, _I, _P, _P, _P, _P]
            lib.cq_bounds_batch.restype = None
            _LIB = lib
    return _LIB


def available():
    return bool(_lib())


class CompiledQuery:
    def __init__(self, ix):
        lib = _lib()
        if not lib:
            raise OSError("no C compiler found for the compiled query")
        self._lib, self._keep = lib, []
        nc = len(ix.cut_members)
        K = np.array([L.shape[1] for L in ix.L], dtype=np.int64)
        Mn = np.array([len(M) for M in ix.cut_members], dtype=np.int64)
        L_off = np.r_[0, np.cumsum(K * Mn)[:-1]].astype(np.int64)
        P_off = np.r_[0, np.cumsum(K)[:-1]].astype(np.int64)
        M_off = np.r_[0, np.cumsum(Mn)[:-1]].astype(np.int64)
        cat = lambda xs, dt: np.concatenate([np.asarray(x, dtype=dt).ravel() for x in xs]) if nc else np.zeros(0, dt)
        #gap padded to one entry per portal (the last one of a cut is never read)
        gap = cat([np.r_[g, np.inf] if k else np.zeros(0) for g, k in zip(ix.cut_gap, K)], np.float64)
        E_n = np.array([len(e[0]) for e in ix.cut_end], dtype=np.int64)
        E_off = np.r_[0, np.cumsum(E_n)[:-1]].astype(np.int64)
        interp = {False: 0, None: 0, "parabola": 1, True: 1, "model": 2}.get(ix.interpolate, 1)
        arrays = dict(anc=(ix.anc, np.int32), row=(ix.row, np.int32), depth_of=(ix.depth_of, np.int32),
                      xy=(ix.xy, np.float64), V3=(ix.V3, np.float64),
                      near_ptr=(ix.near_ptr, np.int64), near_ids=(ix.near_ids, np.int32), near_d=(ix.near_d, np.float32),
                      L_off=(L_off, np.int64), P_off=(P_off, np.int64), M_off=(M_off, np.int64), K=(K, np.int64),
                      L=(cat(ix.L, np.float32), np.float32), gap=(gap, np.float64),
                      arc=(cat(ix.cut_arc, np.float64), np.float64), run=(cat(ix.cut_run, np.int64), np.int64),
                      foot_r=(cat(ix.foot_r, np.float32), np.float32), foot_arc=(cat(ix.foot_arc, np.float32), np.float32),
                      foot_run=(cat(ix.foot_run, np.int64), np.int64),
                      E_off=(E_off, np.int64), E_n=(E_n, np.int64),
                      end_idx=(cat([e[0] for e in ix.cut_end], np.int64), np.int64),
                      end_gap=(cat([e[1] for e in ix.cut_end], np.float64), np.float64),
                      lb_ok=(ix.cut_lb_ok, np.uint8))
        st = _CIndex(n=ix.n, depth=ix.anc.shape[1], exact=int(bool(ix.oracle.exact)), interp=interp,
                     Rn2=ix.Rq * ix.Rq, lam=ix.lam, f32rel=2.0 ** -22)
        for name, (a, dt) in arrays.items():
            a = np.ascontiguousarray(a, dtype=dt)
            self._keep.append(a)
            setattr(st, name, a.ctypes.data)
        self._st = st
        self._ref = ctypes.addressof(st)
        self._out = np.zeros(3)
        self._f = np.zeros(max(1, int(K.max()) if nc else 1))
        self._po, self._pf = self._out.ctypes.data, self._f.ctypes.data

    def bounds(self, s, t):
        self._lib.cq_bounds(self._ref, s, t, self._po, self._pf)
        o = self._out
        return float(o[0]), float(o[1]), float(o[2])

    def query(self, s, t):
        self._lib.cq_bounds(self._ref, s, t, self._po, self._pf)
        return float(self._out[2])

    def query_ub(self, s, t):
        self._lib.cq_bounds(self._ref, s, t, self._po, self._pf)
        return float(self._out[1])

    def bounds_batch(self, S, T):
        S = np.ascontiguousarray(S, dtype=np.int64); T = np.ascontiguousarray(T, dtype=np.int64)
        out = np.zeros((len(S), 3))
        self._lib.cq_bounds_batch(self._ref, len(S), S.ctypes.data, T.ctypes.data, out.ctypes.data, self._pf)
        return out
