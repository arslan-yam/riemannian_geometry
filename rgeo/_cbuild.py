"""Build-on-first-use for the two small C kernels (compiled query, multi-seed Dijkstra).

`load(src)` compiles `src` with the system compiler (`$CC`, default `cc`) into a shared library
named after the source's hash, next to the source (or in the temp dir if that is read-only), and
returns the ctypes handle, or None when no compiler works. Callers keep a pure-Python fallback.
"""
import ctypes
import hashlib
import os
import subprocess
import tempfile

_CACHE = {}


def load(src):
    if src in _CACHE:
        return _CACHE[src]
    lib = None
    code = open(src, "rb").read()
    tag = hashlib.sha1(code).hexdigest()[:12]
    base = os.path.splitext(os.path.basename(src))[0]
    cc = os.environ.get("CC", "cc")
    for d in (os.path.dirname(src), tempfile.gettempdir()):
        so = os.path.join(d, f"{base}_{tag}.so")
        if not os.path.exists(so):
            try:
                subprocess.run([cc, "-O2", "-ffp-contract=off", "-shared", "-fPIC", "-o", so, src],
                               check=True, capture_output=True)
            except (OSError, subprocess.CalledProcessError):
                continue
        try:
            lib = ctypes.CDLL(so)
            break
        except OSError:
            continue
    _CACHE[src] = lib
    return lib
