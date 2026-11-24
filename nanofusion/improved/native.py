"""Load the C force kernel, compiling it on first use; fall back to numpy.

The kernel is compiled with the system C compiler (``cc``) into a per-user
cache directory, keyed by a hash of the source, with OpenMP when available.
Without a compiler the pure-numpy implementation below is used. It computes
the same quantities, only more slowly.
"""

from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from .physics import J_PER_KEV, K_COULOMB, MB_TO_M2, sigma_dd_mb, sigma_dli_mb

SOURCE = Path(__file__).with_name("kernel.c")
NBINS = 64  # pair-exposure histogram bins, 8 per decade from 1e-4 keV (must match kernel.c)
HIST_EDGES_KEV = 10.0 ** (np.arange(NBINS + 1) / 8.0 - 4.0)

StepFn = Callable[..., np.ndarray]


def _cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "nanofusion"


def _compile() -> Optional[Path]:
    compiler = shutil.which(os.environ.get("CC", "cc")) or shutil.which("gcc") or shutil.which("clang")
    if compiler is None:
        return None
    source = SOURCE.read_bytes()
    digest = hashlib.sha256(source + sys.platform.encode()).hexdigest()[:16]
    target = _cache_dir() / f"kernel-{digest}.so"
    if target.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    base = [compiler, "-O3", "-shared", "-fPIC", str(SOURCE), "-lm"]
    for extra in (["-fopenmp"], []):  # try with OpenMP first
        with tempfile.NamedTemporaryFile(suffix=".so", dir=target.parent, delete=False) as tmp:
            tmp_path = Path(tmp.name)
        result = subprocess.run(base[:4] + extra + ["-o", str(tmp_path)] + base[4:], capture_output=True)
        if result.returncode == 0:
            os.replace(tmp_path, target)
            return target
        tmp_path.unlink(missing_ok=True)
    return None


_library: Optional[ctypes.CDLL] = None
_tried = False


def native_available() -> bool:
    return _load() is not None


def _load() -> Optional[ctypes.CDLL]:
    global _library, _tried
    if _tried:
        return _library
    _tried = True
    if os.environ.get("NANOFUSION_NO_NATIVE"):
        return None
    try:
        path = _compile()
        if path is not None:
            lib = ctypes.CDLL(str(path))
            dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
            ip = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
            lib.coulomb_step.argtypes = [ctypes.c_int, dp, dp, dp, dp, ip, ctypes.c_double, ctypes.c_double,
                                         ctypes.c_double, ctypes.c_double, ctypes.c_double, dp, dp, dp, dp]
            lib.coulomb_step.restype = None
            _library = lib
    except OSError:
        _library = None
    return _library


def coulomb_step(x, v, q, mass, kind, soft2, rc2, inv_vc, s_dli, bg_dli):
    """Forces (n, 3), stats [rate_DLi, rate_DD, min_dLi, maxE_DLi, min_tau, potential, min_r],
    and pair-exposure histograms (D+D, D+Li) over relative kinetic energy."""
    lib = _load()
    n = len(q)
    if lib is not None:
        force = np.empty((n, 3))
        stats = np.empty(7)
        hist_dd = np.empty(NBINS)
        hist_dli = np.empty(NBINS)
        lib.coulomb_step(n, np.ascontiguousarray(x), np.ascontiguousarray(v), q, mass, kind,
                         soft2, rc2, inv_vc, s_dli, bg_dli, force, stats, hist_dd, hist_dli)
        return force, stats, hist_dd, hist_dli
    return _coulomb_step_numpy(x, v, q, mass, kind, soft2, rc2, inv_vc, s_dli, bg_dli)


def _coulomb_step_numpy(x, v, q, mass, kind, soft2, rc2, inv_vc, s_dli, bg_dli, chunk=512):
    n = len(q)
    force = np.zeros((n, 3))
    rate_dli = rate_dd = pot = 0.0
    min_dli = min_tau = min_r2 = np.inf
    max_e = 0.0
    hist_dd = np.zeros(NBINS)
    hist_dli = np.zeros(NBINS)
    for start in range(0, n, chunk):
        sl = slice(start, min(n, start + chunk))
        d = x[sl, None, :] - x[None, :, :]
        r2 = np.einsum("ijk,ijk->ij", d, d)
        rows = np.arange(sl.start, sl.stop)
        r2s = r2 + soft2
        inv = 1.0 / np.sqrt(r2s)
        inv[rows - sl.start, rows] = 0.0
        c = K_COULOMB * q[sl, None] * q[None, :] * inv**3
        force[sl] = np.einsum("ij,ijk->ik", c, d)
        upper = rows[:, None] < np.arange(n)[None, :]
        pot += float(np.sum(K_COULOMB * q[sl, None] * q[None, :] * inv * upper))
        if np.any(upper):
            min_r2 = min(min_r2, float(r2[upper].min()))
        close = upper & (r2 <= rc2)
        ii, jj = np.nonzero(close)
        if len(ii) == 0:
            continue
        gi = ii + sl.start
        r = np.sqrt(r2[ii, jj])
        dli = kind[gi] != kind[jj]
        if np.any(dli):
            min_dli = min(min_dli, float(r[dli].min()))
        dv = v[gi] - v[jj]
        v2 = np.einsum("ij,ij->i", dv, dv)
        vrel = np.sqrt(v2)
        approaching = (np.einsum("ij,ij->i", d[ii, jj], dv) < 0) & (vrel > 0)
        if np.any(approaching):
            min_tau = min(min_tau, float((r[approaching] / vrel[approaching]).min()))
        mu = mass[gi] * mass[jj] / (mass[gi] + mass[jj])
        ek = 0.5 * mu * v2 / J_PER_KEV
        e_kev = ek + K_COULOMB * q[gi] * q[jj] / np.maximum(r, 1e-30) / J_PER_KEV
        ok = ek > 1e-4
        bins = np.minimum(((np.log10(np.where(ok, ek, 1.0)) + 4.0) * 8.0).astype(int), NBINS - 1)
        np.add.at(hist_dli, bins[ok & dli], vrel[ok & dli] * inv_vc)
        np.add.at(hist_dd, bins[ok & ~dli & ~((kind[gi] == 1) & (kind[jj] == 1))],
                  vrel[ok & ~dli & ~((kind[gi] == 1) & (kind[jj] == 1))] * inv_vc)
        not_lili = ~((kind[gi] == 1) & (kind[jj] == 1))
        sel = dli & not_lili
        if np.any(sel):
            max_e = max(max_e, float(e_kev[sel].max()))
            rate_dli += float(np.sum(sigma_dli_mb(e_kev[sel], s_dli, bg_dli) * MB_TO_M2 * vrel[sel] * inv_vc))
        sel = ~dli
        sel &= kind[gi] == 0
        if np.any(sel):
            rate_dd += float(np.sum(sigma_dd_mb(e_kev[sel]) * MB_TO_M2 * vrel[sel] * inv_vc))
    stats = np.array([rate_dli, rate_dd, min_dli if np.isfinite(min_dli) else 1e300,
                            max_e, min_tau if np.isfinite(min_tau) else 1e300, pot,
                      np.sqrt(min_r2) if np.isfinite(min_r2) else 1e300])
    return force, stats, hist_dd, hist_dli
