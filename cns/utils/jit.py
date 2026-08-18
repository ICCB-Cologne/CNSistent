"""Optional numba acceleration.

The JIT-compiled kernels in CNSistent are plain Python and numpy, so they also run
without numba installed, at roughly half the aggregation speed. Numba is a default
dependency, but it (and llvmlite) cannot be installed everywhere, so the decorators
fall back to no-ops instead of breaking the import of the whole package.
"""

from cns.utils.logging import log_warn

try:
    from numba import jit, njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def _no_jit(*args, **kwargs):
        # Handles both the bare `@njit` and the called `@jit(nopython=True)` form.
        if len(args) == 1 and not kwargs and callable(args[0]):
            return args[0]
        return lambda func: func

    jit = _no_jit
    njit = _no_jit

    log_warn("numba is not installed, running uncompiled kernels (aggregation is about 2x slower).")
