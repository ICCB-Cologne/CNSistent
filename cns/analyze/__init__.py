import importlib

from .aneuploidy import *
from .breakage import *
from .coverage import *
from .distance import *

# The plotting API is not imported eagerly, as it pulls in matplotlib, which is only
# an optional dependency (the `plot` extra). It is loaded on first attribute access,
# so that `cns.analyze.fig_lines` keeps working as before.


def _plot_module():
    # Uses importlib rather than `from . import plot`, as the latter probes this
    # module for the `plot` attribute and would recurse back into __getattr__.
    return importlib.import_module(__name__ + ".plot")


def __getattr__(name):
    if name.startswith("_"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    plot = _plot_module()
    if not hasattr(plot, name):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(plot, name)


def __dir__():
    names = set(globals())
    try:
        plot = _plot_module()
    except ImportError:  # matplotlib is not installed
        return sorted(names)
    return sorted(names | {n for n in vars(plot) if not n.startswith("_")})
