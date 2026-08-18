from .pipelines import *
from .utils import *
from .analyze import *
from .process import *


def __getattr__(name):
    # Forwards to the lazily loaded plotting API, see cns.analyze.
    if name.startswith("_"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from . import analyze
    try:
        return getattr(analyze, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None


def __dir__():
    from . import analyze
    return sorted(set(globals()) | set(dir(analyze)))
