"""Minimal stub of `celerite` so that pytransit (which hard-imports
celerite.solver.LinAlgError at package import) can be imported on
Python 3.14, where celerite has no wheel and needs MSVC to build.
Only the names pytransit touches at import time exist here; any actual
use raises. Added to sys.path explicitly by the limb_asym validation
code — NOT installed into site-packages.
"""


class GP:  # pragma: no cover - never used, import-time placeholder only
    def __init__(self, *a, **k):
        raise ImportError("celerite stub: real celerite is not installed")
