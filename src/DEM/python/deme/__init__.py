"""Python interface for DEM-Engine.

The implementation lives in a private extension module so this package can
grow Python-level helpers without changing the public import path.
"""

from ._deme import *  # noqa: F401,F403
from ._deme import __version__
