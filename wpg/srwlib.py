"""
Deprecated alias of :mod:`srwpy.srwlib`.

SRW is no longer bundled with WPG, it is installed as the ``srwpy`` package.
This module keeps old code like ``from wpg.srwlib import srwl`` working.
"""

import warnings

from srwpy.srwlib import *  # noqa: F401,F403
from srwpy.srwlib import srwl  # noqa: F401

warnings.warn(
    "wpg.srwlib is deprecated, import from srwpy.srwlib instead",
    DeprecationWarning,
    stacklevel=2,
)
