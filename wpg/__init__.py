# -*- coding: utf-8 -*-
"""
This module contains utils and classes for X-ray wavefront propagation.

The most propagation methods based on SRW library.

.. module:: wpg
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

# Load the SRW extension before numpy/scipy pull in their own FFT libraries
# (fixes segmentation fault with FFTW from numpy MKL)
try:
    from srwpy import srwlpy  # noqa: F401
except ImportError as exc:
    raise ImportError(f'{exc}\nWPG needs SRW: python -m pip install "srwpy>=4.2". '
                      'To report a problem, run "python -m wpg_diagnose" and attach its output.') from exc

# Create aliases for simple importing
from .wavefront import Wavefront  # noqa: E402
from .beamline import Beamline  # noqa: E402
from .glossary import print_glossary  # noqa: E402


def __getattr__(name):
    # old code used ``wpg.srwlib.srwl`` after ``import wpg``
    if name == 'srwlib':
        import importlib
        return importlib.import_module('.srwlib', __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
