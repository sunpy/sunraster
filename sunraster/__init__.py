"""
sunraster
=========

Does stuff.
"""
import sys

from .spectrogram import SpectrogramCube
from .spectrogram_sequence import RasterSequence, SpectrogramSequence
from .version import version as __version__

# Enforce Python version check during package import.
__minimum_python_version__ = "3.7"

__all__ = ["SpectrogramCube", "SpectrogramSequence", "RasterSequence"]


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split(".")):
    # This has to be .format to keep backwards compatibly.
    raise UnsupportedPythonError(f"sunraster does not support Python < {__minimum_python_version__}")
