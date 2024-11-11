"""
This module contains package tests.
"""

from pathlib import Path

import sunraster

__all__ = ["TEST_DATA_PATH"]


TEST_DATA_PATH = Path(sunraster.__file__).parent / "tests" / "data"
