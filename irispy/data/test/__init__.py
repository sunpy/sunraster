"""Irispy test data files"""
from __future__ import absolute_import

import os
import glob

from astropy.utils.data import get_pkg_data_filename

import irispy

__all__ = ['rootdir', 'file_list', 'get_test_filepath']

rootdir = os.path.join(os.path.dirname(irispy.__file__), "data", "test")

def get_test_filepath(filename, **kwargs):
    """
    Return the full path to a test file in the ``data/test`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``data/test`` directory.

    Return
    ------
    filepath : `str`
        The full path to the file.

    See Also
    --------

    astropy.utils.data.get_pkg_data_filename : Get package data filename

    Notes
    -----

    This is a wrapper around `astropy.utils.data.get_pkg_data_filename` which
    sets the ``package`` kwarg to be 'irispy.data.test`.

    """
    return get_pkg_data_filename(filename, package="irispy.data.test", **kwargs)


file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))
