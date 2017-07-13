# -*- coding: utf-8 -*-
import os
import pytest
import copy
import numpy as np

import astropy.wcs as wcs
from astropy.io import fits
import astropy.units as u
import six
import os
import irispy.spectrograph
import irispy.data.test

testpath = irispy.data.test.rootdir


@pytest.fixture
def iris_l2_test_raster():
    return irispy.spectrograph.IRISSpectrograph(os.path.join(testpath, 'iris_l2_20170222_153635_3690215148_raster_t000_r00000.fits'))


def test_fits_data_comparison(iris_l2_test_raster):
    """Make sure the data is the same in pyfits and irispy"""
    hdulist = fits.open(os.path.join(
        testpath, 'iris_l2_20170222_153635_3690215148_raster_t000_r00000.fits'))
    data1 = copy.deepcopy(hdulist[1].data)
    data2 = copy.deepcopy(hdulist[2].data)
    data3 = copy.deepcopy(hdulist[3].data)

    data1[hdulist[1].data == -200.] = np.nan
    data2[hdulist[2].data == -200.] = np.nan
    data3[hdulist[3].data == -200.] = np.nan

    spectral_windows = iris_l2_test_raster.spectral_windows['name']

    np.testing.assert_allclose(iris_l2_test_raster.data[spectral_windows[0]][0].data, data1)
    np.testing.assert_allclose(iris_l2_test_raster.data[spectral_windows[1]][0].data, data2)
    np.testing.assert_allclose(iris_l2_test_raster.data[spectral_windows[2]][0].data, data3)
