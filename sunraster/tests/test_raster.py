
import numpy as np
import pytest
from ndcube.tests.helpers import assert_cubes_equal
from ndcube.utils.wcs import WCS

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunraster.raster import Raster

# Define a sample wcs object
h0 = {
    'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10,
    'NAXIS1': 3,
    'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 2,
    'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 2,
}
wcs0 = WCS(header=h0, naxis=3)

SOURCE_DATA_DN = np.array([[[ 0.563,  1.132, -1.343], [-0.719,  1.441, 1.566]],
                           [[ 0.563,  1.132, -1.343], [-0.719,  1.441, 1.566]]])
SOURCE_UNCERTAINTY_DN = np.sqrt(SOURCE_DATA_DN)

time_dim_len = SOURCE_DATA_DN.shape[0]
single_exposure_time = 2.
EXPOSURE_TIME = u.Quantity(np.zeros(time_dim_len)+single_exposure_time, unit=u.s)

# Define sample extra coords
extra_coords0 = [("time", 0,
                  Time('2017-01-01') + TimeDelta(np.arange(time_dim_len), format='sec')),
                 ("exposure time", 0, EXPOSURE_TIME)]
extra_coords1 = [("time", 0,
                  (Time('2017-01-01') +
                   TimeDelta(np.arange(time_dim_len, time_dim_len*2), format='sec'))),
                 ("exposure time", 0, EXPOSURE_TIME)]

# Define meta data
meta_seq = {"a": 0}

# Define Rasters in various units.
spectrogram_DN0 = Raster(
    SOURCE_DATA_DN, wcs0, extra_coords0, u.ct, SOURCE_UNCERTAINTY_DN)
spectrogram_DN_per_s0 = Raster(
    SOURCE_DATA_DN/single_exposure_time, wcs0, extra_coords0, u.ct/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time)
spectrogram_DN_per_s_per_s0 = Raster(
    SOURCE_DATA_DN/single_exposure_time/single_exposure_time, wcs0, extra_coords0, u.ct/u.s/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time/single_exposure_time)
spectrogram_DN_s0 = Raster(
    SOURCE_DATA_DN*single_exposure_time, wcs0, extra_coords0, u.ct*u.s,
    SOURCE_UNCERTAINTY_DN*single_exposure_time)
spectrogram_DN1 = Raster(
    SOURCE_DATA_DN, wcs0, extra_coords1, u.ct, SOURCE_UNCERTAINTY_DN)
spectrogram_DN_per_s1 = Raster(
    SOURCE_DATA_DN/single_exposure_time, wcs0, extra_coords1, u.ct/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time)
spectrogram_DN_per_s_per_s1 = Raster(
    SOURCE_DATA_DN/single_exposure_time/single_exposure_time, wcs0, extra_coords1, u.ct/u.s/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time/single_exposure_time)
spectrogram_DN_s1 = Raster(
    SOURCE_DATA_DN*single_exposure_time, wcs0, extra_coords1, u.ct*u.s,
    SOURCE_UNCERTAINTY_DN*single_exposure_time)


@pytest.mark.parametrize("input_cube, undo, force, expected_cube", [
    (spectrogram_DN0, False, False, spectrogram_DN_per_s0),
    (spectrogram_DN_per_s0, True, False, spectrogram_DN0),
    (spectrogram_DN_per_s0, False, True, spectrogram_DN_per_s_per_s0),
    (spectrogram_DN0, True, True, spectrogram_DN_s0)
])
def test_apply_exposure_time_correction(input_cube, undo, force, expected_cube):
    output_cube = input_cube.apply_exposure_time_correction(undo=undo, force=force)
    assert_cubes_equal(output_cube, expected_cube)
