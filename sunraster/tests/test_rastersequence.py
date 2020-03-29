
import numpy as np
import pytest
from ndcube.tests.helpers import assert_cubesequences_equal
from ndcube.utils.wcs import WCS

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunraster import Raster, RasterSequence

# Define an sample wcs objects.
h0 = {
    'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10,
    'NAXIS1': 3,
    'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 2,
    'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 2}
wcs0 = WCS(header=h0, naxis=3)

# WCS with no spectral axis.
h_no_wave = {
        'CTYPE1': 'HPLT-TAN', 'CUNIT1': 'deg', 'CDELT1': 0.5, 'CRPIX1': 2, 'CRVAL1': 0.5, 'NAXIS1': 2,
        'CTYPE2': 'HPLN-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.4, 'CRPIX2': 2, 'CRVAL2': 1, 'NAXIS2': 2}
wcs_no_wave = WCS(header=h_no_wave, naxis=2)

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

# Define RasterSequences in various units.
spectrogram_DN0 = Raster(
    SOURCE_DATA_DN, wcs0, extra_coords0, u.ct, SOURCE_UNCERTAINTY_DN)
spectrogram_DN1 = Raster(
    SOURCE_DATA_DN, wcs0, extra_coords1, u.ct, SOURCE_UNCERTAINTY_DN)
sequence_DN = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq)
sequence_DN0 = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=0)
sequence_DN1 = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=1)
sequence_DN2 = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=2)

spectrogram_DN_per_s0 = Raster(
    SOURCE_DATA_DN/single_exposure_time, wcs0, extra_coords0, u.ct/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time)
spectrogram_DN_per_s1 = Raster(
    SOURCE_DATA_DN/single_exposure_time, wcs0, extra_coords1, u.ct/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time)
sequence_DN_per_s = RasterSequence(
        [spectrogram_DN_per_s0, spectrogram_DN_per_s1], meta=meta_seq)

spectrogram_DN_per_s_per_s0 = Raster(
    SOURCE_DATA_DN/single_exposure_time/single_exposure_time, wcs0, extra_coords0, u.ct/u.s/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time/single_exposure_time)
spectrogram_DN_per_s_per_s1 = Raster(
    SOURCE_DATA_DN/single_exposure_time/single_exposure_time, wcs0, extra_coords1, u.ct/u.s/u.s,
    SOURCE_UNCERTAINTY_DN/single_exposure_time/single_exposure_time)
sequence_DN_per_s_per_s = RasterSequence(
    [spectrogram_DN_per_s_per_s0, spectrogram_DN_per_s_per_s1], meta=meta_seq)

spectrogram_DN_s0 = Raster( 
    SOURCE_DATA_DN*single_exposure_time, wcs0, extra_coords0, u.ct*u.s,
    SOURCE_UNCERTAINTY_DN*single_exposure_time)
spectrogram_DN_s1 = Raster(
    SOURCE_DATA_DN*single_exposure_time, wcs0, extra_coords1, u.ct*u.s,
    SOURCE_UNCERTAINTY_DN*single_exposure_time)
sequence_DN_s = RasterSequence([spectrogram_DN_s0, spectrogram_DN_s1], meta=meta_seq)

# Define raster sequence with no spectral axes.
raster_no_wave0 = Raster(SOURCE_DATA_DN[:, :, 0], wcs_no_wave, extra_coords0, u.ct)
raster_no_wave1 = Raster(SOURCE_DATA_DN[:, :, 0], wcs_no_wave, extra_coords0, u.ct)
sequence_DN_no_wave = RasterSequence([raster_no_wave0, raster_no_wave1], 0)

# Define raster sequence with missing slit axes.
raster_no_slit0 = Raster(SOURCE_DATA_DN[:, 0], wcs0, extra_coords0, u.ct,
                         missing_axes=[False, True, False])
raster_no_slit1 = Raster(SOURCE_DATA_DN[:, 0], wcs0, extra_coords0, u.ct,
                         missing_axes=[False, True, False])
sequence_DN_no_slit = RasterSequence([raster_no_slit0, raster_no_slit1], 0)

# Define raster sequence with missing slit step axes.
raster_no_step0 = Raster(SOURCE_DATA_DN[0], wcs0, extra_coords0, u.ct,
                         missing_axes=[True, False, False])
raster_no_step1 = Raster(SOURCE_DATA_DN[0], wcs0, extra_coords0, u.ct,
                         missing_axes=[True, False, False])
sequence_DN_no_step = RasterSequence([raster_no_step0, raster_no_step1], None)


@pytest.mark.parametrize("input_sequence, undo, force, expected_sequence", [
    (sequence_DN, False, False, sequence_DN_per_s),
    (sequence_DN_per_s, True, False, sequence_DN),
    (sequence_DN_per_s, False, True, sequence_DN_per_s_per_s),
    (sequence_DN, True, True, sequence_DN_s)
])
def test_apply_exposure_time_correction(input_sequence, undo, force, expected_sequence):
    output_sequence = input_sequence.apply_exposure_time_correction(undo, copy=True, force=force)
    assert_cubesequences_equal(output_sequence, expected_sequence)


@pytest.mark.parametrize("input_sequence, expected_slit_step_axis", [
    (sequence_DN1, 2),
    #(sequence_DN1[:, 0], 1), # Won't wok until NDCubeSequence._common_axis slicing bug fixed.
    #(sequence_DN2[:, 0, 0], 1), # Won't wok until NDCubeSequence._common_axis slicing bug fixed.
    #(sequence_DN1[:, :, 0], None), # Won't wok until NDCubeSequence._common_axis slicing bug fixed.
    (sequence_DN1[:, :, :, 0], 2),
])
def test_slit_step_axis_index(input_sequence, expected_slit_step_axis):
    assert input_sequence._slit_step_axis_index == expected_slit_step_axis


@pytest.mark.parametrize("input_sequence, expected_slit_axis_raster_index", [
    (sequence_DN0, 2),
    (sequence_DN0[:, 0], 1),
    (sequence_DN0[:, :, 0], None)
])
def test_slit_axis_raster_index(input_sequence, expected_slit_axis_raster_index):
    assert input_sequence._slit_axis_raster_index == expected_slit_axis_raster_index


@pytest.mark.parametrize("input_sequence, expected_slit_axis_SnS_index", [
    (sequence_DN0, 1),
    (sequence_DN0[:, 0], 0),
    (sequence_DN0[:, :, 0], None)
])
def test_slit_axis_SnS_index(input_sequence, expected_slit_axis_SnS_index):
    assert input_sequence._slit_axis_SnS_index == expected_slit_axis_SnS_index


@pytest.mark.parametrize("input_sequence, expected_spectral_axis_raster_index", [
    (sequence_DN0, 3),
    (sequence_DN0[:, 0], 2),
    (sequence_DN0[:, :, :, 0], None)
])
def test_spectral_axis_raster_index(input_sequence, expected_spectral_axis_raster_index):
    assert input_sequence._spectral_axis_raster_index == expected_spectral_axis_raster_index


@pytest.mark.parametrize("input_sequence, expected_spectral_axis_SnS_index", [
    (sequence_DN0, 2),
    (sequence_DN0[:, 0], 1),
    (sequence_DN0[:, :, :, 0], None)
])
def test_spectral_axis_SnS_index(input_sequence, expected_spectral_axis_SnS_index):
    assert input_sequence._spectral_axis_SnS_index == expected_spectral_axis_SnS_index


@pytest.mark.parametrize("input_sequence, expected_raster_axes_types", [
    (sequence_DN0, (sequence_DN0._raster_axis_name, sequence_DN0._slit_step_axis_name,
                    sequence_DN0._slit_axis_name, sequence_DN0._spectral_axis_name)),

    (sequence_DN0[:, :, 0, 0], (sequence_DN0._raster_axis_name, sequence_DN0._slit_step_axis_name)),

    #(sequence_DN0[:, 0], (sequence_DN0._raster_axis_name, sequence_DN0._slit_axis_name,
    #                      sequence_DN0._spectral_axis_name))  # Won't wok until NDCubeSequence._common_axis slicing bug fixed.

    (sequence_DN_no_wave, (sequence_DN_no_wave._raster_axis_name,
                           sequence_DN_no_wave._slit_step_axis_name,
                           sequence_DN_no_wave._slit_axis_name)),

    (sequence_DN_no_slit, (sequence_DN_no_slit._raster_axis_name,
                           sequence_DN_no_slit._slit_step_axis_name,
                           sequence_DN_no_slit._spectral_axis_name)),

    (sequence_DN_no_step, (sequence_DN_no_step._raster_axis_name,
                           sequence_DN_no_step._slit_axis_name,
                           sequence_DN_no_step._spectral_axis_name))
])
def test_raster_axes_types(input_sequence, expected_raster_axes_types):
    assert input_sequence.raster_axes_types == expected_raster_axes_types


@pytest.mark.parametrize("input_sequence, expected_SnS_axes_types", [
    (sequence_DN0, (sequence_DN0._SnS_axis_name, sequence_DN0._slit_axis_name,
                    sequence_DN0._spectral_axis_name)),

    (sequence_DN0[:, :, 0, 0], (sequence_DN0._SnS_axis_name,)),

    #(sequence_DN0[:, 0], (sequence_DN0._SnS_axis_name, sequence_DN0._slit_axis_name,
    #                      sequence_DN0._spectral_axis_name))  # Won't wok until NDCubeSequence._common_axis slicing bug fixed.
])
def test_SnS_axes_types(input_sequence, expected_SnS_axes_types):
    assert input_sequence.SnS_axes_types == expected_SnS_axes_types
