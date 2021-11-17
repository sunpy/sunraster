import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from ndcube.tests.helpers import assert_cubes_equal

import sunraster.spectrogram
from sunraster import SpectrogramCube
from sunraster.extern.meta import Meta

# Define a sample wcs object
H0 = {
    "CTYPE1": "WAVE    ",
    "CUNIT1": "Angstrom",
    "CDELT1": 0.2,
    "CRPIX1": 0,
    "CRVAL1": 10,
    "NAXIS1": 3,
    "CTYPE2": "HPLT-TAN",
    "CUNIT2": "deg",
    "CDELT2": 0.5,
    "CRPIX2": 2,
    "CRVAL2": 0.5,
    "NAXIS2": 2,
    "CTYPE3": "HPLN-TAN",
    "CUNIT3": "deg",
    "CDELT3": 0.4,
    "CRPIX3": 2,
    "CRVAL3": 1,
    "NAXIS3": 2,
}
WCS0 = WCS(header=H0, naxis=3)

H_NO_COORDS = {
    "CTYPE1": "PIX     ",
    "CUNIT1": "",
    "CDELT1": 1,
    "CRPIX1": 0,
    "CRVAL1": 0,
    "NAXIS1": 3,
    "CTYPE2": "PIX     ",
    "CUNIT2": "",
    "CDELT2": 1,
    "CRPIX2": 0,
    "CRVAL2": 0,
    "NAXIS2": 3,
    "CTYPE3": "PIX     ",
    "CUNIT3": "",
    "CDELT3": 1,
    "CRPIX3": 0,
    "CRVAL3": 0,
    "NAXIS3": 3,
}
WCS_NO_COORDS = WCS(header=H_NO_COORDS, naxis=3)

SOURCE_DATA_DN = np.array(
    [
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
    ]
)
SOURCE_UNCERTAINTY_DN = np.sqrt(SOURCE_DATA_DN)
MASK = SOURCE_DATA_DN > 1
TIME_DIM_LEN = SOURCE_DATA_DN.shape[0]
SINGLES_EXPOSURE_TIME = 2.0
EXPOSURE_TIME = u.Quantity(np.zeros(TIME_DIM_LEN) + SINGLES_EXPOSURE_TIME, unit=u.s)
EXTRA_COORDS0 = [
    ("time", 0, Time("2017-01-01") + TimeDelta(np.arange(TIME_DIM_LEN), format="sec")),
]
EXTRA_COORDS1 = [
    (
        "time",
        0,
        (Time("2017-01-01") + TimeDelta(np.arange(TIME_DIM_LEN, TIME_DIM_LEN * 2), format="sec")),
    ),
]
meta_exposure0 = Meta({"exposure time": EXPOSURE_TIME}, axes={"exposure time": 0}, data_shape=SOURCE_DATA_DN.shape)

spectrogram_DN0 = SpectrogramCube(
    SOURCE_DATA_DN, wcs=WCS0, unit=u.ct, uncertainty=SOURCE_UNCERTAINTY_DN, meta=meta_exposure0
)
spectrogram_DN0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_per_s_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct / u.s / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s_per_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_s0 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct * u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN1 = SpectrogramCube(
    SOURCE_DATA_DN, wcs=WCS0, unit=u.ct, uncertainty=SOURCE_UNCERTAINTY_DN, meta=meta_exposure0
)
spectrogram_DN1.extra_coords.add(*EXTRA_COORDS1[0])
spectrogram_DN_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s1.extra_coords.add(*EXTRA_COORDS1[0])
spectrogram_DN_per_s_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct / u.s / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s_per_s1.extra_coords.add(*EXTRA_COORDS1[0])
spectrogram_DN_s1 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    unit=u.ct * u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_s1.extra_coords.add(*EXTRA_COORDS1[0])
spectrogram_NO_COORDS = SpectrogramCube(SOURCE_DATA_DN, WCS_NO_COORDS)
spectrogram_instrument_axes = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs=WCS0,
    unit=u.ct,
    uncertainty=SOURCE_UNCERTAINTY_DN,
    mask=MASK,
    instrument_axes=("a", "b", "c"),
    meta=meta_exposure0,
)
spectrogram_instrument_axes.extra_coords.add(*EXTRA_COORDS0[0])


def test_spectral_axis():
    assert all(spectrogram_DN0.spectral_axis == spectrogram_DN0.axis_world_coords("em.wl")[0])


def test_spectral_axis_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.spectral_axis


def test_time():
    assert all(spectrogram_DN0.time == EXTRA_COORDS0[0][2])


def test_time_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.time


def test_exposure_time():
    assert all(spectrogram_DN0.exposure_time == EXPOSURE_TIME)


def test_exposure_time_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.exposure_time


@pytest.mark.parametrize(
    "input_cube, undo, force, expected_cube",
    [
        (spectrogram_DN0, False, False, spectrogram_DN_per_s0),
        (spectrogram_DN_per_s0, True, False, spectrogram_DN0),
        (spectrogram_DN_per_s0, False, True, spectrogram_DN_per_s_per_s0),
        (spectrogram_DN0, True, True, spectrogram_DN_s0),
    ],
)
def test_apply_exposure_time_correction(input_cube, undo, force, expected_cube):
    output_cube = input_cube.apply_exposure_time_correction(undo=undo, force=force)
    assert_cubes_equal(output_cube, expected_cube)


def test_calculate_exposure_time_correction_error():
    with pytest.raises(ValueError):
        sunraster.spectrogram._calculate_exposure_time_correction(SOURCE_DATA_DN, None, u.s, EXPOSURE_TIME)


def test_uncalculate_exposure_time_correction_error():
    with pytest.raises(ValueError):
        sunraster.spectrogram._uncalculate_exposure_time_correction(SOURCE_DATA_DN, None, u.ct, EXPOSURE_TIME)


@pytest.mark.parametrize(
    "item,expected",
    [
        (0, np.array(["b", "c"])),
        (slice(0, 1), np.array(["a", "b", "c"])),
        ((slice(None), 0), np.array(["a", "c"])),
        ((slice(None), slice(None), slice(0, 1)), np.array(["a", "b", "c"])),
    ],
)
def test_instrument_axes_slicing(item, expected):
    sliced_cube = spectrogram_instrument_axes[item]
    assert all(sliced_cube.instrument_axes == expected)


def test_ndcube_components_after_slicing():
    """
    Tests all cube components are correctly propagated by slicing.
    """
    item = tuple([slice(0, 1)] * 3)
    sliced_cube = spectrogram_instrument_axes[item]

    data = spectrogram_instrument_axes.data[item]
    uncertainty = spectrogram_instrument_axes.uncertainty[item]
    mask = spectrogram_instrument_axes.mask[item]
    extra_coords = list(EXTRA_COORDS0)
    ec_axis = 0
    ec0 = list(extra_coords[0])
    ec0[-1] = ec0[-1][item[ec_axis]]
    wcs = spectrogram_instrument_axes.wcs[item]
    expected_cube = SpectrogramCube(
        data=data,
        wcs=wcs,
        uncertainty=uncertainty,
        mask=mask,
        meta=spectrogram_instrument_axes.meta,
        unit=spectrogram_instrument_axes.unit,
        instrument_axes=spectrogram_instrument_axes.instrument_axes,
    )
    expected_cube.extra_coords.add(*ec0)
    assert str(sliced_cube)
    assert str(expected_cube)
    assert_cubes_equal(sliced_cube, expected_cube)
