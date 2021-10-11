import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from ndcube.tests.helpers import assert_cubes_equal

import sunraster.spectrogram
from sunraster import SpectrogramCube

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

# Define sample extra coords
EXTRA_COORDS0 = [
    ("time", 0, Time("2017-01-01") + TimeDelta(np.arange(TIME_DIM_LEN), format="sec")),
    ("exposure time", 0, EXPOSURE_TIME),
]
EXTRA_COORDS1 = [
    (
        "time",
        0,
        (Time("2017-01-01") + TimeDelta(np.arange(TIME_DIM_LEN, TIME_DIM_LEN * 2), format="sec")),
    ),
    ("exposure time", 0, EXPOSURE_TIME),
]

# Define SpectrogramCubes in various units.
spectrogram_DN0 = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS0,
    unit=u.ct,
    uncertainty=SOURCE_UNCERTAINTY_DN,
)
spectrogram_DN_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS0,
    unit=u.ct / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
)
spectrogram_DN_per_s_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS0,
    unit=u.ct / u.s / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
)
spectrogram_DN_s0 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS0,
    unit=u.ct * u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
)
spectrogram_DN1 = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS1,
    unit=u.ct,
    uncertainty=SOURCE_UNCERTAINTY_DN,
)
spectrogram_DN_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS1,
    unit=u.ct / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
)
spectrogram_DN_per_s_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS1,
    unit=u.ct / u.s / u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
)
spectrogram_DN_s1 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS1,
    unit=u.ct * u.s,
    uncertainty=SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
)
spectrogram_NO_COORDS = SpectrogramCube(SOURCE_DATA_DN, WCS_NO_COORDS)
spectrogram_instrument_axes = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs=WCS0,
    extra_coords=EXTRA_COORDS0,
    unit=u.ct,
    uncertainty=SOURCE_UNCERTAINTY_DN,
    mask=MASK,
    instrument_axes=("a", "b", "c"),
)


def test_spectral_axis():
    assert all(spectrogram_DN0.spectral_axis == spectrogram_DN0.axis_world_coords("em.wl"))


def test_spectral_axis_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.spectral_axis


def test_time():
    assert (spectrogram_DN0.time == EXTRA_COORDS0[0][2]).all()


def test_time_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.time


def test_exposure_time():
    assert (spectrogram_DN0.exposure_time == EXTRA_COORDS0[1][2]).all()


def test_exposure_time_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.exposure_time


def test_lon():
    assert (spectrogram_DN0.lon == spectrogram_DN0.axis_world_coords("lon")).all()


def test_lon_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.lon


def test_lat():
    assert (spectrogram_DN0.lat == spectrogram_DN0.axis_world_coords("lat")).all()


def test_lat_error():
    with pytest.raises(ValueError):
        spectrogram_NO_COORDS.lat


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
        sunraster.spectrogram._calculate_exposure_time_correction(SOURCE_DATA_DN, None, u.s, EXTRA_COORDS0[1][2])


def test_uncalculate_exposure_time_correction_error():
    with pytest.raises(ValueError):
        sunraster.spectrogram._uncalculate_exposure_time_correction(SOURCE_DATA_DN, None, u.ct, EXTRA_COORDS0[1][2])


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
    output = sliced_cube.instrument_axes
    assert all(output == expected)


def test_ndcube_components_after_slicing():
    """Tests all cube components are correctly propagated by slicing."""
    # Slice test object
    item = tuple([slice(0, 1)] * 3)
    sliced_cube = spectrogram_instrument_axes[item]
    # Generate expected result.
    data = spectrogram_instrument_axes.data[item]
    uncertainty = spectrogram_instrument_axes.uncertainty[item]
    mask = spectrogram_instrument_axes.mask[item]
    extra_coords = list(EXTRA_COORDS0)
    ec_axis = 0
    ec0 = list(extra_coords[0])
    ec0[-1] = ec0[-1][item[ec_axis]]
    ec1 = list(extra_coords[1])
    ec1[-1] = ec1[-1][item[ec_axis]]
    extra_coords = (tuple(ec0), tuple(ec1))
    wcs = spectrogram_instrument_axes.wcs[item]
    expected_cube = SpectrogramCube(
        data=data,
        wcs=wcs,
        uncertainty=uncertainty,
        mask=mask,
        meta=spectrogram_instrument_axes.meta,
        unit=spectrogram_instrument_axes.unit,
        extra_coords=extra_coords,
        missing_axes=spectrogram_instrument_axes.missing_axes,
        instrument_axes=spectrogram_instrument_axes.instrument_axes,
    )
    assert_cubes_equal(sliced_cube, expected_cube)
