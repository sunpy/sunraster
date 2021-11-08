import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from ndcube.tests.helpers import assert_cubesequences_equal

from sunraster import RasterSequence, SpectrogramCube, SpectrogramSequence
from sunraster.extern.meta import Meta

# Define an sample wcs objects.
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

h2 = {
    "CTYPE1": "HPLN-TAN",
    "CUNIT1": "deg",
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
    "CTYPE3": "WAVE    ",
    "CUNIT3": "Angstrom",
    "CDELT3": 0.4,
    "CRPIX3": 2,
    "CRVAL3": 1,
    "NAXIS3": 2,
}
wcs2 = WCS(header=h2, naxis=3)

# WCS with no spectral axis.
h_no_wave = {
    "CTYPE1": "HPLT-TAN",
    "CUNIT1": "deg",
    "CDELT1": 0.5,
    "CRPIX1": 2,
    "CRVAL1": 0.5,
    "NAXIS1": 2,
    "CTYPE2": "HPLN-TAN",
    "CUNIT2": "deg",
    "CDELT2": 0.4,
    "CRPIX2": 2,
    "CRVAL2": 1,
    "NAXIS2": 2,
}
wcs_no_wave = WCS(header=h_no_wave, naxis=2)

SOURCE_DATA_DN = np.array(
    [
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
    ]
)
SOURCE_UNCERTAINTY_DN = np.sqrt(SOURCE_DATA_DN)

TIME_DIM_LEN = SOURCE_DATA_DN.shape[0]
SINGLES_EXPOSURE_TIME = 2.0
EXPOSURE_TIME = u.Quantity(np.zeros(TIME_DIM_LEN) + SINGLES_EXPOSURE_TIME, unit=u.s)

# Define sample extra coords
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

extra_coords20 = [
    (
        "time",
        2,
        Time("2017-01-01") + TimeDelta(np.arange(SOURCE_DATA_DN.shape[2]), format="sec"),
    ),
]
extra_coords21 = [
    (
        "time",
        2,
        (
            Time("2017-01-01")
            + TimeDelta(
                np.arange(SOURCE_DATA_DN.shape[2], SOURCE_DATA_DN.shape[2] * 2),
                format="sec",
            )
        ),
    ),
]

# Define meta data
meta_seq = {
    "a": 0,
}
meta_exposure0 = Meta({"exposure time": EXPOSURE_TIME}, axes={"exposure time": 0}, data_shape=SOURCE_DATA_DN.shape)
meta_exposure2 = Meta(
    {"exposure time": u.Quantity(np.zeros(SOURCE_DATA_DN.shape[2]) + SINGLES_EXPOSURE_TIME, unit=u.s)},
    axes={"exposure time": 2},
    data_shape=SOURCE_DATA_DN.shape,
)

# Define RasterSequences in various units.
spectrogram_DN0 = SpectrogramCube(
    SOURCE_DATA_DN,
    WCS0,
    u.ct,
    SOURCE_UNCERTAINTY_DN,
    meta=meta_exposure0,
)
spectrogram_DN0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN1 = SpectrogramCube(
    SOURCE_DATA_DN,
    WCS0,
    u.ct,
    SOURCE_UNCERTAINTY_DN,
    meta=meta_exposure0,
)
spectrogram_DN1.extra_coords.add(*EXTRA_COORDS1[0])
sequence_DN = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=0)
sequence_DN0 = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=0)
sequence_DN1 = RasterSequence([spectrogram_DN0, spectrogram_DN1], meta=meta_seq, common_axis=1)


spectrogram_DN20 = SpectrogramCube(SOURCE_DATA_DN, wcs2, u.ct, SOURCE_UNCERTAINTY_DN, meta=meta_exposure2)
spectrogram_DN20.extra_coords.add(*extra_coords20[0])
spectrogram_DN21 = SpectrogramCube(SOURCE_DATA_DN, wcs2, u.ct, SOURCE_UNCERTAINTY_DN, meta=meta_exposure2)
spectrogram_DN21.extra_coords.add(*extra_coords21[0])
sequence_DN2 = RasterSequence([spectrogram_DN20, spectrogram_DN21], meta=meta_seq, common_axis=2)

spectrogram_DN_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct / u.s,
    SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct / u.s,
    SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s1.extra_coords.add(*EXTRA_COORDS1[0])
sequence_DN_per_s = RasterSequence([spectrogram_DN_per_s0, spectrogram_DN_per_s1], meta=meta_seq, common_axis=0)
spectrogram_DN_per_s_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct / u.s / u.s,
    SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s_per_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_per_s_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct / u.s / u.s,
    SOURCE_UNCERTAINTY_DN / SINGLES_EXPOSURE_TIME / SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_per_s_per_s1.extra_coords.add(*EXTRA_COORDS1[0])
sequence_DN_per_s_per_s = RasterSequence(
    [spectrogram_DN_per_s_per_s0, spectrogram_DN_per_s_per_s1],
    meta=meta_seq,
    common_axis=0,
)
spectrogram_DN_s0 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct * u.s,
    SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_s0.extra_coords.add(*EXTRA_COORDS0[0])
spectrogram_DN_s1 = SpectrogramCube(
    SOURCE_DATA_DN * SINGLES_EXPOSURE_TIME,
    WCS0,
    u.ct * u.s,
    SOURCE_UNCERTAINTY_DN * SINGLES_EXPOSURE_TIME,
    meta=meta_exposure0,
)
spectrogram_DN_s1.extra_coords.add(*EXTRA_COORDS1[0])
sequence_DN_s = RasterSequence([spectrogram_DN_s0, spectrogram_DN_s1], meta=meta_seq, common_axis=0)

# Define raster sequence with no spectral axes.
raster_no_wave0 = SpectrogramCube(SOURCE_DATA_DN[:, :, 0], wcs_no_wave, u.ct, meta=meta_exposure0)
raster_no_wave0.extra_coords.add(*EXTRA_COORDS0[0])
raster_no_wave1 = SpectrogramCube(SOURCE_DATA_DN[:, :, 0], wcs_no_wave, u.ct, meta=meta_exposure0)
raster_no_wave1.extra_coords.add(*EXTRA_COORDS0[0])
sequence_DN_no_wave = RasterSequence([raster_no_wave0, raster_no_wave1], 0)

# Define raster sequence with missing slit axes.
raster_no_slit0 = raster_no_slit1 = spectrogram_DN0[:, 0]
sequence_DN_no_slit = RasterSequence([raster_no_slit0, raster_no_slit1], 0)

# Define raster sequence with missing slit step axes.
raster_no_step0 = raster_no_step1 = spectrogram_DN0[0]
sequence_DN_no_step = RasterSequence([raster_no_step0, raster_no_step1], None)


def test_spectral_axis():
    spectral_axis = u.Quantity([d.spectral_axis for d in sequence_DN.data])
    assert (sequence_DN.spectral_axis == spectral_axis).all()


def test_time():
    times = np.concatenate([d.time for d in sequence_DN.data])
    assert (sequence_DN.time == times).all()


def test_exposure_time():
    exposure_time = np.concatenate([d.exposure_time for d in sequence_DN.data])
    assert (sequence_DN.exposure_time == exposure_time).all()


@pytest.mark.parametrize(
    "input_sequence, undo, force, expected_sequence",
    [
        (sequence_DN, False, False, sequence_DN_per_s),
        (sequence_DN_per_s, True, False, sequence_DN),
        (sequence_DN_per_s, False, True, sequence_DN_per_s_per_s),
        (sequence_DN, True, True, sequence_DN_s),
    ],
)
def test_apply_exposure_time_correction(input_sequence, undo, force, expected_sequence):
    output_sequence = input_sequence.apply_exposure_time_correction(undo, copy=True, force=force)
    assert_cubesequences_equal(output_sequence, expected_sequence)


@pytest.mark.parametrize(
    "input_sequence, expected_raster_axes_types",
    [
        (
            sequence_DN0,
            (
                sequence_DN0._raster_axis_name,
                sequence_DN0._slit_step_axis_name,
                sequence_DN0._slit_axis_name,
                sequence_DN0._spectral_axis_name,
            ),
        ),
        (
            sequence_DN0.slice_as_raster[:, :, 0, 0],
            (sequence_DN0._raster_axis_name, sequence_DN0._slit_step_axis_name),
        ),
        (
            sequence_DN_no_wave,
            (
                sequence_DN_no_wave._raster_axis_name,
                sequence_DN_no_wave._slit_step_axis_name,
                sequence_DN_no_wave._slit_axis_name,
            ),
        ),
        (
            sequence_DN_no_slit,
            (
                sequence_DN_no_slit._raster_axis_name,
                sequence_DN_no_slit._slit_step_axis_name,
                sequence_DN_no_slit._spectral_axis_name,
            ),
        ),
        (
            sequence_DN_no_step,
            (
                sequence_DN_no_step._raster_axis_name,
                sequence_DN_no_step._slit_axis_name,
                sequence_DN_no_step._spectral_axis_name,
            ),
        ),
    ],
)
def test_raster_instrument_axes_types(input_sequence, expected_raster_axes_types):
    print(input_sequence.raster_instrument_axes_types, expected_raster_axes_types)
    assert input_sequence.raster_instrument_axes_types == expected_raster_axes_types


@pytest.mark.parametrize(
    "input_sequence, expected_sns_axes_types",
    [
        (
            sequence_DN0,
            (
                sequence_DN0._sns_axis_name,
                sequence_DN0._slit_axis_name,
                sequence_DN0._spectral_axis_name,
            ),
        ),
        (sequence_DN0.slice_as_raster[:, :, 0, 0], (sequence_DN0._sns_axis_name,)),
    ],
)
def test_sns_instrument_axes_types(input_sequence, expected_sns_axes_types):
    assert input_sequence.sns_instrument_axes_types == expected_sns_axes_types


def test_slice_as_raster():
    assert isinstance(sequence_DN[:, 0], SpectrogramSequence)
