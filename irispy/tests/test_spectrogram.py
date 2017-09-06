import pytest

import numpy as np
import astropy.units as u
import ndcube
from ndcube import NDCube, NDCubeSequence
from ndcube.wcs_util import WCS
from ndcube import SequenceDimensionPair

from irispy.spectrogram import SpectrogramSequence

# sample data for tests
# TODO: use a fixture reading from a test file. file TBD.
data = np.array([[[1, 2, 3, 4], [2, 4, 5, 3], [0, -1, 2, 3]],
                 [[2, 4, 5, 1], [10, 5, 2, 2], [10, 3, 3, 0]]])

data2 = np.array([[[11, 22, 33, 44], [22, 44, 55, 33], [0, -1, 22, 33]],
                  [[22, 44, 55, 11], [10, 55, 22, 22], [10, 33, 33, 0]]])

ht = {'CTYPE3': 'HPLT-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.5, 'CRPIX3': 0, 'CRVAL3': 0, 'NAXIS3': 2,
      'CTYPE2': 'WAVE    ', 'CUNIT2': 'Angstrom', 'CDELT2': 0.2, 'CRPIX2': 0, 'CRVAL2': 0, 'NAXIS2': 3,
      'CTYPE1': 'TIME    ', 'CUNIT1': 'min', 'CDELT1': 0.4, 'CRPIX1': 0, 'CRVAL1': 0, 'NAXIS1': 4}

hm = {
    'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 4,
    'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 3,
    'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 2,
}

wt = WCS(header=ht, naxis=3)
wm = WCS(header=hm, naxis=3)

cube1 = NDCube(data, wcs=wt, missing_axis=[False, False, False, True])
cube2 = NDCube(data, wcs=wm)
cube3 = NDCube(data2, wcs=wt, missing_axis=[False, False, False, True])
cube4 = NDCube(data2, wcs=wm)

seq = SpectrogramSequence([cube1, cube2, cube3, cube4], 0, [2, 2, 2, 2], 2)
seq1 = SpectrogramSequence([cube1, cube2, cube3, cube4], 0, [2, 2, 2, 2], 2)


@pytest.mark.parametrize("test_input,expected", [
    (seq.index_by_raster[0], NDCube),
    (seq.index_by_raster[1], NDCube),
    (seq.index_by_raster[2], NDCube),
    (seq.index_by_raster[3], NDCube),
    (seq.index_by_raster[0:1], NDCubeSequence),
    (seq.index_by_raster[1:3], NDCubeSequence),
    (seq.index_by_raster[0:2], NDCubeSequence),
    (seq.index_by_raster[slice(0, 2)], NDCubeSequence),
    (seq.index_by_raster[slice(0, 3)], NDCubeSequence),
])
def test_index_by_raster(test_input, expected):
    assert isinstance(test_input, expected)


@pytest.mark.parametrize("test_input,expected", [
    (seq.index_by_raster[0:1].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((2, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq.index_by_raster[1:3].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((4, 3, 4), unit=u.pix))), axis_types=('HPLN-TAN', 'HPLT-TAN', 'WAVE'))),
    (seq.index_by_raster[0:2].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((4, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq.index_by_raster[0::].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((8, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
])
def test_index_by_raster_dimensions(test_input, expected):
    for seq_indexed, expected_dim in zip(test_input.shape, expected.shape):
        assert seq_indexed.value == expected_dim.value
        assert seq_indexed.unit == expected_dim.unit
    assert test_input.axis_types == expected.axis_types


@pytest.mark.parametrize("test_input,expected", [
    (seq[0:5].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((6, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq[1:3].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((2, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq[0:6].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((6, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq[0::].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((8, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq[0:5, 0].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((9, 4), unit=u.pix))), axis_types=('WAVE', 'TIME'))),
    (seq[1:3, 0:2].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((2, 3, 4), unit=u.pix))), axis_types=('HPLT-TAN', 'WAVE', 'TIME'))),
    (seq[0:6, 0, 0:1].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((3, 4), unit=u.pix))), axis_types=('WAVE', 'TIME'))),
    (seq[0::, 0, 0].dimensions, SequenceDimensionPair(
        shape=(list(u.Quantity((16,), unit=u.pix))), axis_types=('TIME',))),
])
def test_index_by_getitem(test_input, expected):
    for seq_indexed, expected_dim in zip(test_input.shape, expected.shape):
        assert seq_indexed.value == expected_dim.value
        assert seq_indexed.unit == expected_dim.unit
    assert test_input.axis_types == expected.axis_types
