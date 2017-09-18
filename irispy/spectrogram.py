# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

import numpy as np
import astropy.units as u
from ndcube import NDCube, NDCubeSequence
import ndcube.cube_utils as cu
from ndcube import SequenceDimensionPair

from irispy import iris_tools

__all__ = ['SpectrogramSequence']

class SpectrogramSequence(NDCubeSequence):
    """docstring for SpectrogramSequence"""

    def __init__(self, data_list, common_axis, raster_positions_per_scan, first_exposure_raster_position,
                 meta=None, **kwargs):
        self.raster_positions_per_scan = raster_positions_per_scan
        self.first_exposure_raster_position = first_exposure_raster_position
        super(SpectrogramSequence, self).__init__(
            data_list, meta=meta, common_axis=common_axis, **kwargs)
        self.exposure_axis_extra_coords = self._common_axis_extra_coords

    def __getitem__(self, item):
        if item is None or (isinstance(item, tuple) and None in item):
            raise IndexError("None indices not supported")
        # need to slice the time here
        return self.index_as_cube[item]

    @property
    def index_by_raster(self):
        return _IndexByRasterSlicer(self)

    @property
    def dimensions(self):
        return SequenceDimensionPair(
            shape=tuple([int(sum([d.dimensions.shape[0].value for d in self.data]))] + \
                        list(self.data[0].dimensions.shape[1::])),
            axis_types=tuple(self.data[0].dimensions.axis_types))

    def __repr__(self):
        return(
            """SpectrogramSequence
---------------------
Rasters:  {n_rasters}
Exposures per Raster: {n_steps}
Axis Types: {axis_types}
Sequence Shape: {seq_shape}\n
""".format(n_rasters=int((self.dimensions.shape[0]+self.first_exposure_raster_position)/ \
                         self.raster_positions_per_scan),
           n_steps=self.raster_positions_per_scan, axis_types=self.dimensions.axis_types[::],
           seq_shape=self.dimensions.shape))

class _IndexByRasterSlicer(object):
    """
    Helper class to make slicing in index_as_cube more pythonic.
    Helps to make operations like in numpy array.
    >>> import numpy as np
    >>> data_list = np.array(range(10))
    >>> data_list[3:5]
    >>> [4, 5, 6]
    This makes slicing like this possible for index_as_cube.

    Attributes
    ----------
    seq : `ndcube.ndcube.NDCubeSequence`
        Object of NDCubeSequence.
    """

    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):
        # slice time here
        return cu.get_cube_from_sequence(self.seq, item)

class IRISSpectrogramSequence(SpectrogramSequence):
    """A SpectrogramSequence for IRIS data including additional functionalities."""
    def to_counts(self):
        """Converts data and uncertainty to photon count units."""
        for i, cube in enumerate(self.data):
            if "FUV" in cube.meta["detector type"]:
                DN_unit = iris_tools.DN_UNIT["FUV"]
            elif cube.meta["detector type"] == "NUV":
                DN_unit = iris_tools.DN_UNIT["NUV"]
            else:
                raise ValueError("Detector type in FITS header not recognized.")
            if cube.unit == DN_unit or cube.unit == DN_unit/u.s:
                data = (cube.data*DN_unit).to(u.ct).value
                uncertainty = (cube.uncertainty.array*DN_unit).to(u.ct).value
                self.data[i] = NDCube(
                    data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=u.ct,
                    uncertainty=uncertainty,
                    extra_coords=_extra_coords_to_input_format(cube._extra_coords))
            elif cube.unit != u.ct and cube.unit != u.ct/u.s:
                raise TypeError("Data unit of {0}th cube in sequence incompatible "
                                "({1}) with {2} unit.".format(i, cube.unit, u.ct))

    def to_DN(self):
        """Converts data and uncertainty to photon units."""
        for i, cube in enumerate(self.data):
            if "FUV" in cube.meta["detector type"]:
                DN_unit = iris_tools.DN_UNIT["FUV"]
            elif cube.meta["detector type"] == "NUV":
                DN_unit = iris_tools.DN_UNIT["NUV"]
            else:
                raise ValueError("Detector type in FITS header not recognized.")
            if cube.unit == u.ct or cube.unit == u.ct/u.s:
                data = (cube.data * u.ct).value
                uncertainty = (cube.uncertainty.array * u.ct).value
                self.data[i] = NDCube(
                    data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=DN_unit,
                    uncertainty=uncertainty,
                    extra_coords=_extra_coords_to_input_format(cube._extra_coords))
            elif cube.unit != DN_unit and cube.unit != DN_unit/u.s:
                raise TypeError("Data unit of {0}th cube in sequence incompatible "
                                "({1}) with {2} unit.".format(i, cube.unit, DN_unit))

    def apply_exposure_time_correction(self):
        """Applies or undoes exposure time correction to data and uncertainty."""
        for i, cube in enumerate(self.data):
            exposure_time_s = cube._extra_coords["exposure time"]["value"].to(u.s).value
            data = cube.data / exposure_time_s[:, np.newaxis, np.newaxis]
            uncertainty = cube.uncertainty.array / exposure_time_s[:, np.newaxis, np.newaxis]
            self.data[i] = NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=cube.unit/u.s,
                uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords))

    def undo_exposure_time_correction(self):
        """Removes exposure time correction from data and uncertainty."""
        for i, cube in enumerate(self.data):
            exposure_time_s = cube._extra_coords["exposure time"]["value"].to(u.s).value
            data = cube.data * exposure_time_s[:, np.newaxis, np.newaxis]
            uncertainty = cube.uncertainty.array * exposure_time_s[:, np.newaxis, np.newaxis]
            self.data[i] = NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=cube.unit*u.s,
                uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords))


def _extra_coords_to_input_format(extra_coords):
    """
    Converts NDCube._extra_coords attribute to format required as input for new NDCube.

    Paramaters
    ----------
    extra_coords: dict
        An NDCube._extra_coords instance.
    """
    return [(key, extra_coords[key]["axis"], extra_coords[key]["value"]) for key in extra_coords]