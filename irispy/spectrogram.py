# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

import astropy.units as u
from ndcube import NDCubeSequence
import ndcube.cube_utils as cu
from ndcube import SequenceDimensionPair

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
