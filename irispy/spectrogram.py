# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

#from collections import namedtuple

import astropy.units as u
from ndcube import NDCubeSequence
import ndcube.cube_utils as cu
from ndcube import SequenceDimensionPair

__all__ = ['SpectrogramSequence']

class SpectrogramSequence(NDCubeSequence):
    """docstring for SpectrogramSequence"""

    def __init__(self, data_list, common_axis, raster_positions_per_scan, first_exposure_raster_position, meta=None, **kwargs):
        self.time = kwargs.get('time', None)
        self.raster_positions_per_scan = raster_positions_per_scan
        self.first_exposure_raster_position = first_exposure_raster_position
        super(SpectrogramSequence, self).__init__(
            data_list, meta=meta, common_axis=common_axis, **kwargs)

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
        return SequenceDimensionPair(shape=tuple(
            [int(sum([d.dimensions.shape[0].value for d in self.data]))]+list(self.data[0].dimensions.shape[1::])),
            axis_types=tuple(self.data[0].dimensions.axis_types))

    def axes_to_world(self, origin=0):
        list_arg = []
        indexed_not_as_one = []
        result = []
        quantity_index = 0
        missing_axis = self.data[0].missing_axis
        wcs = self.data[0].wcs
        shape = self.data[0].data.shape
        for i, _ in enumerate(missing_axis):
            # the cases where the wcs dimension was made 1 and the missing_axis is True
            if missing_axis[wcs.naxis-1-i]:
                list_arg.append(wcs.wcs.crpix[wcs.naxis-1-i]-1+origin)
            else:
                # else it is not the case where the dimension of wcs is 1.
                list_arg.append(shape[quantity_index])
                quantity_index += 1
            # appending all the indexes to be returned in the answer
                indexed_not_as_one.append(wcs.naxis-1-i)
        list_arguments = list_arg[::-1]
        pixel_to_world = wcs.all_pix2world(*list_arguments, origin)
        # collecting all the needed answer in this list.
        for index in indexed_not_as_one[::-1]:
            result.append(u.Quantity(pixel_to_world[index], unit=wcs.wcs.cunit[index]))
        return result[::-1]


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
