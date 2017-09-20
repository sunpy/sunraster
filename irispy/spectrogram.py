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
    """Class for holding, slicing and plotting spectrogram data.

    Parameters
    ----------
    data_list: `list`
        List of `ndcube.NDCube` objects holding data.

    common_axis: `int`
        The axis of the NDCubes corresponding to time.

    raster_positions_per_scan: `int`
        Number of slit positions per raster scan.

    first_exposure_raster_position: `int`
        The slit position of the first exposure in the data assuming zero-based counting.
        For example, a raster scan goes from left to right.  The right most position
        is designated 0.  But the data does not include the first 3 (0, 1, 2) raster
        positions in the first scan.  Therefore this variable should be set to 3.
        This enables partial raster scans to be stored in the object.

    meta: `dict` or header object
        Metadata associated with the sequence.

    """

    def __init__(self, data_list, common_axis, raster_positions_per_scan,
                 first_exposure_raster_position, meta=None):
        self.raster_positions_per_scan = raster_positions_per_scan
        self.first_exposure_raster_position = first_exposure_raster_position
        super(SpectrogramSequence, self).__init__(
            data_list, meta=meta, common_axis=common_axis)
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
    """Class for holding, slicing and plotting IRIS spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `ndcube.NDCube` objects holding data.

    common_axis: `int`
        The axis of the NDCubes corresponding to time.

    raster_positions_per_scan: `int`
        Number of slit positions per raster scan.

    first_exposure_raster_position: `int`
        The slit position of the first exposure in the data assuming zero-based counting.
        For example, a raster scan goes from left to right.  The right most position
        is designated 0.  But the data does not include the first 3 (0, 1, 2) raster
        positions in the first scan.  Therefore this variable should be set to 3.
        This enables partial raster scans to be stored in the object.

    meta: `dict` or header object
        Metadata associated with the sequence.

    """

    def to_counts(self, copy=False):
        """Converts data and uncertainty attributes to photon count units.

        Parameters
        ----------
        copy: `bool`
           If True a new instance with the converted data values is return.
           If False, the current instance is overwritten.
           Default=False

        """
        converted_data_list = _convert_iris_sequence(self, u.ct)
        if copy is True:
            return IRISSpectrogramSequence(
                converted_data_list, self._common_axis, self.raster_positions_per_scan,
                self.first_exposure_raster_position, meta=self.meta)
        else:
            self.data = converted_data_list

    def to_DN(self, copy=False):
        """Converts data and uncertainty attributes to data number (DN).

        Parameters
        ----------
        copy: `bool`
           If True a new instance with the converted data values is returned.
           If False, the current instance is overwritten.
           Default=False

        """
        converted_data_list = _convert_iris_sequence(self, "DN")
        if copy is True:
            return IRISSpectrogramSequence(
                converted_data_list, self._common_axis, self.raster_positions_per_scan,
                self.first_exposure_raster_position, meta=self.meta)
        else:
            self.data = converted_data_list

    def apply_exposure_time_correction(self, undo=False, copy=False):
        """Applies or undoes exposure time correction to data and uncertainty.

        Parameters
        ----------
        undo: `bool`
            If False, exposure time correction is applied.
            If True, exposure time correction is removed.
            Default=False

        copy: `bool`
            If True a new instance with the converted data values is returned.
            If False, the current instance is overwritten.
            Default=False

        """
        if undo is True:
            correction_function = _uncalculate_exposure_time_correction
        else:
            correction_function = _calculate_exposure_time_correction
        converted_data_list = _apply_or_undo_exposure_time_correction(
            self, correction_function)
        if copy is True:
            IRISSpectrogramSequence(
                converted_cube_list, self._common_axis, self.raster_positions_per_scan,
                self.first_exposure_raster_position, meta=self.meta)
        else:
            self.data = converted_data_list

def _convert_iris_sequence(sequence, new_unit):
    """Converts data and uncertainty in an IRISSpectrogramSequence between units.

    Parameters
    ----------
    sequence: `NDCubeSequence`, `SpectrogramSequence` or `IRISSpectrogramSequence`
        Sequence whose constituent NDCubes are be converted to new units.

    new_unit: `astropy.units.Unit` or `str`
       Unit to which the data is to be converted.

    Returns
    -------
    converted_data_list: `list` of `NDCube`s.
       List of NDCubes with data and uncertainty attributes converted to new_unit.

    """
    # Define empty list to hold NDCubes with converted data and uncertainty.
    converted_data_list = []
    # Cycle through each NDCube, convert data and uncertainty to new
    # units, and append to list.
    for i, cube in enumerate(sequence.data):
        # Determine what type of DN unit is needed based on detector type.
        detector_type = iris_tools._get_detector_type(cube.meta)
        if new_unit == "DN":
            new_unit = iris_tools.DN_UNIT[detector_type]
        # If NDCube is already in new unit, add NDCube as is to list.
        if cube.unit is new_unit or cube.unit is new_unit / u.s:
            converted_data_list.append(cube)
        # Else convert data and uncertainty to new unit.
        if cube.unit != new_unit or cube.unit != new_unit / u.s:
            # During calculations, the time component due to exposure
            # time correction, if it has been applied, is ignored.
            # Check here whether the time correction is present in the
            # original unit so that is carried through to new unit.
            if u.s not in (cube.unit.decompose() * u.s).bases:
                new_unit_time_accounted = new_unit / u.s
            else:
                new_unit_time_accounted = new_unit
            # Convert data and uncertainty to new unit.
            data = (cube.data * cube.unit).to(new_unit).value
            uncertainty = (cube.uncertainty.array * cube.unit).to(new_unit).value
            # Append new instance of NDCube in new unit to list.
            converted_data_list.append(NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask,
                unit=new_unit_time_accounted, uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords)))
    return converted_data_list

def _extra_coords_to_input_format(extra_coords):
    """
    Converts NDCube._extra_coords attribute to format required as input for new NDCube.

    Paramaters
    ----------
    extra_coords: dict
        An NDCube._extra_coords instance.


    Returns
    -------
    input_format: `list`
        Infomation on extra coords in format required by `ndcube.NDCube.__init__`.

    """
    return [(key, extra_coords[key]["axis"], extra_coords[key]["value"]) for key in extra_coords]

def _apply_or_undo_exposure_time_correction(sequence, correction_function):
    """Applies or undoes exposure time correction to a sequence of NDCubes.

    Correction is applied (or undone) to both data and uncertainty attributes of NDCubes.

    Parameters
    ----------
    sequence: `NDCubeSequence`, `SpectrogramSequence` or `IRISSpectrogramSequence`
        Sequence whose constituent NDCubes are be converted to new units.
        NDCubes with sequence must have an 'exposure time' entry in its extra
        coords attribute.

    correction_function: function
        Function applying or undoing exposure time correction.

    Returns
    -------
    converted_data_list: `list` of `NDCube`s.
       List of NDCubes with data and uncertainty corrected (or uncorrected)
       for exposure time.

    """
    for i, cube in enumerate(sequence.data):
        if u.s not in cube.unit.decompose().bases:
            exposure_time_s = cube._extra_coords["exposure time"]["value"].to(u.s).value
            if len(cube.dimensions.shape) == 1:
                pass
            elif len(cube.dimensions.shape) == 2:
                exposure_time_s = exposure_time_s[:, np.newaxis]
            elif len(cube.dimensions.shape) == 3:
                exposure_time_s = exposure_time_s[:, np.newaxis, np.newaxis]
            else:
                raise ValueError("NDCube dimensions must be 2 or 3. Dimensions={0}".format(
                    len(cube.dimensions.shape)))
            data = correction_function(cube.data, exposure_time_s)
            uncertainty = correction_function(cube.uncertainty.array, exposure_time_s)
            converted_data_list.append(NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=cube.unit / u.s,
                uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords)))
        else:
            converted_data_list.append(cube)
    return converted_data_list

def _calculate_exposure_time_correction(data, exposure_time):
    return data/exposure_time

def _uncalculate_exposure_time_correction(data, exposure_time):
    return data*exposure_time