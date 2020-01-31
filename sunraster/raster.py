# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

import numpy as np
import astropy.units as u
from ndcube import NDCube, NDCubeSequence
from ndcube.utils.cube import convert_extra_coords_dict_to_input_format
import ndcube.utils.sequence

from sunraster import utils

__all__ = ['Raster', 'RasterSequence']

# Define some custom error messages.
APPLY_EXPOSURE_TIME_ERROR = ("Exposure time correction has probably already "
                             "been applied since the unit already includes "
                             "inverse time. To apply exposure time correction "
                             "anyway, set 'force' kwarg to True.")
UNDO_EXPOSURE_TIME_ERROR = ("Exposure time correction has probably already "
                            "been undone since the unit does not include "
                            "inverse time. To undo exposure time correction "
                            "anyway, set 'force' kwarg to True.")
AXIS_NOT_FOUND_ERROR = " axis not found. If in extra_coords, axis name must be supported: "

class RasterSequence(NDCubeSequence):
    """Class for holding, slicing and plotting spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `Raster` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.

    meta: `dict` or header object
        Metadata associated with the sequence.

    slit_step_axis: `int`
        The axis of the Raster instances corresponding to time.

    """
    def __init__(self, data_list, meta=None, slit_step_axis=0):
        # Initialize Sequence.
        super().__init__(data_list, meta=meta, common_axis=slit_step_axis)
        self._slit_step_axis = self._common_axis

    def __repr__(self):
        return (
            """RasterSequence
--------------
Start time: {start_time}
Pixel dimensions (raster scans, slit steps, slit height, spectral): {dimensions}
Longitude range: {lon_range}
Latitude range: {lat_range}
Spectral range: {spectral_range}
Data unit: {unit}
""".format(dimensions=self.dimensions,
           start_time=self.time[0],
           lon_range=u.Quantity([self.lon.min(), self.lon.max()]),
           lat_range=u.Quantity([self.lat.min(), self.lat.max()]),
           spectral_range=u.Quantity([self.spectral_axis.min(), self.spectral_axis.max()]),
           unit=self.data[0].unit)
)

    @property
    def slice_as_SnS(self):
        """
        Method to slice instance as though data were taken as a sit-and-stare,
        i.e. slit position and raster number are combined into a single axis.
        """
        return _SnSSlicer(self)

    @property
    def slice_as_raster(self):
        """
        Method to slice instance as though data were 4D,
        i.e. raster number, slit step position, position along slit, wavelength.
        """
        return _SequenceSlicer(self)

    @property
    def spectral_axis(self):
        return u.Quantity([raster.spectral_axis for raster in self.data])

    @property
    def time(self):
        return np.concatenate([raster.time for raster in self.data])

    @property
    def exposure_time(self):
        return np.concatenate([raster.exposure_time for raster in self.data])

    @property
    def lon(self):
        return u.Quantity([raster.lon for raster in self.data])

    @property
    def lat(self):
        return u.Quantity([raster.lat for raster in self.data])

    def apply_exposure_time_correction(self, undo=False, copy=False, force=False):
        """
        Applies or undoes exposure time correction to data and uncertainty and adjusts unit.

        Correction is only applied (undone) if the object's unit doesn't (does)
        already include inverse time.  This can be overridden so that correction
        is applied (undone) regardless of unit by setting force=True.

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

        force: `bool`
            If not True, applies (undoes) exposure time correction only if unit
            doesn't (does) already include inverse time.
            If True, correction is applied (undone) regardless of unit.  Unit is still
            adjusted accordingly.

        Returns
        -------
        result: `None` or `RasterSequence`
            If copy=False, the original RasterSequence is modified with the
            exposure time correction applied (undone).
            If copy=True, a new RasterSequence is returned with the correction
            applied (undone).

        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.apply_exposure_time_correction(undo=undo,
                                                                           force=force))
        if copy is True:
            return RasterSequence(
                converted_data_list, meta=self.meta, common_axis=self._common_axis)
        else:
            self.data = converted_data_list

class Raster(NDCube):
    """
    Class representing a sit-and-stare or single raster of slit spectrogram data.

    Must be described by a single WCS.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.

    wcs: `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information

    unit : `astropy.unit.Unit` or `str`
        Unit for the dataset. Strings that can be converted to a Unit are allowed.

    meta : dict-like object
        Additional meta information about the dataset.

    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn’t mandatory. If the uncertainty
        has no such attribute the uncertainty is stored as UnknownUncertainty.
        Defaults to None.

    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.

    extra_coords : iterable of `tuple`s, each with three entries
        (`str`, `int`, `astropy.units.quantity` or array-like)
        Gives the name, axis of data, and values of coordinates of a data axis not
        included in the WCS object.

    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every attribute
        before saving it while False tries to save every parameter as reference.
        Note however that it is not always possible to save the input as reference.
        Default is False.
    """
    def __init__(self, data, wcs, extra_coords, unit, uncertainty=None, meta=None,
                 mask=None, copy=False, missing_axes=None):
        # Initialize Raster.
        super().__init__(data, wcs, uncertainty=uncertainty, mask=mask, meta=meta, unit=unit,
                         extra_coords=extra_coords, copy=copy, missing_axes=missing_axes)

        # Determine labels and location of each key real world coordinate.
        supported_spectral_names = ["em.wl", "em.energy", "em.freq",
                                    "wavelength", "energy", "frequency", "freq", "lambda"]
        supported_spectral_names += [name.upper() for name in supported_spectral_names]
        supported_spectral_names += [name.capitalize() for name in supported_spectral_names]
        self._supported_spectral_names = supported_spectral_names
        self._spectral_name = self._find_axis_name(self._supported_spectral_names)

        supported_time_names = ["time"]
        supported_time_names += [name.upper() for name in supported_time_names]
        supported_time_names += [name.capitalize() for name in supported_time_names]
        self._supported_time_names = supported_time_names
        self._time_name = self._find_axis_name(self._supported_time_names)

        supported_exposure_names = [
                "exposure time", "exposure_time", "exposure times", "exposure_times",
                "exp time", "exp_time", "exp times", "exp_times"]
        supported_exposure_names += [name.upper() for name in supported_exposure_names]
        supported_exposure_names += [name.capitalize() for name in supported_exposure_names]
        self._supported_exposure_time_names = supported_exposure_names
        self._exposure_time_name = self._find_axis_name(self._supported_exposure_time_names)

        supported_longitude_names = [".lon", "longitude", "lon"]
        supported_longitude_names += [name.upper() for name in supported_longitude_names]
        supported_longitude_names += [name.capitalize() for name in supported_longitude_names]
        self._supported_longitude_names = supported_longitude_names
        self._longitude_name = self._find_axis_name(self._supported_longitude_names)

        supported_latitude_names = [".lat", "latitude", "lat"]
        supported_latitude_names += [name.upper() for name in supported_latitude_names]
        supported_latitude_names += [name.capitalize() for name in supported_latitude_names]
        self._supported_latitude_names = supported_latitude_names
        self._latitude_name = self._find_axis_name(self._supported_latitude_names)

    def __repr__(self):
        return (
            """Raster
------
Start time: {start_time}
Pixel dimensions (Slit steps, Slit height, Spectral): {dimensions}
Longitude range: {lon_range}
Latitude range: {lat_range}
Spectral range: {spectral_range}
Data unit: {unit}
""".format(dimensions=self.dimensions,
           start_time=self.time[0],
           lon_range=u.Quantity([self.lon.min(), self.lon.max()]),
           lat_range=u.Quantity([self.lat.min(), self.lat.max()]),
           spectral_range=u.Quantity([self.spectral_axis.min(), self.spectral_axis.max()]),
           unit=self.unit)
)



    def __getitem__(self, item):
        result = super().__getitem__(item)
        return Raster(
            result.data, result.wcs,
            convert_extra_coords_dict_to_input_format(result.extra_coords, result.missing_axes),
            result.unit,result.uncertainty, result.meta,
            mask=result.mask, missing_axes=result.missing_axes)

    @property
    def spectral_axis(self):
        if not self._spectral_name:
            raise ValueError("Spectral" + AXIS_NOT_FOUND_ERROR + \
                             "{0}".format(self._supported_spectral_names))
        else:
            return self._get_axis_coord(*self._spectral_name)

    @property
    def time(self):
        if not self._time_name:
            raise ValueError("Time" + AXIS_NOT_FOUND_ERROR + \
                             "{0}".format(self._supported_time_name))
        else:
            return self._get_axis_coord(*self._time_name)

    @property
    def exposure_time(self):
        if not self._exposure_time_name:
            raise ValueError("Exposure time" + AXIS_NOT_FOUND_ERROR + \
                             "{0}".format(self._supported_exposure_time_spectral_name))
        else:
            return self._get_axis_coord(*self._exposure_time_name)

    @property
    def lon(self):
        if not self._longitude_name:
            raise ValueError("Longitude" + AXIS_NOT_FOUND_ERROR + \
                             "{0}".format(self._supported_longitude_spectral_name))
        else:
            return self._get_axis_coord(*self._longitude_name)

    @property
    def lat(self):
        if not self._latitude_name:
            raise ValueError("Latitude" + AXIS_NOT_FOUND_ERROR + \
                             "{0}".format(self._supported_latitude_spectral_name))
        else:
            return self._get_axis_coord(*self._latitude_name)

    def apply_exposure_time_correction(self, undo=False, force=False):
        """
        Applies or undoes exposure time correction to data and uncertainty and adjusts unit.

        Correction is only applied (undone) if the object's unit doesn't (does)
        already include inverse time.  This can be overridden so that correction
        is applied (undone) regardless of unit by setting force=True.

        Parameters
        ----------
        undo: `bool`
            If False, exposure time correction is applied.
            If True, exposure time correction is undone.
            Default=False

        force: `bool`
            If not True, applies (undoes) exposure time correction only if unit
            doesn't (does) already include inverse time.
            If True, correction is applied (undone) regardless of unit.  Unit is still
            adjusted accordingly.

        Returns
        -------
        result: `Raster`
            New Raster in new units.

        """
        # Get exposure time in seconds and change array's shape so that
        # it can be broadcast with data and uncertainty arrays.
        exposure_time_s = self.exposure_time.to(u.s).value
        if not np.isscalar(exposure_time_s):
            if len(self.dimensions) == 1:
                pass
            elif len(self.dimensions) == 2:
                exposure_time_s = exposure_time_s[:, np.newaxis]
            elif len(self.dimensions) == 3:
                exposure_time_s = exposure_time_s[:, np.newaxis, np.newaxis]
            else:
                raise ValueError(
                    "Raster dimensions must be 2 or 3. Dimensions={0}".format(
                        len(self.dimensions.shape)))
        # Based on value on undo kwarg, apply or remove exposure time correction.
        if undo is True:
            new_data_arrays, new_unit = _uncalculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        else:
            new_data_arrays, new_unit = _calculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        # Return new instance of Raster with correction applied/undone.
        return Raster(
            new_data_arrays[0], self.wcs,
            convert_extra_coords_dict_to_input_format(self.extra_coords, self.missing_axes),
            new_unit, new_data_arrays[1], self.meta, mask=self.mask, missing_axes=self.missing_axes)

    def _find_axis_name(self, supported_names):
        axis_name = None
        n_names = len(supported_names)
        i = 0
        while axis_name is None:
            if i >= n_names:
                break
            # Check WCS.
            wcs_name_index = ([supported_names[i] in world_axis_type
                               for world_axis_type in self.world_axis_physical_types])
            if sum(wcs_name_index) == 1:
                wcs_name_index = \
                        int(np.arange(len(self.world_axis_physical_types))[wcs_name_index])
                axis_name = self.world_axis_physical_types[wcs_name_index]
                loc = "wcs"

            # If label not contained in WCS, check extra coords.
            if axis_name is None:
                if supported_names[i] in self.extra_coords.keys():
                    axis_name = supported_names[i]
                    loc = "extra_coords"

        if axis_name is None:
            return axis_name
        else:
            return (axis_name, loc)

    def _get_axis_coord(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            return self.axis_world_coords(axis_name)
        elif coord_loc == "extra_coords":
            return self.extra_coords[axis_name]["value"]


def _calculate_exposure_time_correction(old_data_arrays, old_unit, exposure_time,
                                        force=False):
    """
    Applies exposure time correction to data arrays.
    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.
    old_unit: `astropy.unit.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.
    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays with exposure time corrected for.
    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays after exposure time correction.
    """
    # If force is not set to True and unit already includes inverse time,
    # raise error as exposure time correction has probably already been
    # applied and should not be applied again.
    if force is not True and u.s in old_unit.decompose().bases:
        raise ValueError(APPLY_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does not include inverse time and so
        # exposure does need to be applied, or
        # user has set force=True and wants the correction applied
        # regardless of the unit.
        new_data_arrays = [old_data/exposure_time for old_data in old_data_arrays]
        new_unit = old_unit/u.s
    return new_data_arrays, new_unit


def _uncalculate_exposure_time_correction(old_data_arrays, old_unit,
                                          exposure_time, force=False):
    """
    Removes exposure time correction from data arrays.
    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.
    old_unit: `astropy.unit.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.
    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays with exposure time correction removed.
    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays after exposure time correction removed.
    """
    # If force is not set to True and unit does not include inverse time,
    # raise error as exposure time correction has probably already been
    # undone and should not be undone again.
    if force is not True and u.s in (old_unit*u.s).decompose().bases:
        raise ValueError(UNDO_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does include inverse time and so
        # exposure does need to be removed, or
        # user has set force=True and wants the correction removed
        # regardless of the unit.
        new_data_arrays = [old_data * exposure_time for old_data in old_data_arrays]
        new_unit = old_unit*u.s
    return new_data_arrays, new_unit


class _SnSSlicer:
    """
    Helper class to make slicing in index_as_cube sliceable/indexable like a
    numpy array.
    Parameters
    ----------
    seq : `ndcube.NDCubeSequence`
        Object of NDCubeSequence.
    """

    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):
        return utils.sequence._slice_sequence_as_SnS(self.seq, item)


class _SequenceSlicer:
    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):
        return ndcube.utils.sequence.slice_sequence(self.seq, item)
