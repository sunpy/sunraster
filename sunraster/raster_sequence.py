import textwrap
import numbers

import numpy as np
import astropy.units as u
import ndcube.utils.sequence
from ndcube import NDCubeSequence

from sunraster.raster import SpectrogramABC
from sunraster import utils

__all__ = ['SpectrogramSequence', 'RasterSequence']

RASTER_AXIS_NAME = "raster scan"
SNS_AXIS_NAME = "temporal"
SLIT_STEP_AXIS_NAME = "slit step"
SLIT_AXIS_NAME = "position along slit"
SPECTRAL_AXIS_NAME = "spectral"


class SpectrogramSequence(NDCubeSequence, SpectrogramABC):
    """
    Class for holding, slicing and plotting a sequence of spectrogram cubes.

    Spectrogram cubes can be 2D or higher.

    Parameters
    ----------
    data_list: `list`
        List of `SpectrogramCube` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.

    common_axis: `int` or `None`
        If the sequence axis is aligned with an axis of the component SpectrogramCube
        instances, e.g. Spectrogram cubes have a time dimension and are arranged within
        the sequence in chronological order, set this input to the axis number of the
        time axis within the cubes.
        Default=None implies there is no common axis.

    meta: `dict` or header object (optional)
        Metadata associated with the sequence.
    """
    def __init__(self, data_list, common_axis=None, meta=None):
        # Initialize Sequence.
        super().__init__(data_list, common_axis=common_axis, meta=meta)

    @property
    def spectral(self):
        return u.Quantity([raster.spectral for raster in self.data])

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
        Applies or undoes exposure time correction to data and uncertainty and
        adjusts unit.

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
        result: `None` or `SpectrogramSequence`
            If copy=False, the original SpectrogramSequence is modified with the
            exposure time correction applied (undone).
            If copy=True, a new SpectrogramSequence is returned with the correction
            applied (undone).
        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.apply_exposure_time_correction(undo=undo,
                                                                           force=force))
        if copy is True:
            return self.__class__(
                converted_data_list, meta=self.meta, common_axis=self._common_axis)
        else:
            self.data = converted_data_list

    def __str__(self):
        return (textwrap.dedent(f"""\
                SpectrogramSequence
                -------------------
                {self._str}"""))

    @property
    def _str(self):
        """Derives summary for self.__str__."""
        if self.data[0]._time_name:
            time_period = (self.data[0].time[0].value, self.data[-1].time[-1].value)
        else:
            time_period = None
        if self.data[0]._longitude_name:
            lon_range = u.Quantity([self.lon.min(), self.lon.max()])
        else:
            lon_range = None
        if self.data[0]._latitude_name:
            lat_range = u.Quantity([self.lat.min(), self.lat.max()])
        else:
            lat_range = None
        if self.data[0]._spectral_name:
            spectral_range = u.Quantity([self.spectral.min(), self.spectral.max()])
        else:
            spectral_range = None
        return (textwrap.dedent(f"""\
                Time Range: {time_period}
                Pixel Dimensions {self.raster_instrument_axes_types}: {self.dimensions}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.data[0].unit}"""))

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"


class RasterSequence(NDCubeSequence):
    """
    Class for holding, slicing and plotting spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `Raster` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.

    common_axis: `int`
        The axis of the Raster instances corresponding to time.

    meta: `dict` or header object (optional)
        Metadata associated with the sequence.
    """
    def __init__(self, data_list, common_axis, meta=None):
        # Initialize Sequence.
        super().__init__(data_list, common_axis=common_axis, meta=meta)

        # Determine axis indices of instrument axis types.
        self._raster_axis_name = RASTER_AXIS_NAME
        self._SnS_axis_name = SNS_AXIS_NAME
        self._slit_step_axis_name = SLIT_STEP_AXIS_NAME
        self._slit_axis_name = SLIT_AXIS_NAME
        self._spectral_axis_name = SPECTRAL_AXIS_NAME
        self._single_scan_instrument_axes_types = np.empty(self.data[0].data.ndim, dtype=object)
        # Slit step axis name.
        if self._common_axis is not None:
            self._single_scan_instrument_axes_types[self._common_axis] = self._slit_step_axis_name
        # Spectral axis name.
        spectral_raster_index = np.where(np.array(self.data[0].world_axis_physical_types) ==
                                         self.data[0]._spectral_name)
        if len(spectral_raster_index) == 1:
            self._single_scan_instrument_axes_types[spectral_raster_index] = \
                    self._spectral_axis_name
        # Slit axis name.
        print(self._single_scan_instrument_axes_types)
        w = self._single_scan_instrument_axes_types == None
        if w.sum() > 1:
            print(len(w), w, w.sum())
            raise ValueError("WCS, missing_axes, and/or common_axis not consistent.")
        self._single_scan_instrument_axes_types[w] = self._slit_axis_name
        # Remove any instrument axes types whose axes are missing.
        self._single_scan_instrument_axes_types.astype(str)

    raster_dimensions = NDCubeSequence.dimensions
    SnS_dimensions = NDCubeSequence.cube_like_dimensions
    raster_world_axis_physical_types = NDCubeSequence.world_axis_physical_types
    SnS_world_axis_physical_types = NDCubeSequence.cube_like_world_axis_physical_types
    raster_axis_extra_coords = NDCubeSequence.sequence_axis_extra_coords
    SnS_axis_extra_coords = NDCubeSequence.common_axis_extra_coords
    plot_as_raster = NDCubeSequence.plot
    plot_as_SnS = NDCubeSequence.plot_as_cube

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
        Method to slice instance as though data were 4D, i.e. raster number,
        slit step position, position along slit, wavelength.
        """
        return _SequenceSlicer(self)

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if isinstance(item, tuple) and not isinstance(item[0], numbers.Integral):
            result._single_scan_instrument_axes_types = _slice_scan_axis_types(
                    result._single_scan_instrument_axes_types, item[1:])
        return result

    @property
    def raster_instrument_axes_types(self):
        return tuple([self._raster_axis_name] + list(self._single_scan_instrument_axes_types))

    @property
    def SnS_instrument_axes_types(self):
        return tuple([self._SnS_axis_name] + list(
            self._single_scan_instrument_axes_types[self._single_scan_instrument_axes_types !=
                                                    self._slit_step_axis_name]))

    @property
    def spectral(self):
        return u.Quantity([raster.spectral for raster in self.data])

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
        Applies or undoes exposure time correction to data and uncertainty and
        adjusts unit.

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

    def __str__(self):
        if self.data[0]._time_name:
            time_period = (self.data[0].time[0].value, self.data[-1].time[-1].value)
        else:
            time_period = None
        if self.data[0]._longitude_name:
            lon_range = u.Quantity([self.lon.min(), self.lon.max()])
        else:
            lon_range = None
        if self.data[0]._latitude_name:
            lat_range = u.Quantity([self.lat.min(), self.lat.max()])
        else:
            lat_range = None
        if self.data[0]._spectral_name:
            spectral_range = u.Quantity([self.spectral.min(), self.spectral.max()])
        else:
            spectral_range = None
        return (textwrap.dedent(f"""\
                RasterSequence
                --------------
                Time Range: {time_period}
                Pixel Dimensions {self.raster_instrument_axes_types}: {self.dimensions}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.data[0].unit}"""))

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"


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
        result = utils.sequence._slice_sequence_as_SnS(self.seq, item)
        if isinstance(item, tuple) and not isinstance(item[0], numbers.Integral):
            result._single_scan_instrument_axes_types = _slice_scan_axis_types(
                    self.seq._single_scan_instrument_axes_types, item)
        return result


class _SequenceSlicer:
    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):
        result = ndcube.utils.sequence.slice_sequence(self.seq, item)
        if isinstance(item, tuple) and not isinstance(item[0], numbers.Integral):
            result._single_scan_instrument_axes_types = _slice_scan_axis_types(
                    self.seq._single_scan_instrument_axes_types, item[1:])
        return result


def _slice_scan_axis_types(single_scan_axes_types, item):
    """
    Updates RasterSequence._single_scan_axes_types according to slicing.

    Parameters
    ----------
    single_scan_axes_types: `numpy.ndarray`
        Value of RasterSequence._single_scan_axes_types,
        i.e. array of strings giving type of each axis.
    
    item: `int`, `slice` or `tuple` of `slice`s.
        The slicing item that get applied to the Raster instances within the RasterSequences.

    Returns
    -------
    new_single_scan_axes_types: `numpy.ndarray`
        Update value of axis types with which to replace RasterSequence._single_scan_axes_types.

    """
    # Get boolean axes indices of axis items that aren't int,
    # i.e. axes that are not sliced away.
    not_int_axis_items = [not isinstance(axis_item, numbers.Integral) for axis_item in item]
    # Add boolean indices for axes not included in item.
    not_int_axis_items += [True] * (len(single_scan_axes_types) - len(not_int_axis_items))
    return single_scan_axes_types[np.array(not_int_axis_items)]
