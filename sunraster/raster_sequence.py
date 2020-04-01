import textwrap

import numpy as np
import astropy.units as u
import ndcube.utils.sequence
from ndcube import NDCubeSequence

from sunraster import utils

__all__ = ['RasterSequence']


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

    meta: `dict` or header object
        Metadata associated with the sequence.

    slit_step_axis: `int`
        The axis of the Raster instances corresponding to time.
    """
    def __init__(self, data_list, slit_step_axis=0, meta=None):
        # Initialize Sequence.
        super().__init__(data_list, common_axis=slit_step_axis, meta=meta)
        self._slit_step_axis = self._common_axis

    raster_dimensions = NDCubeSequence.dimensions
    SnS_dimensions = NDCubeSequence.cube_like_dimensions
    raster_world_axis_physical_types = NDCubeSequence.world_axis_physical_types
    SnS_world_axis_physical_types = NDCubeSequence.cube_like_world_axis_physical_types
    raster_axis_extra_coords = NDCubeSequence.sequence_axis_extra_coords
    SnS_axis_extra_coords = NDCubeSequence.common_axis_extra_coords
    plot_as_raster = NDCubeSequence.plot
    plot_as_SnS = NDCubeSequence.plot_as_cube

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
                Pixel Dimensions (raster scans, slit steps, slit height, spectral): {self.dimensions}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.data[0].unit}"""))

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

