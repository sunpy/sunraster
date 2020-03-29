import textwrap

import numpy as np
import astropy.units as u
from ndcube import NDCubeSequence
import ndcube.utils.sequence

import sunraster.utils

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

    common_axis: `int`
        The axis of the Raster instances corresponding to the slit step/time.
    """
    def __init__(self, data_list, common_axis=0, meta=None):
        # Initialize Sequence.
        super().__init__(data_list, common_axis=common_axis, meta=meta)
        self.sequence_axis = 0
        self._raster_axis_physical_type = self.world_axis_physical_types[self.sequence_axis]
        self._raster_axis_name = "raster"
        self._SnS_axis_physical_type = "time"
        self._SnS_axis_name = "temporal"
        self._slit_step_axis_name = "slit step"
        self._slit_axis_name = "position along slit"
        self._spectral_axis_name = "spectral"
        self._raster_axes_names = np.array([self._raster_axis_name, self._slit_step_axis_name,
                                            self._slit_axis_name, self._spectral_axis_name])
        self._SnS_axes_names = np.array([self._SnS_axis_name, self._slit_axis_name,
                                         self._spectral_axis_name])

    raster_dimensions = NDCubeSequence.dimensions
    SnS_dimensions = NDCubeSequence.cube_like_dimensions
    raster_world_axis_physical_types = NDCubeSequence.world_axis_physical_types
    SnS_world_axis_physical_types = NDCubeSequence.cube_like_world_axis_physical_types
    raster_axis_extra_coords = NDCubeSequence.sequence_axis_extra_coords
    SnS_axis_extra_coords = NDCubeSequence.common_axis_extra_coords
    plot_as_raster = NDCubeSequence.plot
    plot_as_SnS = NDCubeSequence.plot_as_cube

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

    @property
    def raster_axes_types(self):
        # Get array of indices for each axis.
        axes_indices = [self._raster_axis_index, self._slit_step_axis_index,
                        self._slit_axis_raster_index, self._spectral_axis_raster_index]
        # Remove any Nones.  These correspond to missing axes.
        axes_indices = np.array(list(filter((None).__ne__, axes_indices)))
        return tuple(self._raster_axes_names[axes_indices])

    @property
    def SnS_axes_types(self):
        # Get array of indices for each axis.
        axes_indices = [self._SnS_axis_index, self._slit_axis_SnS_index,
                        self._spectral_axis_SnS_index]
        # Remove any Nones.  These correspond to missing axes.
        axes_indices = np.array(list(filter((None).__ne__, axes_indices)))
        return tuple(self._SnS_axes_names[axes_indices])

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
                converted_data_list, meta=self.meta, common_axis=self._SnS_axis_index)
        else:
            self.data = converted_data_list

    def __str__(self):
        if self.data[0]._time_name:
            times = self.data[0].time
            if times.isscalar:
                time_period = times.value
            else:
                time_period = (times[0].value, times[-1].value)
        else:
            time_period = None
        if self.data[0]._longitude_name:
            lon = self.lon
            if lon.isscalar:
                lon_range = lon
            else:
                lon_range = u.Quantity([lon.min(), lon.max()])
        else:
            lon_range = None
        if self.data[0]._latitude_name:
            lat = self.lat
            if lat.isscalar:
                lat_range = self.lat
            else:
                lat_range = u.Quantity([lat.min(), lat.max()])
        else:
            lat_range = None
        if self.data[0]._spectral_name:
            spectral = self.spectral
            if spectral.isscalar:
                spectral_range = spectral
            else:
                spectral_range = u.Quantity([spectral.min(), spectral.max()])
        else:
            spectral_range = None
        return (textwrap.dedent(f"""\
                RasterSequence
                --------------
                Time Range: {time_period}
                Axes Types: {self.raster_axes_types}
                Axes Dimensions: {self.dimensions}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.data[0].unit}"""))

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    @property
    def _raster_axis_index(self):
        return self.sequence_axis

    @property
    def _slit_step_axis_index(self):
        if self._common_axis is not None:
            return self._common_axis + 1

    @property
    def _SnS_axis_index(self):
        return self._common_axis

    @property
    def _slit_axis_raster_index(self):
        return self._slit_axis_index_partial(self._world_types_to_raster_axes,
                                             self._slit_step_axis_index)

    @property
    def _slit_axis_SnS_index(self):
        slit_axis_SnS_index = self._slit_axis_index_partial(self._world_types_to_SnS_axes,
                                                            self._SnS_axis_index)
        # If slit step axis is sliced away, slit index has been incorrectly
        # decremented since the raster axis also counts as the 0th SnS axis.
        # Reincrement if necessary.
        if self._slit_step_axis_index is None:
            slit_axis_SnS_index += 1
        return slit_axis_SnS_index

    @property
    def _spectral_axis_raster_index(self):
        spectral_axis_index = self._axis_index_partial(self.data[0]._spectral_name,
                                                       self._world_types_to_raster_axes)
        return spectral_axis_index

    @property
    def _spectral_axis_SnS_index(self):
        spectral_axis_index = self._axis_index_partial(self.data[0]._spectral_name,
                                                       self._world_types_to_SnS_axes)
        # If slit step axis is sliced away, slit index has been incorrectly
        # decremented since the raster axis also counts as the 0th SnS axis.
        # Reincrement if necessary.
        if self._slit_step_axis_index is None:
            spectral_axis_index += 1
        return spectral_axis_index

    def _slit_axis_index_partial(self, world_types_to_axes_func, index_to_remove):
        non_unique_name_error_fragment = "not unique to a physical axis type"
        try:
            slit_axis_index = self._axis_index_partial(self.data[0]._latitude_name,
                                                       world_types_to_axes_func)
        except ValueError as err1:
            try:
                if non_unique_name_error_fragment in err.args[0]:
                    slit_axis_index = self._axis_index_partial(self.data[0]._longitude_name,
                                                               world_types_to_axes_func)
            except ValueError as err2:
                if non_unique_name_error_fragment in err.args[0]:
                    slit_axis_index = None
                else:
                    raise err2
        if isinstance(slit_axis_index, tuple):
            slit_axis_index = list(slit_axis_index)
            slit_axis_index.remove(index_to_remove)
            if slit_axis_index:
                return slit_axis_index[0]
        else:
            return slit_axis_index

    def _axis_index_partial(self, _physical_type_name, world_types_to_axes_func):
        # If _physical_type_name is None, then axis must have been sliced away.
        if _physical_type_name:
            physical_type, location = _physical_type_name
            if location == "extra_coords":
                include_extra_coords = True
            else:
                include_extra_coords = False
            axes = world_types_to_axes_func(physical_type, include_extra_coords=include_extra_coords)
            if len(axes) != 1:
                raise TypeError(f"Unrecognized output of world_types_to_axes_func: {axes}.\n"
                                "Must be a length 1 tuple.")
            return axes[0]

    def _raster_axes_to_world_types(self, *axes, include_extra_coords=True):
        """
        Retrieve the world axis physical types for each pixel axis given in raster format.

        This differs from world_axis_physical_types in that it provides an explicit
        mapping between pixel axes and physical types, including dependent physical
        types.

        Parameters
        ----------
        axes: `int` or multiple `int`
            Axis number in numpy ordering of axes for which real world physical types
            are desired.
            axes=None implies axis names for all axes will be returned.

        include_extra_coords: `bool`
           If True, also search extra_coords for coordinate corresponding to axes.
           Default=True.

        Returns
        -------
        axes_names: `tuple` of `str`
            The world axis physical types corresponding to each axis.
            If more than one physical type found for an axis, that axis's entry will
            be a tuple of `str`.
        """
        # Parse user input.
        if axes == ():
            axes = np.arange(len(self.dimensions))
        elif isinstance(axes, int):
            axes = np.array([axes])
        else:
            axes = np.array(axes)

        n_axes = len(axes)
        axes_names = [None] * n_axes

        # If sequence axis in axes, get names for it separately.
        if self._raster_axis_index in axes:
            raster_names_indices = np.array([axis == self._raster_axis for axis in axes])
            # Get standard sequence axis name from world_axis_physical_types.
            raster_axis_names = [self._raster_axis_physical_type]
            # If desired, get extra coord sequence axis names.
            if include_extra_coords:
                extra_sequence_names = ndcube.utils.sequence._get_axis_extra_coord_names_and_units(
                        self.data, None)[0]
                if extra_sequence_names:
                    raster_axis_names += list(extra_sequence_names)
            # Enter sequence axis names into output.
            # Must use for loop as can't assign tuples to multiple array location
            # with indexing and setitem.
            for i in np.arange(n_axes)[raster_names_indices]:
                axes_names[i] = tuple(raster_axis_names)

            # Get indices of axes numbers associated with cube axes.
            cube_indices = np.invert(raster_names_indices)
            cube_axes = axes[cube_indices] - 1
        else:
            cube_indices = np.ones(n_axes, dtype=bool)
            cube_axes = axes

        # Get world types from cube axes.
        if len(cube_axes) > 0:
            cube_axes_names = self.data[0]._pixel_axes_to_world_types(
                    *cube_axes, include_extra_coords=include_extra_coords)
            for i, name in zip(np.arange(n_axes)[cube_indices], cube_axes_names):
                axes_names[i] = name

        return tuple(axes_names)

    def _SnS_axes_to_world_types(self, *axes, include_extra_coords=True):
        return self.data[0]._pixel_axes_to_world_types(
                *axes, include_extra_coords=include_extra_coords)

    def _world_types_to_raster_axes(self, *axes_names, include_extra_coords=True):
        """
        Retrieve the pixel axes (numpy ordering) corresponding to each world axis physical type.

        Parameters
        ----------
        axes_names: `str` or multiple `str`
            world axis physical types for which the pixel axis numbers are desired.
            axes_names=None implies all axes will be returned.

        include_extra_coords: `bool`
           If True, also search extra_coords for axis name.
           Default=True.

        Returns
        -------
        axes: `tuple` of `int`
            The pixel axis numbers (numpy ordering) that correspond to the supplied
            axis names.
            If more than one axis corresponds to the physical type, that physical type's
            entry in the output will be a tuple of `int`.
            If no axes names supplied, the ordering of the axis numbers returned will
            correspond to the physical types returned by NDCube.world_axis_physical_types.
        """
        # Parse user input.
        if axes_names == ():
            axes_names = np.array(self.world_axis_physical_types)
        elif isinstance(axes_names, str):
            axes_names = np.array([axes_names])
        else:
            axes_names = np.array(axes_names)
        n_names = len(axes_names)
        axes = np.array([None] * n_names, dtype=object)

        # Get world types associated with sequence axis.
        raster_axis_names = [self._raster_axis_physical_type]
        # If desired, also get extra coord sequence axis names.
        if include_extra_coords:
            extra_sequence_names = ndcube.utils.sequence._get_axis_extra_coord_names_and_units(self.data, None)[0]
            if extra_sequence_names is not None:
                raster_axis_names += list(extra_sequence_names)
        raster_axis_names = np.asarray(raster_axis_names)
        # Find indices of axes_names that correspond to sequence axis and
        # and enter axis number to output
        raster_names_indices = np.isin(axes_names, raster_axis_names)
        axes[raster_names_indices] = self._raster_axis_index

        # Get indices of cube axis names and use Raster version of this method to get axis numbers.
        cube_names_indices = np.invert(raster_names_indices)
        if cube_names_indices.any():
            # Get axes from cube.
            cube_axes = self.data[0]._world_types_to_pixel_axes(
                    *axes_names, include_extra_coords=include_extra_coords)
            # Must use loop as can't use array indexing to enter >1 tuples
            # into single array elements.
            for i in np.arange(len(cube_names_indices))[cube_names_indices]:
                if isinstance(cube_axes[i], tuple):
                   this_axis_indices = tuple(np.array(cube_axes[i]) + 1)
                else:
                    this_axis_indices = cube_axes[i] + 1
                axes[i] = this_axis_indices

        return tuple(axes)

    def _world_types_to_SnS_axes(self, *axes_names, include_extra_coords=True):
        return self.data[0]._world_types_to_pixel_axes(
                *axes_names, include_extra_coords=include_extra_coords)


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
        return sunraster.utils.sequence._slice_sequence_as_SnS(self.seq, item)


class _SequenceSlicer:
    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):
        return ndcube.utils.sequence.slice_sequence(self.seq, item)
