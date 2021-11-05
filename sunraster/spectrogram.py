import abc
import numbers
import textwrap
from copy import deepcopy

import astropy.units as u
import ndcube.utils.wcs as nuw
import numpy as np
from astropy.time import Time
from ndcube.ndcube import NDCube

from sunraster.meta import Meta

__all__ = ["SpectrogramCube"]


# Define some custom error messages.
APPLY_EXPOSURE_TIME_ERROR = (
    "Exposure time correction has probably already "
    "been applied since the unit already includes "
    "inverse time. To apply exposure time correction "
    "anyway, set 'force' kwarg to True."
)
UNDO_EXPOSURE_TIME_ERROR = (
    "Exposure time correction has probably already "
    "been undone since the unit does not include "
    "inverse time. To undo exposure time correction "
    "anyway, set 'force' kwarg to True."
)
AXIS_NOT_FOUND_ERROR = "axis not found. If in extra_coords, axis name must be supported:"

# Define supported coordinate names for coordinate properties.
SUPPORTED_LONGITUDE_NAMES = [
    "custom:pos.helioprojective.lon",
    "pos.helioprojective.lon",
    "longitude",
    "lon",
]
SUPPORTED_LONGITUDE_NAMES += [name.upper() for name in SUPPORTED_LONGITUDE_NAMES] + [
    name.capitalize() for name in SUPPORTED_LONGITUDE_NAMES
]
SUPPORTED_LONGITUDE_NAMES = np.array(SUPPORTED_LONGITUDE_NAMES)

SUPPORTED_LATITUDE_NAMES = [
    "custom:pos.helioprojective.lat",
    "pos.helioprojective.lat",
    "latitude",
    "lat",
]
SUPPORTED_LATITUDE_NAMES += [name.upper() for name in SUPPORTED_LATITUDE_NAMES] + [
    name.capitalize() for name in SUPPORTED_LATITUDE_NAMES
]
SUPPORTED_LATITUDE_NAMES = np.array(SUPPORTED_LATITUDE_NAMES)

SUPPORTED_SPECTRAL_NAMES = [
    "em.wl",
    "em.energy",
    "em.freq",
    "wavelength",
    "energy",
    "frequency",
    "freq",
    "lambda",
    "spectral",
]
SUPPORTED_SPECTRAL_NAMES += [name.upper() for name in SUPPORTED_SPECTRAL_NAMES] + [
    name.capitalize() for name in SUPPORTED_SPECTRAL_NAMES
]
SUPPORTED_SPECTRAL_NAMES = np.array(SUPPORTED_SPECTRAL_NAMES)

SUPPORTED_TIME_NAMES = ["time"]
SUPPORTED_TIME_NAMES += [name.upper() for name in SUPPORTED_TIME_NAMES] + [
    name.capitalize() for name in SUPPORTED_TIME_NAMES
]
SUPPORTED_TIME_NAMES = np.array(SUPPORTED_TIME_NAMES)

SUPPORTED_EXPOSURE_NAMES = [
    "exposure time",
    "exposure_time",
    "exposure times",
    "exposure_times",
    "exp time",
    "exp_time",
    "exp times",
    "exp_times",
]
SUPPORTED_EXPOSURE_NAMES += [name.upper() for name in SUPPORTED_EXPOSURE_NAMES] + [
    name.capitalize() for name in SUPPORTED_EXPOSURE_NAMES
]
SUPPORTED_EXPOSURE_NAMES = np.array(SUPPORTED_EXPOSURE_NAMES)


class SpectrogramABC(abc.ABC):

    # Abstract Base Class to define the basic API of Spectrogram classes.

    @abc.abstractproperty
    def spectral_axis(self):
        """
        Return the spectral coordinates for each pixel.
        """

    @abc.abstractproperty
    def time(self):
        """
        Return the time coordinates for each pixel.
        """

    @abc.abstractproperty
    def exposure_time(self):
        """
        Return the exposure time for each exposure.
        """

    @abc.abstractproperty
    def celestial(self):
        """
        Return the celestial coordinates for each pixel.
        """

    @abc.abstractmethod
    def apply_exposure_time_correction(self, undo=False, force=False):
        """
        Applies or undoes exposure time correction to data and uncertainty and
        adjusts unit.

        Correction is only applied (undone) if the object's unit doesn't (does)
        already include inverse time. This can be overridden so that correction
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
            If True, correction is applied (undone) regardless of unit. Unit is still
            adjusted accordingly.

        Returns
        -------
        result: `SpectrogramCube`
            New SpectrogramCube in new units.
        """


class SpectrogramCube(NDCube, SpectrogramABC):
    """
    Class representing a sit-and-stare or single raster of slit spectrogram
    data.

    Must be described by a single WCS.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.
    wcs: `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information
    unit : `astropy.unit.Unit` or `str`, optional
        Unit for the dataset. Strings that can be converted to a Unit are allowed.
    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn't mandatory. If the uncertainty
        has no such attribute the uncertainty is stored as UnknownUncertainty.
        Defaults to None.
    meta : dict-like object, optional
        Additional meta information about the dataset.
    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.
    instrument_axes : list, optional
        This is the relationship between the array axes and the instrument,
        i.e. repeat raster axis, slit position, position along slit, and spectral.
        These are needed because they cannot be inferred simply from the physical types.
    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every attribute
        before saving it while False tries to save every parameter as reference.
        Note however that it is not always possible to save the input as reference.
        Default is False.
    """

    def __init__(
        self,
        data,
        wcs,
        unit=None,
        uncertainty=None,
        meta=None,
        mask=None,
        instrument_axes=None,
        copy=False,
        **kwargs,
    ):
        # Initialize SpectrogramCube.
        super().__init__(
            data,
            wcs=wcs,
            uncertainty=uncertainty,
            mask=mask,
            meta=meta,
            unit=unit,
            copy=copy,
            **kwargs,
        )
        # Determine labels and location of each key real world coordinate.
        self_extra_coords = self.extra_coords
        world_axis_physical_types = np.array(self.wcs.world_axis_physical_types)
        self._longitude_name, self._longitude_loc = _find_axis_name(
            SUPPORTED_LONGITUDE_NAMES,
            world_axis_physical_types,
            self_extra_coords,
            self.meta,
        )
        self._latitude_name, self._latitude_loc = _find_axis_name(
            SUPPORTED_LATITUDE_NAMES,
            world_axis_physical_types,
            self_extra_coords,
            self.meta,
        )
        self._spectral_name, self._spectral_loc = _find_axis_name(
            SUPPORTED_SPECTRAL_NAMES,
            world_axis_physical_types,
            self_extra_coords,
            self.meta,
        )
        self._time_name, self._time_loc = _find_axis_name(
            SUPPORTED_TIME_NAMES,
            world_axis_physical_types,
            self_extra_coords,
            self.meta,
        )
        self._exposure_time_name, self._exposure_time_loc = _find_axis_name(
            SUPPORTED_EXPOSURE_NAMES,
            world_axis_physical_types,
            self_extra_coords,
            self.meta,
        )
        # Set up instrument axes if set.
        if instrument_axes is None:
            self.instrument_axes = instrument_axes
        elif len(instrument_axes) != data.ndim:
            raise ValueError("Length of instrument_axes must match number of data axes.")
        else:
            self.instrument_axes = np.asarray(instrument_axes, dtype=str)
        # TODO: Remove after ndcube 2.1
        if not isinstance(self.meta, Meta):
            self.meta = Meta(self.meta, data_shape=self.data.shape)

    def __str__(self):
        try:
            if self.time.isscalar:
                time_period = self.time
            else:
                times = self.time
                time_period = Time([times.min(), times.max()]).iso
        except ValueError as err:
            if AXIS_NOT_FOUND_ERROR in err.args[0]:
                time_period = None
            else:
                raise err
        try:
            sc = self.celestial
            component_names = dict([(item, key) for key, item in sc.representation_component_names.items()])
            lon = getattr(sc, component_names["lon"])
            lat = getattr(sc, component_names["lat"])
            if sc.isscalar:
                lon_range = lon
                lat_range = lat
            elif sc.size == 0:
                lon_range = None
                lat_range = None
            else:
                lon_range = u.Quantity([lon.min(), lon.max()])
                lat_range = u.Quantity([lat.min(), lat.max()])
        except ValueError as err:
            if AXIS_NOT_FOUND_ERROR in err.args[0]:
                lon_range = None
                lat_range = None
            else:
                raise err
        try:
            if self.spectral_axis.isscalar:
                spectral_range = self.spectral_axis
            else:
                spectral_range = u.Quantity([self.spectral_axis.min(), self.spectral_axis.max()])
        except ValueError as err:
            if AXIS_NOT_FOUND_ERROR in err.args[0]:
                spectral_range = None
            else:
                raise err
        return textwrap.dedent(
            f"""\
                {self.__class__.__name__}
                {"".join(["-"] * len(self.__class__.__name__))}
                Time Period: {time_period}
                Instrument axes: {self.instrument_axes}
                Pixel dimensions: {self.dimensions.astype(int)}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.unit}"""
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __getitem__(self, item):
        # Slice SpectrogramCube using parent slicing.
        result = super().__getitem__(item)
        # Slice instrument_axes if it exists.
        # If item is a slice, cube dimensionality is not reduced
        # so instrument_axes need not be sliced.
        if self.instrument_axes is None or isinstance(item, slice):
            instrument_axes = self.instrument_axes
        else:
            # If item is int, instrument_axes needs slicing.
            if isinstance(item, numbers.Integral):
                instrument_axes = self.instrument_axes[1:]
            # If item is tuple, instrument axes will need to be sliced if tuple contains an int.
            elif isinstance(item, tuple):
                instr_item = [isinstance(i, numbers.Integral) for i in item] + [False] * (
                    len(self.instrument_axes) - len(item)
                )
                instrument_axes = self.instrument_axes[np.invert(instr_item)]
            else:
                raise TypeError("Unrecognized slice item. Must be int, slice or tuple.")
            # If slicing causes cube to be a scalar, set instrument_axes to None.
            if len(instrument_axes) == 0:
                instrument_axes = None
        result.instrument_axes = instrument_axes
        # TODO: Remove for ndcube 2.1
        # Slice metadata if possible.
        try:
            result.meta = self.meta[item]
        except TypeError as err:
            if "unhashable type" not in err.args[0]:
                raise err
        except KeyError:
            pass
        return result

    @property
    def spectral_axis(self):
        if not self._spectral_name:
            self._spectral_name, self._spectral_loc = _find_axis_name(
                SUPPORTED_SPECTRAL_NAMES,
                self.wcs.world_axis_physical_types,
                self.extra_coords,
                self.meta,
            )
        if not self._spectral_name:
            raise ValueError("Spectral" + AXIS_NOT_FOUND_ERROR + f"{SUPPORTED_SPECTRAL_NAMES}")
        return self._get_axis_coord(self._spectral_name, self._spectral_loc)

    @property
    def time(self):
        if not self._time_name:
            self._time_name, self._time_loc = _find_axis_name(
                SUPPORTED_TIME_NAMES,
                self.wcs.world_axis_physical_types,
                self.extra_coords,
                self.meta,
            )
            if not self._time_name:
                raise ValueError(f"Time {AXIS_NOT_FOUND_ERROR} {SUPPORTED_TIME_NAMES}")
        return Time(self._get_axis_coord(self._time_name, self._time_loc))

    @property
    def exposure_time(self):
        if not self._exposure_time_name or not hasattr(self, "_exposure_time_loc"):
            self._exposure_time_name, self._exposure_time_loc = _find_axis_name(
                SUPPORTED_EXPOSURE_NAMES,
                self.wcs.world_axis_physical_types,
                self.extra_coords,
                self.meta,
            )
            if not self._exposure_time_name:
                raise ValueError(f"Exposure time {AXIS_NOT_FOUND_ERROR} {SUPPORTED_EXPOSURE_NAMES}")
        return self._get_axis_coord(self._exposure_time_name, self._exposure_time_loc)

    @property
    def celestial(self):
        if not self._longitude_name:
            self._longitude_name, self._longitude_loc = _find_axis_name(
                SUPPORTED_LONGITUDE_NAMES,
                self.wcs.world_axis_physical_types,
                self.extra_coords,
                self.meta,
            )
        if not self._latitude_name:
            self._latitude_name, self._latitude_loc = _find_axis_name(
                SUPPORTED_LATITUDE_NAMES,
                self.wcs.world_axis_physical_types,
                self.extra_coords,
                self.meta,
            )
        if self._longitude_name:
            celestial_name = self._longitude_name
            celestial_loc = self._longitude_loc
        elif self._latitude_name:
            celestial_name = self._latitude_name
            celestial_loc = self._latitude_loc
        else:
            raise ValueError(
                f"Celestial {AXIS_NOT_FOUND_ERROR} "
                f"{np.concatenate([SUPPORTED_LONGITUDE_NAMES, SUPPORTED_LATITUDE_NAMES])}"
            )
        return self._get_axis_coord(celestial_name, celestial_loc)

    def apply_exposure_time_correction(self, undo=False, force=False):
        # Get exposure time in seconds.
        exposure_time_s = self.exposure_time.to(u.s).value
        # If exposure time is not scalar, change array's shape so that
        # it can be broadcast with data and uncertainty arrays.
        if not np.isscalar(exposure_time_s):
            (exposure_axis,) = self._get_axis_coord_index(self._exposure_time_name, self._exposure_time_loc)
            # Change array shape for broadcasting
            item = [np.newaxis] * self.data.ndim
            item[exposure_axis] = slice(None)
            exposure_time_s = exposure_time_s[tuple(item)]
        # Based on value on undo kwarg, apply or remove exposure time correction.
        if undo is True:
            new_data, new_uncertainty, new_unit = _uncalculate_exposure_time_correction(
                self.data, self.uncertainty, self.unit, exposure_time_s, force=force
            )
        else:
            new_data, new_uncertainty, new_unit = _calculate_exposure_time_correction(
                self.data, self.uncertainty, self.unit, exposure_time_s, force=force
            )
        # Return new instance of SpectrogramCube with correction applied/undone.
        new_cube = deepcopy(self)
        new_cube._data = new_data
        new_cube._uncertainty = new_uncertainty
        new_cube._extra_coords = self.extra_coords
        new_cube._unit = new_unit
        return new_cube

    def _get_axis_coord(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            return self.axis_world_coords(axis_name)[0]
        elif coord_loc == "extra_coords":
            return self.axis_world_coords(wcs=self.extra_coords[axis_name])[0]
        elif coord_loc == "global_coords":
            return self.global_coords[axis_name]
        elif coord_loc == "meta":
            return self.meta[axis_name]
        else:
            raise ValueError(f"{coord_loc} is not a valid coordinate location.")

    def _get_axis_coord_index(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            coord_pix_axes = nuw.physical_type_to_pixel_axes(axis_name, self.wcs)
            coord_array_axes = nuw.convert_between_array_and_pixel_axes(coord_pix_axes, len(self.dimensions))
            return coord_array_axes.tolist()[0]
        elif coord_loc == "extra_coords":
            return self.extra_coords[axis_name].mapping[0]
        elif coord_loc == "meta":
            return self.meta.axes[axis_name]
        else:
            raise ValueError(f"{coord_loc} is not a valid coordinate location.")

    def _get_axis_coord_values(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            return self.axis_world_coords_values(axis_name)[0]
        elif coord_loc == "extra_coords":
            return self.axis_world_coords_values(wcs=self.extra_coords[axis_name])[0]
        elif coord_loc == "global_coords":
            return self.global_coords[axis_name]
        else:
            raise ValueError(f"{coord_loc} is not a valid coordinate location.")


def _find_axis_name(supported_names, world_axis_physical_types, extra_coords, meta):
    """
    Finds name of a SpectrogramCube axis type from WCS and extra coords.

    Parameters
    ----------
    supported_names: 1D `numpy.ndarray`
        The names for the axis supported by `SpectrogramCube`.
    world_axis_physical_types: 1D `numpy.ndarray`
        Output of SpectrogramCube.world_axis_physical_types converted to an array.
    extra_coords: `ndcube.ExtraCoords` or `None`
        Output of SpectrogramCube.extra_coords
    meta: Meta or `None`
        Output of SpectrogramCube.meta

    Returns
    -------
    axis_name: `str`
        The coordinate name of the axis.
    loc: `str`
        The location where the coordinate is stored: "wcs" or "extra_coords".
    """
    axis_name = None
    loc = None
    # Check WCS for axis name.
    axis_name = _find_name_in_array(supported_names, world_axis_physical_types)
    if axis_name:
        loc = "wcs"
    elif extra_coords:  # If axis name not in WCS, check extra_coords.
        axis_name = _find_name_in_array(supported_names, np.array(list(extra_coords.keys())))
        if axis_name:
            loc = "extra_coords"
    if axis_name is None and meta:  # If axis name not in WCS, check meta.
        axis_name = _find_name_in_array(supported_names, np.array(list(meta.keys())))
        if axis_name:
            loc = "meta"
    return axis_name, loc


def _find_name_in_array(supported_names, names_array):
    name_index = np.isin(names_array, supported_names)
    if name_index.sum() > 0:
        name_index = int(np.arange(len(names_array))[name_index])
        return names_array[name_index]
    return None


def _calculate_exposure_time_correction(data, uncertainty, unit, exposure_time, force=False):
    """
    Applies exposure time correction to data and uncertainty arrays.

    Parameters
    ----------
    data: `numpy.ndarray`
        Data array to be converted.
    uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in data.
    old_unit: `astropy.unit.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data: `numpy.ndarray`
        Data array with exposure time corrected for.
    new_uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in new_data.
    new_unit: `astropy.unit.Unit`
        Unit of new_data array after exposure time correction.
    """

    # If force is not set to True and unit already includes inverse time,
    # raise error as exposure time correction has probably already been
    # applied and should not be applied again.
    if force is not True and u.s in unit.decompose().bases:
        raise ValueError(APPLY_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does not include inverse time and so
        # exposure does need to be applied, or
        # user has set force=True and wants the correction applied
        # regardless of the unit.
        new_data = data / exposure_time
        if uncertainty:
            uncertainty_unit = uncertainty.unit / u.s if uncertainty.unit else uncertainty.unit
            new_uncertainty = uncertainty.__class__(uncertainty.array / exposure_time, unit=uncertainty_unit)
        else:
            new_uncertainty = uncertainty
        new_unit = unit / u.s
    return new_data, new_uncertainty, new_unit


def _uncalculate_exposure_time_correction(data, uncertainty, unit, exposure_time, force=False):
    """
    Removes exposure time correction from data and uncertainty arrays.

    Parameters
    ----------
    data: `numpy.ndarray`
        Data array to be converted.
    uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in data.
    old_unit: `astropy.unit.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data: `numpy.ndarray`
        Data array with exposure time corrected for.
    new_uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in new_data.
    new_unit: `astropy.unit.Unit`
        Unit of new_data array after exposure time correction.
    """

    # If force is not set to True and unit does not include inverse time,
    # raise error as exposure time correction has probably already been
    # undone and should not be undone again.
    if force is not True and u.s in (unit * u.s).decompose().bases:
        raise ValueError(UNDO_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does include inverse time and so
        # exposure does need to be removed, or
        # user has set force=True and wants the correction removed
        # regardless of the unit.
        new_data = data * exposure_time
        if uncertainty:
            uncertainty_unit = uncertainty.unit * u.s if uncertainty.unit else uncertainty.unit
            new_uncertainty = uncertainty.__class__(uncertainty.array * exposure_time, unit=uncertainty_unit)
        else:
            new_uncertainty = uncertainty
        new_unit = unit * u.s
    return new_data, new_uncertainty, new_unit
