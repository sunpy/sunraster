import abc
import textwrap
from copy import deepcopy

import numpy as np

import astropy.units as u
from astropy.time import Time

import ndcube.utils.wcs as nuw
from ndcube import NDMeta
from ndcube.ndcube import NDCube

__all__ = ["SpectrogramABC","SpectrogramCube"]


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

    @property
    @abc.abstractmethod
    def spectral_axis(self):
        """
        Return the spectral coordinates for each pixel.
        """

    @property
    @abc.abstractmethod
    def time(self):
        """
        Return the time coordinates for each pixel.
        """

    @property
    @abc.abstractmethod
    def exposure_time(self):
        """
        Return the exposure time for each exposure.
        """

    @property
    @abc.abstractmethod
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
        result: `sunraster.SpectrogramCube`
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
    wcs: `astropy.wcs.WCS`
        The WCS object containing the axes' information
    unit : `astropy.units.Unit` or `str`, optional
        Unit for the dataset. Strings that can be converted to a Unit are allowed.
    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn't mandatory. If the uncertainty
        has no such attribute the uncertainty is stored as UnknownUncertainty.
        Defaults to None.
    meta : `dict` object, optional
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

    Attributes
    ----------
    array_axis_physical_types
    axis_world_coords
    dimensions
    extra_coords
    meta
    pixel_to_world
    plot
    spectral
    uncertainty
    world_to_pixel
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
        if instrument_axes is not None:
            if len(instrument_axes) != data.ndim:
                raise ValueError("Length of instrument_axes must match number of data axes.")
            if meta is None:
                meta = NDMeta()
            if not isinstance(meta, NDMeta):
                meta = NDMeta(meta)
            meta.add("instrument_axes", np.asarray(instrument_axes, dtype=str), axes=np.arange(data.ndim, dtype=int), overwrite=True)

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
            component_names = {item: key for key, item in sc.representation_component_names.items()}
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
            if AXIS_NOT_FOUND_ERROR not in err.args[0]:
                raise err
            lon_range = None
            lat_range = None
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
                Pixel dimensions: {self.shape if hasattr(self, "shape") else self.dimensions.astype(int)}
                Longitude range: {lon_range}
                Latitude range: {lat_range}
                Spectral range: {spectral_range}
                Data unit: {self.unit}"""
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{self!s}"

    @property
    def instrument_axes(self):
        """
        The relationship between the array axes and the instrument,
        i.e. repeat raster axis, slit position, position along slit, and spectral.
        """
        return self.meta.get("instrument_axes")

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
            raise ValueError(f"Spectral{AXIS_NOT_FOUND_ERROR}" + f"{SUPPORTED_SPECTRAL_NAMES}")
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
        if coord_loc == "extra_coords":
            return self.axis_world_coords(wcs=self.extra_coords[axis_name])[0]
        if coord_loc == "global_coords":
            return self.global_coords[axis_name]
        if coord_loc == "meta":
            return self.meta[axis_name]
        raise ValueError(f"{coord_loc} is not a valid coordinate location.")

    def _get_axis_coord_index(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            coord_pix_axes = nuw.physical_type_to_pixel_axes(axis_name, self.wcs)
            coord_array_axes = nuw.convert_between_array_and_pixel_axes(coord_pix_axes, len(self.dimensions))
            return coord_array_axes.tolist()[0]
        if coord_loc == "extra_coords":
            return self.extra_coords[axis_name].mapping[0]
        if coord_loc == "meta":
            return self.meta.axes[axis_name]
        raise ValueError(f"{coord_loc} is not a valid coordinate location.")

    def _get_axis_coord_values(self, axis_name, coord_loc):
        if coord_loc == "wcs":
            return self.axis_world_coords_values(axis_name)[0]
        if coord_loc == "extra_coords":
            return self.axis_world_coords_values(wcs=self.extra_coords[axis_name])[0]
        if coord_loc == "global_coords":
            return self.global_coords[axis_name]
        raise ValueError(f"{coord_loc} is not a valid coordinate location.")


def _find_axis_name(supported_names, world_axis_physical_types, extra_coords, meta):
    """
    Finds name of a SpectrogramCube axis type from WCS and extra coords.

    Parameters
    ----------
    supported_names: 1D `numpy.ndarray`
        The names for the axis supported by `sunraster.SpectrogramCube`.
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
        name_index = np.arange(len(names_array))[name_index][0]
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
    old_unit: `astropy.units.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data: `numpy.ndarray`
        Data array with exposure time corrected for.
    new_uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in new_data.
    new_unit: `astropy.units.Unit`
        Unit of new_data array after exposure time correction.
    """
    if force is not True and u.s in unit.decompose().bases:
        raise ValueError(APPLY_EXPOSURE_TIME_ERROR)
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
    old_unit: `astropy.units.Unit`
        Unit of data arrays.
    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data: `numpy.ndarray`
        Data array with exposure time corrected for.
    new_uncertainty: `astropy.nddata.nduncertainty.NDUncertainty`
        The uncertainty of each element in new_data.
    new_unit: `astropy.units.Unit`
        Unit of new_data array after exposure time correction.
    """
    if force is not True and u.s in (unit * u.s).decompose().bases:
        raise ValueError(UNDO_EXPOSURE_TIME_ERROR)
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
