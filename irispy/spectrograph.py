# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

import copy
import datetime

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from ndcube import NDCube, NDCubeSequence, _extra_coords_to_input_format
from ndcube.wcs_util import WCS
from sunpy.time import parse_time

from irispy.spectrogram import SpectrogramSequence
from irispy import iris_tools

__all__ = ['IRISSpectrograph']

class IRISSpectrograph(object):
    """
    An object to hold data from multiple IRIS raster scans.

    The Interface Region Imaging Spectrograph (IRIS) small explorer spacecraft
    provides simultaneous spectra and images of the photosphere, chromosphere,
    transition region, and corona with 0.33 to 0.4 arcsec spatial resolution,
    2-second temporal resolution and 1 km/s velocity resolution over a
    field-of- view of up to 175 arcsec by 175 arcsec.  IRIS consists of a 19-cm
    UV telescope that feeds a slit-based dual-bandpass imaging spectrograph.

    IRIS was launched into a Sun-synchronous orbit on 27 June 2013.

    References
    ----------
    * `IRIS Mission Page <http://iris.lmsal.com>`_
    * `IRIS Analysis Guide <https://iris.lmsal.com/itn26/itn26.pdf>`_
    * `IRIS Instrument Paper <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0196&file_type=pdf>`_
    * `IRIS FITS Header keywords <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0077&file_type=pdf>`_
    """

    def __init__(self, filenames, spectral_windows="All", common_axis=0):
        """Initializes an IRISSpectrograph object from IRIS level 2 files."""
        # default common axis is 0.
        if type(filenames) is str:
            filenames = [filenames]
        for f, filename in enumerate(filenames):
            hdulist = fits.open(filename)
            hdulist.verify('fix')
            if f == 0:
                # Determine number of raster positions in a scan
                raster_positions_per_scan = int(hdulist[0].header["NRASTERP"])
                # collecting the window observations.
                windows_in_obs = np.array([hdulist[0].header["TDESC{0}".format(i)]
                                           for i in range(1, hdulist[0].header["NWIN"]+1)])
                # if spectral window is All then get every window
                # else take the appropriate windows
                if spectral_windows == "All":
                    spectral_windows_req = windows_in_obs
                    window_fits_indices = range(1, len(hdulist)-2)
                else:
                    if type(spectral_windows) is str:
                        spectral_windows_req = [spectral_windows]
                    spectral_windows_req = np.asarray(spectral_windows_req, dtype="U")
                    window_is_in_obs = np.asarray(
                        [window in windows_in_obs for window in spectral_windows_req])
                    if not all(window_is_in_obs):
                        missing_windows = window_is_in_obs == False
                        raise ValueError("Spectral windows {0} not in file {1}".format(
                            spectral_windows[missing_windows], filenames[0]))
                    window_fits_indices = np.nonzero(np.in1d(windows_in_obs,
                                                             spectral_windows))[0]+1
                # Generate top level meta dictionary from first file
                # main header.
                self.meta = {"TELESCOP": hdulist[0].header["TELESCOP"],
                             "INSTRUME": hdulist[0].header["INSTRUME"],
                             "DATA_LEV": hdulist[0].header["DATA_LEV"],
                             "OBSID": hdulist[0].header["OBSID"],
                             "OBS_DESC": hdulist[0].header["OBS_DESC"],
                             "STARTOBS": parse_time(hdulist[0].header["STARTOBS"]),
                             "ENDOBS": parse_time(hdulist[0].header["ENDOBS"]),
                             "SAT_ROT": hdulist[0].header["SAT_ROT"] * u.deg,
                             "AECNOBS": int(hdulist[0].header["AECNOBS"]),
                             "FOVX": hdulist[0].header["FOVX"] * u.arcsec,
                             "FOVY": hdulist[0].header["FOVY"] * u.arcsec,
                             "SUMSPTRN": hdulist[0].header["SUMSPTRN"],
                             "SUMSPTRF": hdulist[0].header["SUMSPTRF"],
                             "SUMSPAT": hdulist[0].header["SUMSPAT"],
                             "NEXPOBS": hdulist[0].header["NEXPOBS"],
                             "NRASTERP": hdulist[0].header["NRASTERP"],
                             "KEYWDDOC": hdulist[0].header["KEYWDDOC"]}
                # Initialize meta dictionary for each spectral_window
                window_metas = {}
                for i, window_name in enumerate(spectral_windows_req):
                    if "FUV" in hdulist[0].header["TDET{0}".format(window_fits_indices[i])]:
                        spectral_summing = hdulist[0].header["SUMSPTRF"]
                    else:
                        spectral_summing = hdulist[0].header["SUMSPTRN"]
                    window_metas[window_name] = {
                        "detector type":
                            hdulist[0].header["TDET{0}".format(window_fits_indices[i])],
                        "spectral window":
                            hdulist[0].header["TDESC{0}".format(window_fits_indices[i])],
                        "brightest wavelength":
                            hdulist[0].header["TWAVE{0}".format(window_fits_indices[i])],
                        "min wavelength":
                            hdulist[0].header["TWMIN{0}".format(window_fits_indices[i])],
                        "max wavelength":
                            hdulist[0].header["TWMAX{0}".format(window_fits_indices[i])],
                        "SAT_ROT": hdulist[0].header["SAT_ROT"],
                        "spatial summing": hdulist[0].header["SUMSPAT"],
                        "spectral summing": spectral_summing
                    }
                # creating a empty list for every spectral window and each spectral window
                # is a key for the dictionary.
                data_dict = dict([(window_name, list())
                                  for window_name in spectral_windows_req])
            # Determine extra coords for this raster.
            times = np.array(
                [parse_time(hdulist[0].header["STARTOBS"]) + datetime.timedelta(seconds=s)
                 for s in hdulist[-2].data[:,hdulist[-2].header["TIME"]]])
            raster_positions = np.arange(int(hdulist[0].header["NRASTERP"]))
            pztx = hdulist[-2].data[:, hdulist[-2].header["PZTX"]] * u.arcsec
            pzty = hdulist[-2].data[:, hdulist[-2].header["PZTY"]] * u.arcsec
            xcenix = hdulist[-2].data[:, hdulist[-2].header["XCENIX"]] * u.arcsec
            ycenix = hdulist[-2].data[:, hdulist[-2].header["YCENIX"]] * u.arcsec
            obs_vrix = hdulist[-2].data[:, hdulist[-2].header["OBS_VRIX"]] * u.m/u.s
            ophaseix = hdulist[-2].data[:, hdulist[-2].header["OPHASEIX"]]
            exposure_times_fuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEF"]] * u.s
            exposure_times_nuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEN"]] * u.s
            general_extra_coords = [("time", 0, times), ("raster position", 0, raster_positions),
                                    ("pztx", 0, pztx), ("pzty", 0, pzty),
                                    ("xcenix", 0, xcenix), ("ycenix", 0, ycenix),
                                    ("obs_vrix", 0, obs_vrix), ("ophaseix", 0, ophaseix)]
            for i, window_name in enumerate(spectral_windows_req):
                # Determine values of properties dependent on detector type.
                if "FUV" in hdulist[0].header["TDET{0}".format(window_fits_indices[i])]:
                    exposure_times = exposure_times_fuv
                    DN_unit = iris_tools.DN_UNIT["FUV"]
                    readout_noise = iris_tools.READOUT_NOISE["FUV"]
                elif "NUV" in hdulist[0].header["TDET{0}".format(window_fits_indices[i])]:
                    exposure_times = exposure_times_nuv
                    DN_unit = iris_tools.DN_UNIT["NUV"]
                    readout_noise = iris_tools.READOUT_NOISE["NUV"]
                else:
                    raise ValueError("Detector type in FITS header not recognized.")
                # Derive WCS, data and mask for NDCube from file.
                wcs_ = WCS(hdulist[window_fits_indices[i]].header)
                data_nan_masked = copy.deepcopy(hdulist[window_fits_indices[i]].data)
                data_nan_masked[hdulist[window_fits_indices[i]].data == -200.] = np.nan
                data_mask = hdulist[window_fits_indices[i]].data == -200.
                # Derive extra coords for this spectral window.
                window_extra_coords = copy.deepcopy(general_extra_coords)
                window_extra_coords.append(("exposure time", 0, exposure_times))
                # Collect metadata relevant to single files.
                meta = {"SAT_ROT": hdulist[0].header["SAT_ROT"] * u.deg,
                        "DATE_OBS": parse_time(hdulist[0].header["DATE_OBS"]),
                        "DATE_END": parse_time(hdulist[0].header["DATE_END"]),
                        "HLZ": bool(int(hdulist[0].header["HLZ"])),
                        "SAA": bool(int(hdulist[0].header["SAA"])),
                        "DSUN_OBS": hdulist[0].header["DSUN_OBS"] * u.m,
                        "IAECEVFL": hdulist[0].header["IAECEVFL"],
                        "IAECFLAG": hdulist[0].header["IAECFLAG"],
                        "IAECFLFL": hdulist[0].header["IAECFLFL"],
                        "KEYWDDOC": hdulist[0].header["KEYWDDOC"],
                        "detector type":
                            hdulist[0].header["TDET{0}".format(window_fits_indices[i])],
                        "spectral window": window_name,
                        "OBSID": hdulist[0].header["OBSID"]
                        }
                # Derive uncertainty of data
                uncertainty = u.Quantity(np.sqrt(
                    (data_nan_masked*DN_unit).to(u.photon).value + readout_noise.to(u.photon).value**2),
                    unit=u.photon).to(DN_unit).value
                # Appending NDCube instance to the corresponding window key in dictionary's list.
                data_dict[window_name].append(
                    IRISSpectrogram(data_nan_masked, wcs_, uncertainty, DN_unit, meta,
                                    window_extra_coords, mask=data_mask))
            hdulist.close()
        # Attach dictionary containing level 1 and wcs info for each file used.
        # making a NDCubeSequence of every dictionary key window.
        self.data = dict([(window_name,
                           IRISSpectrogramSequence(data_dict[window_name], common_axis,
                                               raster_positions_per_scan,
                                               first_exposure_raster_position=0,
                                               meta=window_metas[window_name]))
                          for window_name in spectral_windows_req])

    def __repr__(self):
        spectral_window = self.spectral_windows["spectral window"][0]
        spectral_windows_info = "".join(
            ["\n    {0}\n        (raster axis, slit axis, spectral axis) {1}".format(
                name,
                self.data[name].dimensions[1::])
                for name in self.spectral_windows["spectral window"]])
        return "<iris.IRISSpectrograph instance\nOBS ID: {0}\n".format(self.meta["OBSID"]) + \
               "OBS Description: {0}\n".format(self.meta["OBS_DESC"]) + \
               "OBS period: {0} -- {1}\n".format(self.meta["STARTOBS"], self.meta["ENDOBS"]) + \
               "Instance period: {0} -- {1}\n".format(
                   self.data[spectral_window][0]._extra_coords["time"]["value"],
                   self.data[spectral_window][-1]._extra_coords["time"]["value"]) + \
               "Number unique raster positions: {0}\n".format(self.meta["NRASTERP"]) + \
               "Spectral windows{0}>".format(spectral_windows_info)

    @property
    def spectral_windows(self):
        """Returns a table of info on the spectral windows."""
        colnames = ("spectral window", "detector type", "brightest wavelength", "min wavelength",
                    "max wavelength")
        spectral_window_list = []
        for key in list(self.data.keys()):
            if type(self.data[key]) == IRISSpectrogramSequence:
                spectral_window_list.append([self.data[key].meta[colname] for colname in colnames])
        return Table(rows=spectral_window_list, names=colnames)

class IRISSpectrogramSequence(SpectrogramSequence):
    """Class for holding, slicing and plotting IRIS spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `IRISSpectrogram` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.

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
                 first_exposure_raster_position, meta):
        detector_type_key = "detector type"
        # Check that meta contains required keys.
        required_meta_keys = [detector_type_key, "spectral window",
                              "brightest wavelength", "min wavelength", "max wavelength"]
        if not all([key in list(meta) for key in required_meta_keys]):
            raise ValueError("Meta must contain following keys: {0}".format(required_meta_keys))
        # Check that all spectrograms are from same specral window and OBS ID.
        if len(np.unique([cube.meta["OBSID"] for cube in data_list])) != 1:
            raise ValueError("Constituent IRISSpectrogram objects must have same "
                             "value of 'OBSID' in its meta.")
        if len(np.unique([cube.meta["spectral window"] for cube in data_list])) != 1:
            raise ValueError("Constituent IRISSpectrogram objects must have same "
                             "value of 'spectral window' in its meta.")
        # Initialize Sequence.
        super(IRISSpectrogramSequence, self).__init__(
            data_list, common_axis, raster_positions_per_scan,
            first_exposure_raster_position, meta=meta)

    def __repr__(self):
        number_of_rasters = int(
            (self.dimensions.shape[0] + self.first_exposure_raster_position) / \
                self.raster_positions_per_scan)
        return """IRISSpectrogramSequence
---------------------
{obs_repr}

Sequence period: {inst_start} -- {inst_end}
Rasters:  {n_rasters}
Exposures per Raster: {n_steps}
Axis Types: {axis_types}
Sequence Shape: {seq_shape}

""".format(obs_repr=_produce_obs_repr_string(self.meta),
           inst_start=self[0]._extra_coords["time"]["value"][0],
           inst_end=self[-1]._extra_coords["time"]["value"][-1],
           n_rasters=number_of_rasters, n_steps=self.raster_positions_per_scan,
           axis_types=self.dimensions.axis_types[::], seq_shape=self.dimensions.shape)

    def convert_to(self, new_unit_type, copy=False):
        """
        Converts data, uncertainty and unit of each spectrogram in sequence to new unit.

        Parameters
        ----------
        new_unit_type: `str`
           Unit type to convert data to.  Three values are accepted:
           "DN": Relevant IRIS data number based on detector type.
           "photons": photon counts
           "radiance": Perorms radiometric calibration conversion.

        copy: `bool`
            If True a new instance with the converted data values is return.
            If False, the current instance is overwritten.
            Default=False

        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.convert_to(new_unit_type))
        if copy is True:
            return IRISSpectrogramSequence(
                converted_data_list, self._common_axis, self.raster_positions_per_scan,
                self.first_exposure_raster_position, meta=self.meta)
        else:
            self.data = converted_data_list

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
        result: `None` or `IRISSpectrogramSequence`
            If copy=False, the original IRISSpectrogramSequence is modified with the
            exposure time correction applied (undone).
            If copy=True, a new IRISSpectrogramSequence is returned with the correction
            applied (undone).

        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.apply_exposure_time_correction(undo=undo,
                                                                           force=force))
        if copy is True:
            return IRISSpectrogramSequence(
                converted_data_list, self._common_axis, self.raster_positions_per_scan,
                self.first_exposure_raster_position, meta=self.meta)
        else:
            self.data = converted_data_list

class IRISSpectrogram(NDCube):
    """
    Class representing IRISSpectrogram data described by a single WCS.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.

    wcs: `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information

    unit : `astropy.unit.Unit` or `str`
        Unit for the dataset. Strings that can be converted to a Unit are allowed.

    meta : dict-like object
        Additional meta information about the dataset. Must contain at least the
        following keys:
            detector type: str, (FUV1, FUV2 or NUV)
            OBSID: int
            spectral window: str

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

    def __init__(self, data, wcs, uncertainty, unit, meta, extra_coords,
                 mask=None, copy=False, missing_axis=None):
        # Check required meta data is provided.
        required_meta_keys = ["detector type"]
        if not all([key in list(meta) for key in required_meta_keys]):
                raise ValueError("Meta must contain following keys: {0}".format(required_meta_keys))
        # Check extra_coords contains required coords.
        required_extra_coords_keys = ["time", "exposure time"]
        extra_coords_keys = [coord[0] for coord in extra_coords]
        if not all([key in extra_coords_keys for key in required_extra_coords_keys]):
            raise ValueError("The following extra coords must be supplied: {0} vs. {1} from {2}".format(
                required_extra_coords_keys, extra_coords_keys, extra_coords))
        # Initialize IRISSpectrogram.
        super(IRISSpectrogram, self).__init__(
            data, wcs, uncertainty=uncertainty, mask=mask, meta=meta,
            unit=unit, extra_coords=extra_coords, copy=copy, missing_axis=missing_axis)

    def __getitem__(self, item):
        result = super(IRISSpectrogram, self).__getitem__(item)
        return IRISSpectrogram(
            result.data, result.wcs, result.uncertainty, result.unit, result.meta,
            _extra_coords_to_input_format(result._extra_coords, result.missing_axis),
            mask=result.mask, missing_axis=result.missing_axis)

    def __repr__(self):
        if self.missing_axis[::-1][self._extra_coords["time"]["axis"]]:
            instance_start = instance_end = self._extra_coords["time"]["value"]
        else:
            instance_start = self._extra_coords["time"]["value"][0],
            instance_end = self._extra_coords["time"]["value"][-1]
        return """IRISSpectrogram
---------------------
{obs_repr}

Spectrogram period: {inst_start} -- {inst_end}
Data shape: {shape}
Axis Types: {axis_types}
""".format(obs_repr=_produce_obs_repr_string(self.meta),
           inst_start=instance_start, inst_end=instance_end,
           shape=self.dimensions, axis_types=self.dimensions.axis_types[::])

    def convert_to(self, new_unit_type):
        """
        Converts data, unit and uncertainty attributes to new unit type.

        The presence or absence of the exposure time correction is
        preserved in the conversions.

        Parameters
        ----------
        new_unit_type: `str`
           Unit type to convert data to.  Three values are accepted:
           "DN": Relevant IRIS data number based on detector type.
           "photons": photon counts
           "radiance": Perorms radiometric calibration conversion.

        Returns
        -------
        result: `IRISSpectrogram`
            New IRISSpectrogram in new units.

        """
        detector_type = iris_tools.get_detector_type(self.meta)
        if new_unit_type == "radiance" or self.unit.is_equivalent(iris_tools.RADIANCE_UNIT):
            # Get spectral dispersion per pixel.
            spectral_wcs_index = np.where(np.array(self.wcs.wcs.ctype) == "WAVE")[0][0]
            spectral_dispersion_per_pixel = self.wcs.wcs.cdelt[spectral_wcs_index] * \
                                            self.wcs.wcs.cunit[spectral_wcs_index]
            # Get solid angle from slit width for a pixel.
            lat_wcs_index = ["HPLT" in c for c in self.wcs.wcs.ctype]
            lat_wcs_index = np.arange(len(self.wcs.wcs.ctype))[lat_wcs_index]
            lat_wcs_index = lat_wcs_index[0]
            solid_angle = self.wcs.wcs.cdelt[lat_wcs_index] * \
                          self.wcs.wcs.cunit[lat_wcs_index] * iris_tools.SLIT_WIDTH
            # Get wavelength for each pixel.
            spectral_data_index = (-1) * (np.arange(len(self.dimensions)) + 1)[spectral_wcs_index]
            obs_wavelength = self.pixel_to_world([
                np.zeros(int(self.dimensions.shape[spectral_data_index].value)) * u.pix,
                np.zeros(int(self.dimensions.shape[spectral_data_index].value)) * u.pix,
                np.arange(int(self.dimensions.shape[spectral_data_index].value)) * u.pix])[-1]
        if new_unit_type == "DN" or new_unit_type == "photons":
            if self.unit.is_equivalent(iris_tools.RADIANCE_UNIT):
                # Convert from radiance to counts/s
                new_data_quantities = iris_tools.convert_or_undo_photons_per_sec_to_radiance(
                    (self.data * self.unit, self.uncertainty.array * self.unit),
                    obs_wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle,
                    undo=True)
                new_data = new_data_quantities[0].value
                new_uncertainty = new_data_quantities[1].value
                new_unit = new_data_quantities[0].unit
                self = IRISSpectrogram(
                    new_data, self.wcs, new_uncertainty, new_unit, self.meta,
                    _extra_coords_to_input_format(self._extra_coords, self.missing_axis),
                    mask=self.mask, missing_axis=self.missing_axis)
            if new_unit_type == "DN":
                new_unit = iris_tools.DN_UNIT[detector_type]
            else:
                new_unit = u.photon
            new_data_arrays, new_unit = iris_tools.convert_between_DN_and_photons(
                (self.data, self.uncertainty.array), self.unit, new_unit)
            new_data = new_data_arrays[0]
            new_uncertainty = new_data_arrays[1]
        elif new_unit_type == "radiance":
            if self.unit.is_equivalent(iris_tools.RADIANCE_UNIT):
                new_data = self.data
                new_uncertainty = self.uncertainty
                new_unit = self.unit
            else:
                # Ensure spectrogram is in units of counts/s.
                cube = self.convert_to("photons")
                try:
                    cube = cube.apply_exposure_time_correction()
                except ValueError(iris_tools.APPLY_EXPOSURE_TIME_ERROR):
                    pass
                # Convert to radiance units.
                new_data_quantities = iris_tools.convert_or_undo_photons_per_sec_to_radiance(
                    (cube.data*cube.unit, cube.uncertainty.array*cube.unit),
                    obs_wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle)
                new_data = new_data_quantities[0].value
                new_uncertainty = new_data_quantities[1].value
                new_unit = new_data_quantities[0].unit
        else:
            raise ValueError("Input unit type not recognized.")
        return IRISSpectrogram(new_data, self.wcs, new_uncertainty, new_unit, self.meta,
                               _extra_coords_to_input_format(self._extra_coords, self.missing_axis),
                               mask=self.mask, missing_axis=self.missing_axis)

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
        result: `IRISSpectrogram`
            New IRISSpectrogram in new units.

        """
        # Get exposure time in seconds and change array's shape so that
        # it can be broadcast with data and uncertainty arrays.
        exposure_time_s = self._extra_coords["exposure time"]["value"].to(u.s).value
        if not self._extra_coords["exposure time"]["value"].isscalar:
            if len(self.dimensions.shape) == 1:
                pass
            elif len(self.dimensions.shape) == 2:
                exposure_time_s = exposure_time_s[:, np.newaxis]
            elif len(self.dimensions.shape) == 3:
                exposure_time_s = exposure_time_s[:, np.newaxis, np.newaxis]
            else:
                raise ValueError(
                    "IRISSpectrogram dimensions must be 2 or 3. Dimensions={0}".format(
                        len(self.dimensions.shape)))
        # Based on value on undo kwarg, apply or remove exposure time correction.
        if undo is True:
            new_data_arrays, new_unit = iris_tools.uncalculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        else:
            new_data_arrays, new_unit = iris_tools.calculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        # Return new instance of IRISSpectrogram with correction applied/undone.
        return IRISSpectrogram(
            new_data_arrays[0], self.wcs, new_data_arrays[1], new_unit, self.meta,
            _extra_coords_to_input_format(self._extra_coords, self.missing_axis),
            mask=self.mask, missing_axis=self.missing_axis)


def _produce_obs_repr_string(meta):
    obs_info = [meta.get(key, "Unknown") for key in ["OBSID", "OBS_DESC", "STARTOBS", "ENDOBS"]]
    return """OBS ID: {obs_id}
OBS Description: {obs_desc}
OBS period: {obs_start} -- {obs_end}""".format(obs_id=obs_info[0], obs_desc=obs_info[1],
                                               obs_start=obs_info[2], obs_end=obs_info[3])