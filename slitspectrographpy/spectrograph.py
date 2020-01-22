# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

import copy

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time, TimeDelta
from ndcube import NDCube, NDCubeSequence
from ndcube.utils.wcs import WCS
from ndcube.utils.cube import convert_extra_coords_dict_to_input_format

__all__ = ['SlitSpectrogramCube', 'SlitSpectrogramCubeSequence']

class SlitSpectrogramCubeSequence(NDCubeSequence):
    """Class for holding, slicing and plotting spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `SlitSpectrogramCube` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.

    meta: `dict` or header object
        Metadata associated with the sequence.

    common_axis: `int`
        The axis of the NDCubes corresponding to time.

    """
    def __init__(self, data_list, meta=None, common_axis=0):
        # Initialize Sequence.
        super(SlitSpectrogramCubeSequence, self).__init__(
            data_list, meta=meta, common_axis=common_axis)

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
        result: `None` or `SlitSpectrogramCubeSequence`
            If copy=False, the original SlitSpectrogramCubeSequence is modified with the
            exposure time correction applied (undone).
            If copy=True, a new SlitSpectrogramCubeSequence is returned with the correction
            applied (undone).

        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.apply_exposure_time_correction(undo=undo,
                                                                           force=force))
        if copy is True:
            return SlitSpectrogramCubeSequence(
                converted_data_list, meta=self.meta, common_axis=self._common_axis)
        else:
            self.data = converted_data_list

class SlitSpectrogramCube(NDCube):
    """
    Class representing SlitSpectrogramCube data described by a single WCS.

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
    def __init__(self, data, wcs, uncertainty, unit, meta, extra_coords,
                 mask=None, copy=False, missing_axes=None):
        # Check extra_coords contains required coords.
        required_extra_coords_keys = ["time", "exposure time"]
        extra_coords_keys = [coord[0] for coord in extra_coords]
        if not all([key in extra_coords_keys for key in required_extra_coords_keys]):
            raise ValueError("The following extra coords must be supplied: {0} vs. {1} from {2}".format(
                required_extra_coords_keys, extra_coords_keys, extra_coords))
        # Initialize SlitSpectrogramCube.
        super(SlitSpectrogramCube, self).__init__(
            data, wcs, uncertainty=uncertainty, mask=mask, meta=meta,
            unit=unit, extra_coords=extra_coords, copy=copy, missing_axes=missing_axes)

    def __getitem__(self, item):
        result = super(SlitSpectrogramCube, self).__getitem__(item)
        return SlitSpectrogramCube(
            result.data, result.wcs, result.uncertainty, result.unit, result.meta,
            convert_extra_coords_dict_to_input_format(result.extra_coords, result.missing_axes),
            mask=result.mask, missing_axes=result.missing_axes)

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
        result: `SlitSpectrogramCube`
            New SlitSpectrogramCube in new units.

        """
        # Get exposure time in seconds and change array's shape so that
        # it can be broadcast with data and uncertainty arrays.
        exposure_time_s = self.extra_coords["exposure time"]["value"].to(u.s).value
        if not self.extra_coords["exposure time"]["value"].isscalar:
            if len(self.dimensions) == 1:
                pass
            elif len(self.dimensions) == 2:
                exposure_time_s = exposure_time_s[:, np.newaxis]
            elif len(self.dimensions) == 3:
                exposure_time_s = exposure_time_s[:, np.newaxis, np.newaxis]
            else:
                raise ValueError(
                    "SlitSpectrogramCube dimensions must be 2 or 3. Dimensions={0}".format(
                        len(self.dimensions.shape)))
        # Based on value on undo kwarg, apply or remove exposure time correction.
        if undo is True:
            new_data_arrays, new_unit = iris_tools.uncalculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        else:
            new_data_arrays, new_unit = iris_tools.calculate_exposure_time_correction(
                (self.data, self.uncertainty.array), self.unit, exposure_time_s, force=force)
        # Return new instance of SlitSpectrogramCube with correction applied/undone.
        return SlitSpectrogramCube(
            new_data_arrays[0], self.wcs, new_data_arrays[1], new_unit, self.meta,
            convert_extra_coords_dict_to_input_format(self.extra_coords, self.missing_axes),
            mask=self.mask, missing_axes=self.missing_axes)


def read_iris_spectrograph_level2_fits(filenames, spectral_windows=None, uncertainty=True, memmap=False):
    """
    Reads IRIS level 2 spectrograph FITS from an OBS into an IRISSpectrograph instance.

    Parameters
    ----------
    filenames: `list` of `str` or `str`
        Filename of filenames to be read.  They must all be associated with the same
        OBS number.

    spectral_windows: iterable of `str` or `str`
        Spectral windows to extract from files.  Default=None, implies, extract all
        spectral windows.

    Returns
    -------
    result: `irispy.spectrograph.IRISSpectrograph`

    """
    if type(filenames) is str:
        filenames = [filenames]
    for f, filename in enumerate(filenames):
        hdulist = fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap)
        hdulist.verify('fix')
        if f == 0:
            # Determine number of raster positions in a scan
            raster_positions_per_scan = int(hdulist[0].header["NRASTERP"])
            # Collecting the window observations.
            windows_in_obs = np.array([hdulist[0].header["TDESC{0}".format(i)]
                                       for i in range(1, hdulist[0].header["NWIN"]+1)])
            # If spectral_window is not set then get every window.
            # Else take the appropriate windows
            if not spectral_windows:
                spectral_windows_req = windows_in_obs
                window_fits_indices = range(1, len(hdulist)-2)
            else:
                if type(spectral_windows) is str:
                    spectral_windows_req = [spectral_windows]
                else:
                    spectral_windows_req = spectral_windows
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
            top_meta = {"TELESCOP": hdulist[0].header["TELESCOP"],
                        "INSTRUME": hdulist[0].header["INSTRUME"],
                        "DATA_LEV": hdulist[0].header["DATA_LEV"],
                        "OBSID": hdulist[0].header["OBSID"],
                        "OBS_DESC": hdulist[0].header["OBS_DESC"],
                        "STARTOBS": Time(hdulist[0].header["STARTOBS"]),
                        "ENDOBS": Time(hdulist[0].header["ENDOBS"]),
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
            # Create a empty list for every spectral window and each
            # spectral window is a key for the dictionary.
            data_dict = dict([(window_name, list())
                              for window_name in spectral_windows_req])
        # Determine extra coords for this raster.
        times = (Time(hdulist[0].header["STARTOBS"]) +
                 TimeDelta(hdulist[-2].data[:, hdulist[-2].header["TIME"]], format='sec'))
        raster_positions = np.arange(int(hdulist[0].header["NRASTERP"]))
        pztx = hdulist[-2].data[:, hdulist[-2].header["PZTX"]] * u.arcsec
        pzty = hdulist[-2].data[:, hdulist[-2].header["PZTY"]] * u.arcsec
        xcenix = hdulist[-2].data[:, hdulist[-2].header["XCENIX"]] * u.arcsec
        ycenix = hdulist[-2].data[:, hdulist[-2].header["YCENIX"]] * u.arcsec
        obs_vrix = hdulist[-2].data[:, hdulist[-2].header["OBS_VRIX"]] * u.m/u.s
        ophaseix = hdulist[-2].data[:, hdulist[-2].header["OPHASEIX"]]
        exposure_times_fuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEF"]] * u.s
        exposure_times_nuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEN"]] * u.s
        # If OBS is raster, include raster positions.  Otherwise don't.
        if top_meta["NRASTERP"] > 1:
            general_extra_coords = [("time", 0, times),
                                    ("raster position", 0, np.arange(top_meta["NRASTERP"])),
                                    ("pztx", 0, pztx), ("pzty", 0, pzty),
                                    ("xcenix", 0, xcenix), ("ycenix", 0, ycenix),
                                    ("obs_vrix", 0, obs_vrix), ("ophaseix", 0, ophaseix)]
        else:
            general_extra_coords = [("time", 0, times),
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
            # Sit-and-stare have a CDELT of 0 which causes issues in astropy WCS.
            # In this case, set CDELT to a tiny non-zero number.
            if hdulist[window_fits_indices[i]].header["CDELT3"] == 0:
                hdulist[window_fits_indices[i]].header["CDELT3"] = 1e-10
            wcs_ = WCS(hdulist[window_fits_indices[i]].header)
            if not memmap:
                data_mask = hdulist[window_fits_indices[i]].data == -200.
            else:
                data_mask = None
            # Derive extra coords for this spectral window.
            window_extra_coords = copy.deepcopy(general_extra_coords)
            window_extra_coords.append(("exposure time", 0, exposure_times))
            # Collect metadata relevant to single files.
            try:
                date_obs = Time(hdulist[0].header["DATE_OBS"])
            except ValueError:
                date_obs = None
            try:
                date_end = Time(hdulist[0].header["DATE_END"])
            except ValueError:
                date_end = None
            single_file_meta = {"SAT_ROT": hdulist[0].header["SAT_ROT"] * u.deg,
                                "DATE_OBS": date_obs,
                                "DATE_END": date_end,
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
                                "OBSID": hdulist[0].header["OBSID"],
                                "OBS_DESC": hdulist[0].header["OBS_DESC"],
                                "STARTOBS": Time(hdulist[0].header["STARTOBS"]),
                                "ENDOBS": Time(hdulist[0].header["ENDOBS"])
                                }
            # Derive uncertainty of data
            if uncertainty:
                out_uncertainty = u.Quantity(np.sqrt(
                    (hdulist[window_fits_indices[i]].data*DN_unit).to(u.photon).value +
                    readout_noise.to(u.photon).value**2), unit=u.photon).to(DN_unit).value
            else:
                out_uncertainty = None
            # Appending NDCube instance to the corresponding window key in dictionary's list.
            data_dict[window_name].append(
                SlitSpectrogramCube(hdulist[window_fits_indices[i]].data, wcs_, out_uncertainty,
                                    DN_unit, single_file_meta, window_extra_coords,
                                    mask=data_mask))
        hdulist.close()
    # Construct dictionary of SlitSpectrogramCubeSequences for spectral windows
    data = dict([(window_name, SlitSpectrogramCubeSequence(data_dict[window_name],
                                                           window_metas[window_name],
                                                           common_axis=0))
                 for window_name in spectral_windows_req])
    # Initialize an IRISSpectrograph object.
    return IRISSpectrograph(data, meta=top_meta)
