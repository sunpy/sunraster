# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

from datetime import timedelta
from collections import OrderedDict
import copy

import numpy as np
import xarray
from astropy.io import fits
import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.table import Table, vstack
from astropy import wcs
from astropy import constants
from scipy import interpolate
from sunpy.time import parse_time

import iris_tools


class IRISRaster(object):
    """An object to hold data from multiple IRIS raster scans."""
    def __init__(self, filenames, spectral_windows="All"):
        """Initializes an IRISRaster object from IRIS level 2 files."""
        # If a single filename has been entered as a string, convert
        # to a list of length 1 for consistent syntax below.
        if type(filenames) is str:
            filenames = [filenames]
        # Define some empty variables.
        wcs_celestial_objects = dict()
        wcs_spectral_objects = dict()
        raster_index_to_file = []
        raster_positions = []
        # Open files and extract data.
        for f, filename in enumerate(filenames):
            hdulist = fits.open(filename)
            # If this is the first file, extract some common metadata.
            if f == 0:
                # Check user desired spectral windows are in file and
                # find corresponding indices of HDUs.
                n_win = int(hdulist[0].header["NWIN"])
                windows_in_obs = np.array([hdulist[0].header["TDESC{0}".format(i)] for i in range(1, n_win+1)])
                if spectral_windows == "All":
                    spectral_windows = windows_in_obs
                    window_fits_indices = range(1, len(hdulist)-2)
                else:
                    if type(spectral_windows) is str:
                        spectral_windows = [spectral_windows]
                    spectral_windows = np.asarray(spectral_windows, dtype="U")
                    window_is_in_obs = np.asarray([window in windows_in_obs for window in spectral_windows])
                    if not all(window_is_in_obs):
                        missing_windows = window_is_in_obs == False
                        raise ValueError(
                            "Spectral windows {0} not in file {1}".format(spectral_windows[missing_windows],
                                                                          filenames[0]))
                    window_fits_indices = np.nonzero(np.in1d(windows_in_obs, spectral_windows))[0]+1
                # Create table of spectral window info in OBS.
                self.spectral_windows = Table([
                    [hdulist[0].header["TDESC{0}".format(i)] for i in window_fits_indices],
                    [hdulist[0].header["TDET{0}".format(i)] for i in window_fits_indices],
                    Quantity([hdulist[0].header["TWAVE{0}".format(i)] for i in window_fits_indices], unit="angstrom"),
                    Quantity([hdulist[0].header["TWMIN{0}".format(i)] for i in window_fits_indices], unit="angstrom"),
                    Quantity([hdulist[0].header["TWMAX{0}".format(i)] for i in window_fits_indices], unit="angstrom")],
                    names=("name", "detector type", "brightest wavelength", "min wavelength", "max wavelength"))
                # Set spectral window name as table index
                self.spectral_windows.add_index("name")
                # Find wavelength represented by each pixel in the
                # spectral dimension by using a WCS object for each spectral
                # window.
                spectral_coords = dict()
                for i, window_name in enumerate(self.spectral_windows["name"]):
                    wcs_spectral = wcs.WCS(hdulist[window_fits_indices[i]].header).sub(1)
                    spectral_coords[window_name] = Quantity(wcs_spectral.all_pix2world(np.arange(
                        hdulist[window_fits_indices[i]].header["NAXIS1"]), 0),
                        unit=wcs_spectral.wcs.cunit[0]).to("Angstrom")[0]
                    wcs_spectral_objects[window_name] = wcs_spectral
                # Put useful metadata into meta attribute.
                self.meta = {"date data created": parse_time(hdulist[0].header["DATE"]),
                             "telescope": hdulist[0].header["TELESCOP"],
                             "instrument": hdulist[0].header["INSTRUME"],
                             "data level": hdulist[0].header["DATA_LEV"],
                             "level 2 reformatting version": hdulist[0].header["VER_RF2"],
                             "level 2 reformatting date": parse_time(hdulist[0].header["DATE_RF2"]),
                             "DATA_SRC": hdulist[0].header["DATA_SRC"],
                             "origin": hdulist[0].header["origin"],
                             "build version": hdulist[0].header["BLD_VERS"],
                             "look-up table ID": hdulist[0].header["LUTID"],
                             "observation ID": int(hdulist[0].header["OBSID"]),
                             "observation description": hdulist[0].header["OBS_DESC"],
                             "observation label": hdulist[0].header["OBSLABEL"],
                             "observation title": hdulist[0].header["OBSTITLE"],
                             "observation start": hdulist[0].header["STARTOBS"],
                             "observation end": hdulist[0].header["ENDOBS"],
                             "observation repetitions": hdulist[0].header["OBSREP"],
                             "camera": hdulist[0].header["CAMERA"],
                             "status": hdulist[0].header["STATUS"],
                             "data quantity": hdulist[0].header["BTYPE"],
                             "data unit": hdulist[0].header["BUNIT"],
                             "BSCALE": hdulist[0].header["BSCALE"],
                             "BZERO": hdulist[0].header["BZERO"],
                             "high latitude flag": hdulist[0].header["HLZ"],
                             "SAA": bool(int(hdulist[0].header["SAA"])),
                             "satellite roll angle": Quantity(float(hdulist[0].header["SAT_ROT"]), unit=u.deg),
                             "AEC exposures in OBS": hdulist[0].header["AECNOBS"],
                             "dsun": Quantity(hdulist[0].header["DSUN_OBS"], unit="m"),
                             "IAECEVFL": bool(),
                             "IAECFLAG": bool(),
                             "IAECFLFL": bool(),
                             "FOV Y axis": Quantity(float(hdulist[0].header["FOVY"]), unit="arcsec"),
                             "FOV X axis": Quantity(float(hdulist[0].header["FOVX"]), unit="arcsec"),
                             "FOV center Y axis": Quantity(float(hdulist[0].header["YCEN"]), unit="arcsec"),
                             "FOV center X axis": Quantity(float(hdulist[0].header["XCEN"]), unit="arcsec"),
                             "spectral summing NUV": hdulist[0].header["SUMSPTRN"],
                             "spectral summing FUV": hdulist[0].header["SUMSPTRF"],
                             "spatial summing": hdulist[0].header["SUMSPAT"],
                             "exposure time mean": hdulist[0].header["EXPTIME"],
                             "exposure time min": hdulist[0].header["EXPMIN"],
                             "exposure time max": hdulist[0].header["EXPMAX"],
                             "total exposures in OBS": hdulist[0].header["NEXPOBS"],
                             "number unique raster positions": hdulist[0].header["NRASTERP"],
                             "raster step size mean": hdulist[0].header["STEPS_AV"],
                             "raster step size sigma": hdulist[0].header["STEPS_DV"],
                             "time step size mean": hdulist[0].header["STEPT_AV"],
                             "time step size sigma": hdulist[0].header["STEPT_DV"],
                             "spectral windows in OBS": windows_in_obs,
                             "spectral windows in object": spectral_windows,
                             "detector gain": iris_tools.DETECTOR_GAIN,
                             "detector yield": iris_tools.DETECTOR_YIELD,
                             "readout noise": iris_tools.READOUT_NOISE}
                # Translate some metadata to be more helpful.
                if hdulist[0].header["IAECEVFL"] == "YES":
                    self.meta["IAECEVFL"] = True
                if hdulist[0].header["IAECFLAG"] == "YES":
                    self.meta["IAECFLAG"] = True
                if hdulist[0].header["IAECFLFL"] == "YES":
                    self.meta["IAECFLFL"] = True
                if self.meta["data level"] == 2.:
                    if self.meta["camera"] == 1:
                        self.meta["camera"] = "spectra"
                    elif self.meta["camera"] == 2:
                        self.meta["camera"] = "SJI"
                # Define empty dictionary with keys corresponding to
                # spectral windows.  The value of each key will be a
                # list of xarray data arrays, one for each raster scan.
                data_dict = dict([(window_name, None) for window_name in self.spectral_windows["name"]])
                # Record header info of auxiliary data.  Should be
                # consistent between files of same OBS.
                auxiliary_header = hdulist[-2].header
            # Extract the data and meta/auxiliary data.
            # Create WCS object from FITS header and add WCS object
            # wcs dictionary.
            wcs_celestial = wcs.WCS(hdulist[1].header).celestial
            scan_label = "scan{0}".format(f)
            wcs_celestial_objects[scan_label] = wcs_celestial
            # Append to list representing the scan labels of each
            # spectrum.
            len_raster_axis = hdulist[1].header["NAXIS3"]
            raster_index_to_file = raster_index_to_file+[scan_label]*len_raster_axis
            # Append to list representing the raster positions of each
            # spectrum.
            raster_positions = raster_positions+list(range(len_raster_axis))
            # Concatenate auxiliary data arrays from each file.
            try:
                auxiliary_data = np.concatenate((auxiliary_data, np.array(hdulist[-2].data)), axis=0)
            except UnboundLocalError as e:
                if e.args[0] == "local variable 'auxiliary_data' referenced before assignment":
                    auxiliary_data = np.array(hdulist[-2].data)
                else:
                    raise e
            # For each spectral window, concatenate data from each file.
            for i, window_name in enumerate(self.spectral_windows["name"]):
                # Set invalid data values to NaN.
                data_nan_masked = copy.deepcopy(hdulist[window_fits_indices[i]].data)
                data_nan_masked[hdulist[window_fits_indices[i]].data == -200.] = np.nan
                try:
                    data_dict[window_name] = np.concatenate((data_dict[window_name], data_nan_masked))
                except ValueError as e:
                    if e.args[0] == "zero-dimensional arrays cannot be concatenated":
                        data_dict[window_name] = data_nan_masked
                    else:
                        raise e
            # Close file.
            hdulist.close()
        # Having combined various data from files into common objects,
        # convert into final data formats and attach to class.
        # Convert auxiliary data into Table and attach to class.
        self.auxiliary_data = Table()
        # Enter certain properties into auxiliary data table as
        # quantities with units.
        auxiliary_colnames = [key for key in auxiliary_header.keys()][7:]
        quantity_colnames = [("TIME", "s"), ("PZTX", "arcsec"), ("PZTY", "arcsec"),
                             ("EXPTIMEF", "s"), ("EXPTIMEN", "s"), ("XCENIX", "arcsec"),
                             ("YCENIX", "arcsec"), ("OBS_VRIX", "m/s")]
        for col in quantity_colnames:
            self.auxiliary_data[col[0]] = _enter_column_into_table_as_quantity(
                col[0], auxiliary_header, auxiliary_colnames, auxiliary_data, col[1])
        # Enter remaining properties into table without units/
        for i, colname in enumerate(auxiliary_colnames):
            self.auxiliary_data[colname] = auxiliary_data[:, auxiliary_header[colname]]
        # Reorder columns so they reflect order in data file.
        self.auxiliary_data = self.auxiliary_data[[key for key in auxiliary_header.keys()][7:]]
        # Rename some columns to be more user friendly.
        rename_colnames = [("EXPTIMEF", "FUV EXPOSURE TIME"), ("EXPTIMEN", "NUV EXPOSURE TIME")]
        for col in rename_colnames:
            self.auxiliary_data.rename_column(col[0], col[1])
        # Add column designating what scan/file number each spectra
        # comes from.  This can be used to determine the corresponding
        # wcs object and level 1 info.
        self.auxiliary_data["scan"] = raster_index_to_file
        # Attach dictionary containing level 1 and wcs info for each file used.
        self.wcs = {"spectral": wcs_spectral_objects, "celestial": wcs_celestial_objects}
        # Calculate measurement time of each spectrum.
        times = [parse_time(self.meta["observation start"])+timedelta(seconds=s) for s in self.auxiliary_data["TIME"]]
        # Convert data for each spectral window into an an
        # xarray.DataArray and enter into data dictionary.
        self.data = dict([(window_name, xarray.DataArray(data=data_dict[window_name],
                                                         dims=["raster_axis", "slit_axis", "spectral_axis"],
                                                         coords={"wavelength": ("spectral_axis",
                                                                                spectral_coords[window_name].value),
                                                                 "raster_position": ("raster_axis", raster_positions),
                                                                 "time": ("raster_axis", times)},
                                                         name="Intensity [DN]",
                                                         attrs=OrderedDict([(
                                                             "units", {"wavelength": spectral_coords[window_name].unit,
                                                                       "intensity": "DN"})])))
                          for window_name in self.spectral_windows["name"]])


    def __repr__(self):
        spectral_window = self.spectral_windows["name"][0]
        spectral_windows_info = "".join(
            ["\n    {0}\n        (raster axis: {1}, slit axis: {2}, spectral axis: {3})".format(
            name, len(self.data[name].raster_axis), len(self.data[name].slit_axis),
            len(self.data[name].spectral_axis)) for name in self.spectral_windows["name"]])
        return "<iris.IRISRaster instance\nOBS ID: {0}\n".format(self.meta["observation ID"]) + \
               "OBS Description: {0}\n".format(self.meta["observation description"]) + \
               "OBS period: {0} -- {1}\n".format(self.meta["observation start"], self.meta["observation end"]) + \
               "Instance period: {0} -- {1}\n".format(self.data[spectral_window].time.values[0],
                                                    self.data[spectral_window].time.values[-1]) + \
               "Number unique raster positions: {0}\n".format(self.meta["number unique raster positions"]) + \
               "Spectral windows{0}>".format(spectral_windows_info)


    def convert_DN_to_photons(self, spectral_window):
        """Converts DataArray from DN to photon counts."""
        # Check that DataArray is in units of DN.
        if "DN" not in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError("Intensity units of DataArray are not DN.")
        self.data[spectral_window].data = iris_tools.convert_DN_to_photons(spectral_window)
        self.data[spectral_window].name = "Intensity [photons]"
        self.data[spectral_window].atrrs["units"]["intensity"] = "photons"


    def convert_photons_to_DN(self, spectral_window):
        """Converts DataArray from DN to photon counts."""
        # Check that DataArray is in units of DN.
        if "photons" not in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError("Intensity units of DataArray are not DN.")
        self.data[spectral_window].data = iris_tools.convert_photons_to_DN(spectral_window)
        self.data[spectral_window].name = "Intensity [DN]"
        self.data[spectral_window].atrrs["units"]["intensity"] = "DN"


    def apply_exposure_time_correction(self, spectral_window):
        """Converts DataArray from DN or photons to DN or photons per second."""
        # Check that DataArray is in units of DN.
        if "/s" in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError("Data seems to already be in units per second. '/s' in intensity unit string.")
        detector_type = self.spectral_windows.loc[spectral_window]["detector type"][:3]
        exp_time_s = self.auxiliary_data["{0} EXPOSURE TIME".format(detector_type)].to("s").value
        for i in self.data[spectral_window].raster_axis.values:
            self.data[spectral_window].data[i, :, :] = self.data[spectral_window].data[i, :, :]/exp_time_s[i]
        # Make new unit reflecting the division by time.
        unit_str = self.data[spectral_window].attrs["units"]["intensity"]+"/s"
        self.data[spectral_window].attrs["units"]["intensity"] = unit_str
        name_split = self.data[spectral_window].name.split("[")
        self.data[spectral_window].name = "{0}[{1}]".format(name_split[0], unit_str)


    def calculate_intensity_fractional_uncertainty(self, spectral_window):
        return iris_tools.calculate_intensity_fractional_uncertainty(
            self.data[spectral_window].data, self.data[spectral_window].attrs["units"]["intensity"],
            self.spectral_windows.loc[spectral_window]["detector type"][:3])


    def convert_DN_to_radiance(self, spectral_window):
        """Convert DataArray from count rate [DN/s/pixel] to radiance [ergs/s/cm^2/sr/]."""
        # Get spectral dispersion per pixel.
        ######## Must check if data is in units of DN/s ############
        spectral_dispersion_per_pixel = self.wcs["spectral"][spectral_window].cdelt[0] * u.Angstrom
        ########### The definition of solid angle is very rough as defined. #############
        ########### Needs to be made into an array for each scan and also the 0.33"
        ########### width of the slit must be found in the data file somewhere. ############
        solid_angle = self.wcs["celestial"]["scan0"].cdelt[0] * u.steradian
        # Get effective area
        ########### This needs to be generalized to the time of OBS once that functionality is written#########
        iris_response = iris_tools.iris_get_response(pre_launch=True)
        lam = iris_response["LAMBDA"]
        if self.spectral_windows[spectral_window]["detector type"][:3] == "FUV":
            eff_area = iris_response.area_sg[0, :]
            dn2phot = iris_response["DN2PHOT_SG"][0]
        elif self.spectral_windows[spectral_window]["detector type"][:3] == "NUV":
            eff_area = iris_response.area_sg[1, :]
            dn2phot = iris_response["DN2PHOT_SG"][1]
        else:
            raise ValueError("Detector type of spectral window not recognized.")
        # Iterate through raster and slit pixel and make conversion to
        # physical units.
        for i in range(len(self.data[spectral_window].raster_axis)):
            # Interpolate the effective areas to cover the wavelengths
            # at which the data is recorded:
            tck = interpolate.splrep(lam, eff_area.to(u.Angstrom).value, s=0)
            eff_area_interp = interpolate.splev(wave[i, :], tck)
            # Iterate over pixels in slit.
            for k in range(len(self.data[spectral_window].slit_axis)):
                data_rad[i, k, :] = constants.h * constants.c / wave[i, :] * u.Angstrom * \
                                    self.data[spectral_window].data.isel(raster_axis=i).isel(slit_axis=k) * \
                                    dn2phot / solid_angle / (eff_area_interp * spectral_dispersion_per_pixel)


    def calculate_orbital_wavelength_variation(self, slit_pixel_range=None, spline_smoothing=False,
                                               fit_individual_profiles=False):
        if date_created < date_new_pipeline:
            spacecraft_velocity = raster.auxiliary_data["OBS_VRIX"]
            orbital_phase = 2. * np.pi * raster.auxiliary_data["OPHASEIX"]
            roll_angle = raster.meta["satellite roll angle"]
            # Check that there are measurement times with good values of
            # spacecraft velocity and orbital phase.
            bad_aux = np.asarray(np.isfinite(spacecraft_velocity) * np.isfinite(orbital_phase) * (-1), dtype=bool)
        else:
            spacecraft_velocity = None
            orbital_phase = None
            roll_angle = None


def _enter_column_into_table_as_quantity(header_property_name, header, header_colnames, data, unit):
    """Used in initiation of IRISRaster to convert auxiliary data to Quantities."""
    index = np.where(np.array(header_colnames) == header_property_name)[0]
    if len(index) == 1:
        index = index[0]
    else:
        raise ValueError("Multiple property names equal to {0}".format(header_property_name))
    pop_colname = header_colnames.pop(index)
    return Quantity(data[:, header[pop_colname]], unit=unit)