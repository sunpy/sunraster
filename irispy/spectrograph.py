# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

from datetime import timedelta
from collections import OrderedDict
from astropy.units.quantity import Quantity
from astropy.table import Table
from astropy import constants
from scipy import interpolate
from sunpy.time import parse_time
from astropy.io import fits
from sunpycube.cube.datacube import Cube, CubeSequence
from sunpycube.wcs_util import WCS

import copy
import numpy as np
import astropy.units as u
import irispy.iris_tools as iris_tools


__all__ = ['IRISSG']


class IRISSG(object):
    """An object to hold data from multiple IRIS raster scans."""

    def __init__(self, filenames, spectral_windows="All", common_axis=0):
        """Initializes an IRISSG object from IRIS level 2 files."""
        # default common axis is 0.
        if type(filenames) is str:
            filenames = [filenames]

        for f, filename in enumerate(filenames):
            hdulist = fits.open(filename)
            if f == 0:
                # collecting the window observations.
                windows_in_obs = np.array([hdulist[0].header["TDESC{0}".format(i)]
                                           for i in range(1, hdulist[0].header["NWIN"]+1)])
                # if spectral window is All then get every window else take the appropriate windows
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
                        raise ValueError(
                            "Spectral windows {0} not in file {1}".format(spectral_windows[missing_windows],
                                                                          filenames[0]))
                    window_fits_indices = np.nonzero(np.in1d(windows_in_obs, spectral_windows))[0]+1
                # table to store hdulist headers of 0th index.
                self.spectral_windows = Table([
                    [hdulist[0].header["TDESC{0}".format(i)] for i in window_fits_indices],
                    [hdulist[0].header["TDET{0}".format(i)] for i in window_fits_indices],
                    Quantity([hdulist[0].header["TWAVE{0}".format(i)]
                              for i in window_fits_indices], unit="angstrom"),
                    Quantity([hdulist[0].header["TWMIN{0}".format(i)]
                              for i in window_fits_indices], unit="angstrom"),
                    Quantity([hdulist[0].header["TWMAX{0}".format(i)] for i in window_fits_indices], unit="angstrom")],
                    names=("name", "detector type", "brightest wavelength", "min wavelength", "max wavelength"))
                # all spectral "required/all" are given name index
                self.spectral_windows.add_index("name")
                # creating a empty list for every spectral window and each spectral window
                # is a key for the dictionary.
                data_dict = dict([(window_name, list())
                                  for window_name in self.spectral_windows["name"]])

            # the unchanged header of the hdulist indexed 0.
            self.meta = hdulist[0].header
            for i, window_name in enumerate(self.spectral_windows["name"]):
                wcs_ = WCS(hdulist[window_fits_indices[i]].header)
                data_ = hdulist[window_fits_indices[i]].data
                try:
                    # appending Cube instance to the corresponding window key in dictionary's list.
                    data_dict[window_name].append(Cube(data_, wcs_, meta=dict(self.meta)))
                except ValueError as e:
                    raise e
            hdulist.close()
        # making a CubeSequence of every dictionary key window.
        self.data = dict([(window_name, CubeSequence(data_dict[window_name], meta=self.meta, common_axis=common_axis))
                          for window_name in self.spectral_windows['name']])

    def __repr__(self):
        spectral_window = self.spectral_windows["name"][0]
        spectral_windows_info = "".join(
            ["\n    {0}\n        (raster axis: {1}, slit axis: {2}, spectral axis: {3})".format(
                name,
                self.data[name][0].shape[0],
                self.data[name][0].shape[1],
                self.data[name][0].shape[2])
                for name in self.spectral_windows["name"]])
        # Removed the instance Period. if required let me know
        return "<iris.IRISSG instance\nOBS ID: {0}\n".format(self.meta["OBSID"]) + \
               "OBS Description: {0}\n".format(self.meta["OBS_DESC"]) + \
               "OBS period: {0} -- {1}\n".format(self.meta["STARTOBS"], self.meta["ENDOBS"]) + \
               "Number unique raster positions: {0}\n".format(self.meta["NRASTERP"]) + \
               "Spectral windows{0}>".format(spectral_windows_info)

    def convert_DN_to_photons(self, spectral_window):
        """Converts DataArray from DN to photon counts."""
        # Check that DataArray is in units of DN.
        if "DN" not in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError("Intensity units of DataArray are not DN.")
        self.data[spectral_window].data = iris_tools.convert_DN_to_photons(
            self.data[spectral_window], self.spectral_windows.loc[spectral_window]['detector type'])
        self.data[spectral_window].name = "Intensity [photons]"
        self.data[spectral_window].attrs["units"]["intensity"] = "photons"

    def convert_photons_to_DN(self, spectral_window):
        """Converts DataArray from DN to photon counts."""
        # Check that DataArray is in units of DN.
        if "photons" not in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError("Intensity units of DataArray are not DN.")
        self.data[spectral_window].data = iris_tools.convert_photons_to_DN(
            self.data[spectral_window], self.spectral_windows.loc[spectral_window]['detector type'])
        self.data[spectral_window].name = "Intensity [DN]"
        self.data[spectral_window].attrs["units"]["intensity"] = "DN"

    def apply_exposure_time_correction(self, spectral_window):
        """Converts DataArray from DN or photons to DN or photons per second."""
        # Check that DataArray is in units of DN.
        if "/s" in self.data[spectral_window].attrs["units"]["intensity"]:
            raise ValueError(
                "Data seems to already be in units per second. '/s' in intensity unit string.")
        detector_type = self.spectral_windows.loc[spectral_window]["detector type"][:3]
        exp_time_s = self.auxiliary_data["{0} EXPOSURE TIME".format(detector_type)].to("s").value
        for i in self.data[spectral_window].raster_axis.values:
            self.data[spectral_window].data[i, :, :] = self.data[
                spectral_window].data[i, :, :]/exp_time_s[i]
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
        # Needs to be made into an array for each scan and also the 0.33"
        ########### width of the slit must be found in the data file somewhere. ############
        solid_angle = self.wcs["celestial"]["scan0"].cdelt[0] * u.steradian
        # Get effective area
        # This needs to be generalized to the time of OBS once that func
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
            bad_aux = np.asarray(np.isfinite(spacecraft_velocity) *
                                 np.isfinite(orbital_phase) * (-1), dtype=bool)
        else:
            spacecraft_velocity = None
            orbital_phase = None
            roll_angle = None


def _enter_column_into_table_as_quantity(header_property_name, header, header_colnames, data, unit):
    """Used in initiation of IRISSG to convert auxiliary data to Quantities."""
    index = np.where(np.array(header_colnames) == header_property_name)[0]
    if len(index) == 1:
        index = index[0]
    else:
        raise ValueError("Multiple property names equal to {0}".format(header_property_name))
    pop_colname = header_colnames.pop(index)
    return Quantity(data[:, header[pop_colname]], unit=unit)
