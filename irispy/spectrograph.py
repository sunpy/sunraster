# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>


import copy
from datetime import timedelta
from collections import namedtuple

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from ndcube import NDCube, NDCubeSequence
from ndcube.wcs_util import WCS
from sunpy.time import parse_time

from irispy.spectrogram import SpectrogramSequence

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
        raster_index_to_file = []
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
                    window_fits_indices = np.nonzero(np.in1d(
                        windows_in_obs, spectral_windows))[0]+1
                # Create table of spectral window info in OBS.
                self.spectral_windows = Table([
                    [hdulist[0].header["TDESC{0}".format(i)] for i in window_fits_indices],
                    [hdulist[0].header["TDET{0}".format(i)] for i in window_fits_indices],
                    u.Quantity([hdulist[0].header["TWAVE{0}".format(i)]
                                for i in window_fits_indices], unit="angstrom"),
                    u.Quantity([hdulist[0].header["TWMIN{0}".format(i)]
                                for i in window_fits_indices], unit="angstrom"),
                    u.Quantity([hdulist[0].header["TWMAX{0}".format(i)]
                                for i in window_fits_indices], unit="angstrom")],
                    names=("name", "detector type", "brightest wavelength",
                           "min wavelength", "max wavelength"))
                # Set spectral window name as table index.
                self.spectral_windows.add_index("name")
                # creating a empty list for every spectral window and each spectral window
                # is a key for the dictionary.
                data_dict = dict([(window_name, list())
                                  for window_name in self.spectral_windows["name"]])
                auxiliary_header = hdulist[-2].header
            # the unchanged header of the hdulist indexed 0.
            self.meta = hdulist[0].header
            # Calculate measurement time of each spectrum.
            times = np.array([parse_time(self.meta["STARTOBS"]) + timedelta(seconds=s)
                              for s in hdulist[-2].data[:,hdulist[-2].header["TIME"]]])
            raster_positions = np.arange(int(hdulist[0].header["NRASTERP"]))
            extra_coords = [("time", 0, times), ("raster position", 0, raster_positions)]
            for i, window_name in enumerate(self.spectral_windows["name"]):
                wcs_ = WCS(hdulist[window_fits_indices[i]].header)
                data_nan_masked = copy.deepcopy(hdulist[window_fits_indices[i]].data)
                data_nan_masked[hdulist[window_fits_indices[i]].data == -200.] = np.nan
                data_mask = hdulist[window_fits_indices[i]].data == -200.
                # appending NDCube instance to the corresponding
                # window key in dictionary's list.
                data_dict[window_name].append(
                    NDCube(data_nan_masked, wcs=wcs_, meta=dict(self.meta),
                           mask=data_mask, extra_coords=extra_coords))

            scan_label = "scan{0}".format(f)
            # Append to list representing the scan labels of each
            # spectrum.
            len_raster_axis = hdulist[1].header["NAXIS3"]
            raster_index_to_file = raster_index_to_file+[scan_label]*len_raster_axis
            # Concatenate auxiliary data arrays from each file.
            try:
                auxiliary_data = np.concatenate(
                    (auxiliary_data, np.array(hdulist[-2].data)), axis=0)
            except UnboundLocalError as e:
                if e.args[0] == "local variable 'auxiliary_data' referenced before assignment":
                    auxiliary_data = np.array(hdulist[-2].data)
                else:
                    raise e
            hdulist.close()

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
        # making a NDCubeSequence of every dictionary key window.
        self.data = dict([(window_name,
                           SpectrogramSequence(data_dict[window_name], common_axis,
                                               raster_positions_per_scan,
                                               first_exposure_raster_position=0,
                                               meta=self.meta))
                          for window_name in self.spectral_windows['name']])

    def __repr__(self):
        spectral_window = self.spectral_windows["name"][0]
        spectral_windows_info = "".join(
            ["\n    {0}\n        (raster axis, slit axis, spectral axis) {1}".format(
                name,
                self.data[name].dimensions[1::])
                for name in self.spectral_windows["name"]])
        return "<iris.IRISSpectrograph instance\nOBS ID: {0}\n".format(self.meta["OBSID"]) + \
               "OBS Description: {0}\n".format(self.meta["OBS_DESC"]) + \
               "OBS period: {0} -- {1}\n".format(self.meta["STARTOBS"], self.meta["ENDOBS"]) + \
               "Instance period: {0} -- {1}\n".format(
                   self.data[spectral_window][0]._extra_coords["time"]["value"],
                   self.data[spectral_window][-1]._extra_coords["time"]["value"]) + \
               "Number unique raster positions: {0}\n".format(self.meta["NRASTERP"]) + \
               "Spectral windows{0}>".format(spectral_windows_info)

    # A tuple giving coordinate names of axes in NDCubeSequences
    coord_names = ("raster number", "x", "y", "wavelength")
    coord_names_index_as_cube = ("exposure number", "y", "wavelength")


def _enter_column_into_table_as_quantity(header_property_name, header, header_colnames,
                                         data, unit):
    """Used in initiation of IRISSpectrograph to convert auxiliary data to Quantities."""
    index = np.where(np.array(header_colnames) == header_property_name)[0]
    if len(index) == 1:
        index = index[0]
    else:
        raise ValueError("Multiple property names equal to {0}".format(header_property_name))
    pop_colname = header_colnames.pop(index)
    return u.Quantity(data[:, header[pop_colname]], unit=unit)