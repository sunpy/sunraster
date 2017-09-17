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
            times = np.array([parse_time(hdulist[0].header["STARTOBS"]) + timedelta(seconds=s)
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
                # Derive WCS, data and mask for NDCube from file.
                wcs_ = WCS(hdulist[window_fits_indices[i]].header)
                data_nan_masked = copy.deepcopy(hdulist[window_fits_indices[i]].data)
                data_nan_masked[hdulist[window_fits_indices[i]].data == -200.] = np.nan
                data_mask = hdulist[window_fits_indices[i]].data == -200.
                # Derive extra coords for this spectral window.
                window_extra_coords = copy.deepcopy(general_extra_coords)
                if "FUV" in hdulist[0].header["TDET{0}".format(window_fits_indices[i])]:
                    exposure_times = exposure_times_fuv
                else:
                    exposure_times = exposure_times_nuv
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
                        "HISTORY": hdulist[0].header["HISTORY"]
                        }
                # Appending NDCube instance to the corresponding window key in dictionary's list.
                data_dict[window_name].append(
                    NDCube(data_nan_masked, wcs=wcs_, meta=meta, mask=data_mask,
                           extra_coords=window_extra_coords))
            hdulist.close()
        # Attach dictionary containing level 1 and wcs info for each file used.
        # making a NDCubeSequence of every dictionary key window.
        self.data = dict([(window_name,
                           SpectrogramSequence(data_dict[window_name], common_axis,
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
            if type(self.data[key]) == SpectrogramSequence:
                spectral_window_list.append([self.data[key].meta[colname] for colname in colnames])
        return Table(rows=spectral_window_list, names=colnames)
