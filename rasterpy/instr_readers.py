import numpy as np
from astropy.io import fits
from astropy.time import Time

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
                IRISSpectrogramCube(hdulist[window_fits_indices[i]].data, wcs_, out_uncertainty,
                                    DN_unit, single_file_meta, window_extra_coords,
                                    mask=data_mask))
        hdulist.close()
    # Construct dictionary of IRISSpectrogramCubeSequences for spectral windows
    data = dict([(window_name, IRISSpectrogramCubeSequence(data_dict[window_name],
                                                           window_metas[window_name],
                                                           common_axis=0))
                 for window_name in spectral_windows_req])
    # Initialize an IRISSpectrograph object.
    return IRISSpectrograph(data, meta=top_meta)
