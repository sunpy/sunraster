# -*- coding: utf-8 -*-
# Author: Daniel Ryan <ryand5@tcd.ie>

"""Some IRIS instrument tools."""

import datetime
import warnings
import os.path

import numpy as np
import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.modeling import fitting
from astropy.modeling.models import custom_model
from astropy import constants
import scipy.io
from scipy import interpolate
from sunpy.time import parse_time
import sunpy.util.config
from sunpy.util.net import check_download_file
from ndcube import NDCube

# Define some properties of IRIS detectors.  Source: IRIS instrument
# paper.
DETECTOR_GAIN = {"NUV": 18., "FUV": 6., "SJI": 18.}
DETECTOR_YIELD = {"NUV": 1., "FUV": 1.5, "SJI": 1.}
DN_UNIT = {
    "NUV": u.def_unit("DN_IRIS_NUV", DETECTOR_GAIN["NUV"]/DETECTOR_YIELD["NUV"]*u.photon),
    "FUV": u.def_unit("DN_IRIS_FUV", DETECTOR_GAIN["FUV"]/DETECTOR_YIELD["FUV"]*u.photon),
    "SJI": u.def_unit("DN_IRIS_SJI", DETECTOR_GAIN["SJI"]/DETECTOR_YIELD["SJI"]*u.photon)}
READOUT_NOISE = {"NUV": 1.2*DN_UNIT["NUV"], "FUV": 3.1*DN_UNIT["FUV"],
                 "SJI": 1.2*DN_UNIT["SJI"]}
RADIANCE_UNIT = u.erg / u.cm ** 2 / u.s / u.steradian / u.Angstrom
SLIT_WIDTH = 0.33*u.arcsec

IRIS_RESPONSE_REMOTE_PATH = "https://sohowww.nascom.nasa.gov/solarsoft/iris/response/"
RESPONSE_VERSION_FILENAMES = {"1": "iris_sra_20130211.geny", "2": "iris_sra_20130715.geny",
                              "3": "iris_sra_c_20150331.geny", "4": "iris_sra_c_20161022.geny"}

# Define some custom error messages.
APPLY_EXPOSURE_TIME_ERROR = "Exposure time correction has probably already " + \
    "been applied since the unit already includes inverse time.  " + \
    "To apply exposure time correction anyway, set 'force' kwarg to True."
UNDONE_EXPOSURE_TIME_ERROR = "Exposure time correction has probably already " + \
    "been undone since the unit does not include inverse time.  " + \
    "To undo exposure time correction anyway, set 'force' kwarg to True."

# Define whether IRIS WCS is 0 or 1 origin based.
WCS_ORIGIN = 1

def get_iris_response(pre_launch=False, response_file=None, response_version=None,
                      force_download=False):
    """Returns IRIS response structure.

    One and only one of pre_launch, response_file and response_version must be set.

    Parameters
    ----------
    pre_launch: `bool`
        Equivalent to setting response_version=2.  Cannot be set
        simultaneously with response_file kwarg. Default=False
    response_file: `str`
        Version number of effective area file to be used.  Cannot be set
        simultaneously with pre_launch kwarg.  Default=latest
    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=latest

    Returns
    -------
    iris_response: `dict`
        Various parameters regarding IRIS response.  The following keys:
        date_obs: `datetime.datetime`
        lambda: `astropy.units.Quantity`
        area_sg: `astropy.units.Quantity`
        name_sg: `str`
        dn2phot_sg: `tuple` of length 2
        area_sji: `astropy.units.Quantity`
        name_sji: `str`
        dn2phot_sji:  `tuple` of length 4
        comment: `str`
        version: `int`
        version_date: `datetime.datetime`

    Notes
    -----
    This routine does not calculate time dependent effective areas using
    version 3 and above of the response functions as is done in the SSW version
    of this code.  Therefore, asking it to read a version 3 or above response
    function will result in an error.  This code should be updated in future
    versions to calculate time dependent effective areas.

    """
    # Ensures the file exits in the path given.
    if response_file is not None:
        if not(os.path.isfile(response_file)):
            raise KeyError("Not a valid file path")

    
    # Ensure conflicting kwargs are not set.
    if response_file:
        response_file_set = True
    else:
        response_file_set = False
    if response_version:
        response_version_set = True
    else:
        response_version_set = False
    if response_file_set+pre_launch+response_version_set != 1:
        raise ValueError("One and only one of kwargs pre_launch, response_file "
                         "and response_version must be set.")
    # If pre_launch set, define response_version to 2.
    if pre_launch:
        response_version = 2
    # If response_file not set, define appropriate response file
    # based on version.
    if not response_file:
        try:
            response_filename = RESPONSE_VERSION_FILENAMES[str(response_version)]
        except KeyError:
            raise KeyError("Version number not recognized.")
        if response_version > 2:
            warnings.warn("Effective areas are not available (i.e. set  to zero).  "
                  "For response file versions > 2 time dependent effective "
                  "areas must be calculated via fitting, which is not supported "
                  "by this function at this time. "
                  "Version of this response file = {0}".format(response_version))
        # Define the directory in which the response file should exist
        # to be the sunpy download directory.
        config = sunpy.util.config.load_config()
        download_dir = config.get('downloads', 'download_dir')
        # Check response file exists in download_dir.  If not, download it.
        check_download_file(response_filename, IRIS_RESPONSE_REMOTE_PATH, download_dir,
                            replace=force_download)
        # Define response file as path + filename.
        response_file = os.path.join(download_dir, response_filename)

    # Read response file and store in a dictionary.
    raw_response_data = scipy.io.readsav(response_file)
    iris_response = dict([(name, raw_response_data["p0"][name][0])
                          for name in raw_response_data["p0"].dtype.names])
    # Convert some properties to more convenient types.
    iris_response["LAMBDA"] = Quantity(iris_response["LAMBDA"], unit=u.nm)
    iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)
    iris_response["GEOM_AREA"] = Quantity(iris_response["GEOM_AREA"], unit=u.cm**2)
    iris_response["VERSION"] = int(iris_response["VERSION"])
    # Convert some properties not found in version below version 3 to
    # more convenient types.
    if iris_response["VERSION"] > 2:
        # If DATE_OBS has a value, convert to datetime, else set to
        # None.
        try:
            iris_response["DATE_OBS"] = parse_time(iris_response["DATE_OBS"])
        except:
            iris_response["DATE_OBS"] = None
        # Convert C_F_TIME to array of datetime objects while
        # conserving shape.
        c_f_time = np.empty(iris_response["C_F_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_F_TIME"]):
            for j, t in enumerate(row):
                c_f_time[i][j] = parse_time(float(t))
        iris_response["C_F_TIME"] = c_f_time
        # Convert C_F_LAMBDA to Quantity.
        iris_response["C_F_LAMBDA"] = Quantity(iris_response["C_F_LAMBDA"], unit="nm")
        # Convert C_N_TIME to array of datetime objects while
        # conserving shape.
        c_n_time = np.empty(iris_response["C_N_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_N_TIME"]):
            for j, t in enumerate(row):
                c_n_time[i][j] = parse_time(float(t))
        iris_response["C_N_TIME"] = c_n_time
        # Convert C_N_LAMBDA to Quantity.
        iris_response["C_N_LAMBDA"] = Quantity(iris_response["C_N_LAMBDA"], unit="nm")
        # Convert C_S_TIME to array of datetime objects while
        # conserving shape.
        c_s_time = np.empty(iris_response["C_S_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_S_TIME"]):
            for j, column in enumerate(row):
                for k, t in enumerate(column):
                    c_s_time[i][j][k] = parse_time(float(t))
        iris_response["C_S_TIME"] = c_s_time
        # Convert DATE in ELEMENTS array to array of datetime objects.
        for i, t in enumerate(iris_response["ELEMENTS"]["DATE"]):
            iris_response["ELEMENTS"]["DATE"][i] = parse_time(t.decode())
        # Convert VERSION_DATE to datetime object.
        iris_response["VERSION_DATE"] = parse_time(iris_response["VERSION_DATE"].decode())
    else:
        # Change DATE tag in data with version < 2 to VERSION_DATE to
        # be consistent with more recent versions.
        iris_response["VERSION_DATE"] = datetime.datetime(int(iris_response["DATE"][0:4]),
                                                          int(iris_response["DATE"][4:6]),
                                                          int(iris_response["DATE"][6:8]))
        del(iris_response["DATE"])
    return iris_response


@custom_model
def _gaussian1d_on_linear_bg(x, amplitude=None, mean=None, standard_deviation=None,
                             constant_term=None, linear_term=None):
    return amplitude * np.exp(-((x - mean) / standard_deviation) ** 2) + constant_term + linear_term * x


def _calculate_orbital_wavelength_variation(data_array, date_data_created, slit_pixel_range=None,
                                            spline_smoothing=False, fit_individual_profiles=False,
                                            spacecraft_velocity=None, orbital_phase=None, roll_angle=None):
    """Calculates orbital corrections of spectral line positions using level 2 files.

    For data generated from the April 2014 pipeline, thermal and spacecraft velocity components
    have both been subtracted in the level 2 files.  Therefore, this routine calculates the
    residual orbital (thermal) variation.  For data generated from the Oct 2013 pipeline,
    this routine calculates the total of thermal and spacecraft velocity components.

    Parameters
    ----------
    data_array: `xarray.DataArray`
        IRIS spectrograph data from spectral window Mg II k 2796 as generated by
        `sunpy.spectra.sources.IRISRaster.`
    date_data_created: `datetime.datetime`
        Date the data was created by IRIS pipeline.  Used to determine where spacecraft
        velocity etc. needs to be accounted for.
    spacecraft_velocity: `astropy.units.quantity.Quantity`
        Velocity of spacecraft at each exposure in data_array.
        Must be set if date_data_created < 1 April 2014.
    orbital_phase: `numpy.array`
        Orbital phase of spacecraft at each exposure in data_array.  Available from
        auxiliary data in IRIS spectrograph fits files.
        Must be set if date_data_created < 1 April 2014.
    roll_angle: `astropy.units.quantity.Quantity`
        Roll angle of spacecraft. Must be set if date_data_created < 1 April 2014.

    Returns
    -------
    orbital_wavelength_variation: `astropy.table.Table`
        Contains the following columns:
        time: `datetime.datetime` objects
            Observation times of wavelength variations.
        FUV: `astropy.quantity.Quantity`
            Wavelength variation in the FUV.
        NUV: `astropy.quantity.Quantity`
            Wavelength variation in the NUV.

    """
    # Define vacuum rest wavelength of Ni I 2799 line.
    wavelength_nii = 2799.474 * u.Angstrom
    # Define factor converting NUV spectral pixel size to Angstrom
    specsize = 0.0255
    # Define date of new pipeline.
    date_new_pipeline = datetime.datetime(2014, 4, 1)
    if date_data_created < date_new_pipeline:
        # Check that there are measurement times with good values of
        # spacecraft velocity and orbital phase.
        bad_aux = np.asarray(np.isfinite(spacecraft_velocity) * np.isfinite(orbital_phase) * (-1), dtype=bool)
    # Generate wavelength vector containing only Ni I line.
    wavelength_window = Quantity(data_array.coords["wavelength"].values,
                                 unit=data_array.attrs["units"]["wavelength"])
    wavelength_roi_index = np.arange(len(wavelength_window))[
        np.logical_and(wavelength_window >= 2799.3 * u.Angstrom, wavelength_window <= 2799.8 * u.Angstrom)]
    # Check that there are at least 5 points in wavelength region.
    # Must have at least this many for a gaussian fit.
    if len(wavelength_roi_index) < 5:
        wavelength_roi_index = np.arange(5) + wavelength_roi_index[0]
    # Extract wavelength of region around Ni I line as array in units
    # of Angstroms.
    wavelength_roi = wavelength_window.to(u.Angstrom).value[wavelength_roi_index]
    # Keep only data within wavelength region of interest.
    data_array = data_array.isel(spectral_axis=slice(wavelength_roi_index[0], wavelength_roi_index[-1] + 1))
    # If user selected a sub-region of the slit, reduce data to just
    # that region.
    if slit_pixel_range:
        if len(slit_pixel_range) == 2:
            data_array = data_array.isel(slit_axis, slice(slit_pixel_range[0], slit_pixel_range[1]))
        else:
            raise TypeError("slit_pixel_range must be tuple of length 2 giving lower and " +
                            "upper bounds of section of slit over which to average line fits.")

    # Derive residual orbital variation.
    # Define array to hold averaged position of Ni I line at different
    # times.
    mean_line_wavelengths = np.empty(len(data_array.time)) * np.nan
    # Define initial guess for gaussian model.
    g_init = _gaussian1d_on_linear_bg(amplitude=-2., mean=wavelength_nii.value,
                                      standard_deviation=2., constant_term=50., linear_term=1.5)
    # Define fitting method.
    fit_g = fitting.LevMarLSQFitter()
    # Depending on user choice, either fit line as measured by each
    # pixel then average line position, or fit average line spectrum
    # from all slit pixels.
    if fit_individual_profiles:
        pixels_in_slit = len(raster.slit_axis)
        for k in range(len(raster.time)):
            pixel_line_wavelengths = np.empty(pixels_in_slit)*np.nan
            data_single_time = raster.isel(raster_axis=k)
            # Iterate through each pixel along slit and perform fit to
            # Ni I line.
            for j in range(2, pixels_in_slit-2):
                # Average over 5 pixels to improve signal-to-noise.
                intensity_mean_5pix = data_single_time.isel(slit_axis=slice(j-2, j+3)).mean(axis=0)
                # Fit gaussian to Ni I line.
                g = fit_g(g_init, wavelength_roi, intensity_mean_5pix)
                # Check that fit is within physically reasonable
                # limits.  If so, store line center wavelength in
                # mean_line_wavelengths array. Else leave element as
                # defined, i.e. NaN.
                if np.isfinite(g.amplitude) and g.amplitude < 0. and \
                            wavelength_roi[0] < g.mean < wavelength_roi[-1]:
                    pixel_line_wavelengths[j] = g.mean
            # Take average of Ni I line position from fits in each
            # pixel.
            mean_line_wavelengths[k] = np.nanmean(pixel_line_wavelengths)
    else:
        # Else average all line profiles then perform fit.
        # Iterate through each measurement time and fit a gaussian to
        # Ni I line.
        for k in range(len(raster.time)):
            # Get data averaged over slit.
            data_single_time = raster.isel(raster_axis=k)
            data_slit_averaged = data_single_time.to_masked_array().mean(axis=0).data
            # Fit Ni I line with a gaussian.
            # Perform fit.
            g = fit_g(g_init, wavelength_roi, data_slit_averaged)
            # Check that fit is within physically reasonable limits.
            # If so, store line center wavelength in
            # mean_line_wavelengths array. Else leave element as
            # defined, i.e. NaN.
            if np.isfinite(g.amplitude) and g.amplitude < 0. and \
                        wavelength_roi[0] < g.mean < wavelength_roi[-1]:
                mean_line_wavelengths[k] = g.mean
            # If data produced by old pipeline, subtract spacecraft velocity
            # from the line position.
            if date_created < date_new_pipeline:
                mean_line_wavelengths[k] = \
                    mean_line_wavelengths[k]-spacecraft_velocity[k]/3e8*wavelength_nii.to(u.Angstrom).value

    # Mark abnormal values.  Thermal drift is of the order of 2
    # unsummed wavelength pixels peak-to-peak.
    w_abnormal = np.where(np.abs(mean_line_wavelengths-np.nanmedian(mean_line_wavelengths)) >= specsize*2)[0]
    if len(w_abnormal) > 0:
        mean_line_wavelengths[w_abnormal] = np.nan
    # Further data reduction required for files from old pipeline.
    if date_created < date_new_pipeline:
        dw_th_A = mean_line_wavelengths - np.nanmean(mean_line_wavelengths)
        # Change the unit from Angstrom into unsummed wavelength pixel.
        dw_th_p = dw_th_A/specsize
        # Adjust reference wavelength using orbital phase information.
        if not(np.isfinite(orbital_phase)).all():
            warnings.warn("Orbital phase values are invalid.  Thermal drift may be offset by at most one pixel.")
            dw_th = dw_th
            # For absolute wavelength calibration of NUV, the
            # following amount (unit Angstrom) has to be
            # subtracted from the wavelengths.
            abswvl_nuv = np.nanmean(mean_line_wavelengths)-wavelength_nii.to(u.Angstrom).value
        else:
            # Define empirical sine fitting at 0 roll angle shifted by
            # different phase.
            sine_params = [-0.66615146, -1.0, 53.106583-roll_angle/360.*2*np.pi]
            phase_adj=np.nanmean(sine_params[0]*np.sin(sine_params[1]*orbital_phase+sine_params[2]))
            # thermal component of the orbital variation, in the unit of unsummed wavelength pixel
            dw_th=dw_th_p+phase_adj
            # For absolute wavelength calibration of NUV the following
            # amount (unit Angstrom) has to be subtracted from the
            # wavelengths.
            abswvl_nuv = np.nanmean(mean_line_wavelengths)-wavelength_nii.to(u.Angstrom).value-phase_adj*specsize
    else:
        # Calculate relative variation of the line position.
        dw_th = mean_line_wavelengths-np.nanmean(mean_line_wavelengths)

    # If spline_smoothing=True, perform spline fit a smoothing to
    # eliminate the 5 minute photospheric oscillation.
    if spline_smoothing:
        # Define spacing of spline knots in seconds.
        spline_knot_spacing = 300.
        # Create array of time in seconds from first time and
        # calculate duration of fitting period.
        time_s = np.asarray(x.coords["time"]-x.coords["time"][0], dtype=float)/1e9
        duration = time_s[-1]-time_s[0]
        # Check whether there is enough good data for a spline fit.
        if duration < spline_knot_spacing:
            raise ValueError("Not enough data for spline fit.")
        # Check whether there is enough good data for a spline fit.
        wgood = np.isfinite(mean_line_wavelengths)
        ngood = float(sum(wgood))
        wbad = not(np.isfinite(mean_line_wavelengths))
        nbad = float(sum(wbad))
        if nbad/ngood > 0.25:
            raise ValuError("Not enough good data for spline fit.")
        # Smooth residual thermal variation curve to eliminate the
        # 5-min photospheric oscillation.
        # Determine number of smoothing point using 3 point
        # lagrangian derivative.
        deriv_time = np.array([(time_s[i+1]-time_s[i-1])/2. for i in range(1,len(time_s)-1)])
        deriv_time = np.insert(deriv_time, 0, (-3*time_s[0]+4*time_s[1]-time_s[2])/2)
        deriv_time = np.insert(deriv_time, -1, (3*time_s[-1]-4*time_s[-2]+time_s[-3])/2)
        n_smooth = int(spline_knot_spacing/deriv_time.mean())
        if n_smooth < len(wgood):
            dw_good = convolve(dw_th[good], Box1DKernel(n_smooth))
        else:
            dw_good = dw_th[good]
        time_good = time_s[good]
        # Fit spline.
        tck = interpolate.splrep(time_good, dw_good, s=0)
        dw_th = interpolate.splev(time_s, tck)

    # Derive residual orbital curves in FUV and NUV and store
    # in a table.
    times = [datetime.datetime.utcfromtimestamp(t/1e9) for t in raster.coords["time"].values.tolist()]
    # Depeding on which pipeline produced the files...
    if date_created < date_new_pipeline:
        dw_orb_fuv = dw_th * (-0.013) + spacecraft_velocity.to(u.km/u.s).value / (3.e5) * 1370. * u.Angstrom
        dw_orb_nuv = dw_th * 0.0255 + spacecraft_velocity.to(u.km/u.s).value / (3.e5) * 2800. * u.Angstrom
    else:
        dw_orb_fuv = dw_th*(-1)*u.Angstrom
        dw_orb_nuv = dw_th*u.Angstrom

    orbital_wavelength_variation = Table([times, dw_orb_fuv, dw_orb_nuv],
                                         names=("time", "wavelength variation FUV", "wavelength variation NUV"))
    return orbital_wavelength_variation

def get_detector_type(meta):
    """
    Gets the IRIS detector type from a meta dictionary.

    In this function, FUV1 and FUV2 are just assigned as FUV.

    Parameters
    ----------
    meta: dict-like
        Dictionary-like object containing entry for "detector type"

    Returns
    -------
    detector_type: `str`
       Detector type.

    """
    if "FUV" in meta["detector type"]:
        detector_type = "FUV"
    else:
        detector_type = meta["detector type"]
    return detector_type

def convert_between_DN_and_photons(old_data_arrays, old_unit, new_unit):
    """Converts arrays from IRIS DN to photons or vice versa.

    In this function, an inverse time component due to exposure time
    correction is ignored during calculations but preserved in final unit.

    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.

    old_unit: `astropy.unit.Unit`
        Unit of data arrays.

    new_unit: `astropy.unit.Unit`
        Unit to convert data arrays to.

    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays converted to new_unit.

    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays with any inverse time component preserved.

    """
    if old_unit == new_unit or old_unit == new_unit / u.s:
        new_data_arrays = [data for data in old_data_arrays]
        new_unit_time_accounted = old_unit
    else:
        # During calculations, the time component due to exposure
        # time correction, if it has been applied, is ignored.
        # Check here whether the time correction is present in the
        # original unit so that is carried through to new unit.
        if u.s not in (old_unit * u.s).decompose().bases:
            old_unit_without_time = old_unit * u.s
            new_unit_time_accounted = new_unit / u.s
        else:
            old_unit_without_time = old_unit
            new_unit_time_accounted = new_unit
        # Convert data and uncertainty to new unit.
        new_data_arrays = [(data * old_unit_without_time).to(new_unit).value
                           for data in old_data_arrays]
    return new_data_arrays, new_unit_time_accounted

def calculate_exposure_time_correction(old_data_arrays, old_unit, exposure_time,
                                       force=False):
    """
    Applies exposure time correction to data arrays.

    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.

    old_unit: `astropy.unit.Unit`
        Unit of data arrays.

    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays with exposure time corrected for.

    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays after exposure time correction.

    """
    # If force is not set to True and unit already includes inverse time,
    # raise error as exposure time correction has probably already been
    # applied and should not be applied again.
    if force is not True and u.s in old_unit.decompose().bases:
        raise ValueError(APPLY_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does not include inverse time and so
        # exposure does need to be applied, or
        # user has set force=True and wants the correction applied
        # regardless of the unit.
        new_data_arrays = [old_data/exposure_time for old_data in old_data_arrays]
        new_unit = old_unit/u.s
    return new_data_arrays, new_unit

def uncalculate_exposure_time_correction(old_data_arrays, old_unit, exposure_time,
                                         force=False):
    """
    Removes exposure time correction from data arrays.

    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`s
        Arrays of data to be converted.

    old_unit: `astropy.unit.Unit`
        Unit of data arrays.

    exposure_time: `numpy.ndarray`
        Exposure time in seconds for each exposure in data arrays.

    Returns
    -------
    new_data_arrays: `list` of `numpy.ndarray`s
        Data arrays with exposure time correction removed.

    new_unit_time_accounted: `astropy.unit.Unit`
        Unit of new data arrays after exposure time correction removed.

    """
    # If force is not set to True and unit does not include inverse time,
    # raise error as exposure time correction has probably already been
    # undone and should not be undone again.
    if force is not True and u.s in (old_unit*u.s).decompose().bases:
        raise ValueError(UNDO_EXPOSURE_TIME_ERROR)
    else:
        # Else, either unit does include inverse time and so
        # exposure does need to be removed, or
        # user has set force=True and wants the correction removed
        # regardless of the unit.
        new_data_arrays = [old_data * exposure_time for old_data in old_data_arrays]
        new_unit = old_unit*u.s
    return new_data_arrays, new_unit

def convert_or_undo_photons_per_sec_to_radiance(
        data_quantities, obs_wavelength, detector_type,
        spectral_dispersion_per_pixel, solid_angle, undo=False):
    """
    Converts data quantities from counts/s to radiance (or vice versa).

    Parameters
    ----------
    data_quantities: iterable of `astropy.units.Quantity`s
        Quantities to be converted.  Must have units of counts/s or
        radiance equivalent counts, e.g. erg / cm**2 / s / sr / Angstrom.

    obs_wavelength: `astropy.units.Quantity`
        Wavelength at each element along spectral axis of data quantities.

    detector_type: `str`
        Detector type: 'FUV', 'NUV', or 'SJI'.

    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.

    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    undo: `bool`
        If False, converts counts/s to radiance.
        If True, converts radiance to counts/s.
        Default=False

    Returns
    -------
    new_data_quantities: `list` of `astropy.units.Quantity`s
        Data quantities converted to radiance or counts/s
        depending on value of undo kwarg.

    """
    # Check data quantities are in the right units.
    if undo is True:
        for i, data in enumerate(data_quantities):
            if not data.unit.is_equivalent(RADIANCE_UNIT):
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=True, "
                    "unit must be equivalent to {0}.  Error found for {1}th element "
                    "of data_quantities. Unit: {2}".format(RADIANCE_UNIT, i, data.unit))
    else:
        for data in data_quantities:
            if data.unit != u.photon/u.s:
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=False, "
                    "unit must be equivalent to {0}.  Error found for {1}th element "
                    "of data_quantities. Unit: {2}".format(u.photon/u.s, i, data.unit))
    photons_per_sec_to_radiance_factor = calculate_photons_per_sec_to_radiance_factor(
        obs_wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle)
    # Change shape of arrays so they are compatible for broadcasting
    # with data and uncertainty arrays.
    photons_per_sec_to_radiance_factor = \
        _reshape_1D_wavelength_dimensions_for_broadcast(photons_per_sec_to_radiance_factor,
                                                        data_quantities[0].ndim)
    # Perform (or undo) radiometric conversion.
    if undo is True:
        new_data_quantities = [(data / photons_per_sec_to_radiance_factor).to(u.photon/u.s)
                               for data in data_quantities]
    else:
        new_data_quantities = [(data*photons_per_sec_to_radiance_factor).to(RADIANCE_UNIT)
                               for data in data_quantities]
    return new_data_quantities

def calculate_photons_per_sec_to_radiance_factor(
        wavelength, detector_type, spectral_dispersion_per_pixel, solid_angle):
    """
    Calculates multiplicative factor that converts counts/s to radiance for given wavelengths.

    Parameters
    ----------
    wavelength: `astropy.units.Quantity`
        Wavelengths for which counts/s-to-radiance factor is to be calculated

    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.

    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.

    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    Returns
    -------
    radiance_factor: `astropy.units.Quantity`
        Mutliplicative conversion factor from counts/s to radiance units
        for input wavelengths.

    """
    # Get effective area and interpolate to observed wavelength grid.
    eff_area_interp = _get_interpolated_effective_area(detector_type, wavelength)
    # Return radiometric conversed data assuming input data is in units of photons/s.
    return constants.h * constants.c / wavelength / u.photon / \
           spectral_dispersion_per_pixel / eff_area_interp / solid_angle

def _get_interpolated_effective_area(detector_type, obs_wavelength):
    # Get effective area
    ########### This needs to be generalized to the time of OBS once that functionality is written #########
    iris_response = get_iris_response(pre_launch=True)
    if detector_type == "FUV":
        detector_type_index = 0
    elif detector_type == "NUV":
        detector_type_index = 1
    else:
        raise ValueError("Detector type not recognized.")
    eff_area = iris_response["AREA_SG"][detector_type_index, :]
    response_wavelength = iris_response["LAMBDA"]
    # Interpolate the effective areas to cover the wavelengths
    # at which the data is recorded:
    eff_area_interp_base_unit = u.Angstrom
    tck = interpolate.splrep(response_wavelength.to(eff_area_interp_base_unit).value,
                             eff_area.to(eff_area_interp_base_unit ** 2).value, s=0)
    eff_area_interp = interpolate.splev(
        obs_wavelength.to(eff_area_interp_base_unit).value, tck) * eff_area_interp_base_unit ** 2
    return eff_area_interp

def _reshape_1D_wavelength_dimensions_for_broadcast(wavelength, n_data_dim):
    if n_data_dim == 1:
        pass
    elif n_data_dim == 2:
        wavelength = wavelength[np.newaxis, :]
    elif n_data_dim == 3:
        wavelength = wavelength[np.newaxis, np.newaxis, :]
    else:
        raise ValueError("IRISSpectrogram dimensions must be 2 or 3.")
    return wavelength

def _convert_iris_sequence(sequence, new_unit):
    """Converts data and uncertainty in an IRISSpectrogramSequence between units.

    Parameters
    ----------
    sequence: `NDCubeSequence`, `SpectrogramSequence` or `IRISSpectrogramSequence`
        Sequence whose constituent NDCubes are be converted to new units.

    new_unit: `astropy.units.Unit` or `str`
       Unit to which the data is to be converted.

    Returns
    -------
    converted_data_list: `list` of `NDCube`s.
       List of NDCubes with data and uncertainty attributes converted to new_unit.

    """
    # Define empty list to hold NDCubes with converted data and uncertainty.
    converted_data_list = []
    # Cycle through each NDCube, convert data and uncertainty to new
    # units, and append to list.
    for i, cube in enumerate(sequence.data):
        # Determine what type of DN unit is needed based on detector type.
        detector_type = _get_detector_type(cube.meta)
        if new_unit == "DN":
            new_unit = DN_UNIT[detector_type]
        # If NDCube is already in new unit, add NDCube as is to list.
        if cube.unit is new_unit or cube.unit is new_unit / u.s:
            converted_data_list.append(cube)
        # Else convert data and uncertainty to new unit.
        if cube.unit != new_unit or cube.unit != new_unit / u.s:
            # During calculations, the time component due to exposure
            # time correction, if it has been applied, is ignored.
            # Check here whether the time correction is present in the
            # original unit so that is carried through to new unit.
            if u.s not in (cube.unit.decompose() * u.s).bases:
                new_unit_time_accounted = new_unit / u.s
            else:
                new_unit_time_accounted = new_unit
            # Convert data and uncertainty to new unit.
            data = (cube.data * cube.unit).to(new_unit).value
            uncertainty = (cube.uncertainty.array * cube.unit).to(new_unit).value
            # Append new instance of NDCube in new unit to list.
            converted_data_list.append(NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask,
                unit=new_unit_time_accounted, uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords)))
    return converted_data_list

def _apply_or_undo_exposure_time_correction(sequence, correction_function):
    """Applies or undoes exposure time correction to a sequence of NDCubes.

    Correction is applied (or undone) to both data and uncertainty attributes of NDCubes.

    Parameters
    ----------
    sequence: `NDCubeSequence`, `SpectrogramSequence` or `IRISSpectrogramSequence`
        Sequence whose constituent NDCubes are be converted to new units.
        NDCubes with sequence must have an 'exposure time' entry in its extra
        coords attribute.

    correction_function: function
        Function applying or undoing exposure time correction.

    Returns
    -------
    converted_data_list: `list` of `NDCube`s.
       List of NDCubes with data and uncertainty corrected (or uncorrected)
       for exposure time.

    """
    converted_data_list = []
    for i, cube in enumerate(sequence.data):
        if u.s not in cube.unit.decompose().bases:
            exposure_time_s = cube._extra_coords["exposure time"]["value"].to(u.s).value
            if len(cube.dimensions.shape) == 1:
                pass
            elif len(cube.dimensions.shape) == 2:
                exposure_time_s = exposure_time_s[:, np.newaxis]
            elif len(cube.dimensions.shape) == 3:
                exposure_time_s = exposure_time_s[:, np.newaxis, np.newaxis]
            else:
                raise ValueError("NDCube dimensions must be 2 or 3. Dimensions={0}".format(
                    len(cube.dimensions.shape)))
            data = correction_function(cube.data, exposure_time_s)
            uncertainty = correction_function(cube.uncertainty.array, exposure_time_s)
            converted_data_list.append(NDCube(
                data, wcs=cube.wcs, meta=cube.meta, mask=cube.mask, unit=cube.unit / u.s,
                uncertainty=uncertainty,
                extra_coords=_extra_coords_to_input_format(cube._extra_coords)))
        else:
            converted_data_list.append(cube)
    return converted_data_list
