"""IRIS instrument tools."""

import datetime
import warnings
import os.path

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.units.quantity import Quantity
from astropy import wcs
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import scipy.io
from sunpy.time import parse_time

# Define some properties of IRIS detectors.  Source: IRIS instrument
# paper.
DETECTOR_GAIN = {"NUV": 18., "FUV": 6., "SJI": 18.}
DETECTOR_YIELD = {"NUV": 1., "FUV": 1.5, "SJI": 1.}
READOUT_NOISE = {"NUV": {"value": 1.2, "unit": "DN"}, "FUV": {"value": 3.1, "unit": "DN"},
                 "SJI": {"value": 1.2, "unit": "DN"}}


def convert_DN_to_photons(data, detector_type):
    """Converts IRIS data number to photon counts depending on which CCD is being used.

    Parameters
    ----------
    data: array-like
        IRIS data in units of DN.
    detector_type: `str`
        Detector/CCD upon which data was recorded.
        Can take values 'NUV', 'FUV' or 'SJI'.

    Returns
    -------
    data_photons: array-like
        IRIS data in units of photon counts.

    """
    return DETECTOR_GAIN[detector_type]/DETECTOR_YIELD[detector_type]*data


def convert_photons_to_DN(data, detector_type):
    """Converts photon counts to IRIS data number depending on which CCD is being used.

    Parameters
    ----------
    data: array-like
        IRIS data in units of photon counts.
    detector_type: `str`
        Detector/CCD upon which data was recorded.
        Can take values 'NUV', 'FUV' or 'SJI'.

    Returns
    -------
    data_dn: array-like
        IRIS data in units of data number.

    """
    return DETECTOR_YIELD[detector_type]/DETECTOR_GAIN[detector_type]*data


def calculate_intensity_fractional_uncertainty(data, data_unit, detector_type):
    """Calculates fractional uncertainty of IRIS data.

    Parameters
    ----------
    data: array-like
        IRIS data.
    data_unit: `str`
        Unit of data.  Must be either 'DN' or 'photons'.
    detector_type: `str`
        Detector/CCD upon which data was recorded.
        Can take values 'NUV', 'FUV' or 'SJI'.

    Returns
    -------
    fractional_uncertainty: array-like
        Fractional uncertainty of each element in data array.
        Same shape as data.

    """
    photons_per_dn = DETECTOR_GAIN[detector_type]/DETECTOR_YIELD[detector_type]
    if data_unit == "DN":
        intensity_ph = photons_per_dn*convert_DN_to_photons(data, detector_type)
    elif data_unit == "photons":
        intensity_ph = data
    else:
        raise ValueError("Data not in recognized units: {0}".format(data_unit))
    readout_noise_ph = READOUT_NOISE[detector_type]["value"]*photons_per_dn
    uncertainty_ph = np.sqrt(intensity_ph+readout_noise_ph**2.)
    return uncertainty_ph/intensity_ph


def get_iris_response(pre_launch=False, response_file=None, response_version=None):
    """Returns IRIS response structure.

    One and only one of pre_launch, response_file and response_version must be set.

    Parameters
    ----------
    pre_launch: `bool`
        Equivalent to setting response_version=2.  Cannot be set simultaneously
        with response_file kwarg. Default=False
    response_file: `int`
        Version number of effective area file to be used.  Cannot be set simultaneously with
        pre_launch kwarg.  Default=latest
    response_version : `int`
        Version number of effective area file to be used. Cannot be set simultaneously with
        response_file or pre_launch kwarg. Default=latest

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
        respdir = os.path.expanduser(os.path.join("~", "ssw", "iris", "response"))
        if response_version == 1:
            response_file = os.path.join(respdir, "iris_sra_20130211.geny")
        elif response_version == 2:
            response_file = os.path.join(respdir, "iris_sra_20130715.geny")
        elif response_version == 3:
            response_file = os.path.join(respdir, "iris_sra_c_20150331.geny")
            warnings.warn("Effective areas are not available (i.e. set  to zero).  "
                          "For response file versions > 2 time dependent effective areas must be"
                          " calculated via fitting, which is not supported by this function at"
                          " this time. Version of this response file = {0}".format(response_version))
        else:
            raise ValueError("Version number not recognized.")
    # Read response file and store in a dictionary.
    raw_response_data = scipy.io.readsav(response_file)
    iris_response = dict([(name, raw_response_data["p0"][name][0])
                          for name in raw_response_data["p0"].dtype.names])
    # Convert some properties to more convenient types.
    iris_response["LAMBDA"] = Quantity(iris_response["LAMBDA"], unit="nm")
    iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit="cm")
    iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit="cm")
    iris_response["GEOM_AREA"] = Quantity(iris_response["GEOM_AREA"], unit="cm")
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