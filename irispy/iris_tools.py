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