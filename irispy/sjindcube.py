# -*- coding: utf-8 -*-

from datetime import timedelta

import astropy.units as u
from astropy.io import fits as pyfits
from ndcube import NDCube
from ndcube.utils.wcs import WCS
import numpy as np
from sunpy.time import parse_time

from irispy import iris_tools

__all__ = ['SJI_NDCube']

# the following value is only appropriate for byte scaled images
BAD_PIXEL_VALUE = -200.
# the following value is only appropriate for unscaled images
BAD_PIXEL_VALUE_UNSCALED = -32768


class SJI_NDCube(NDCube):

    """
    SJICube

    Class representing SJI Images described by a single WCS

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.

    wcs: `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information

    unit : `astropy.unit.Unit` or `str`
        Unit for the dataset.
        Strings that can be converted to a Unit are allowed.

    meta : dict-like object
        Additional meta information about the dataset.

    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn’t mandatory. If the
        uncertainty has no such attribute the uncertainty is stored as
        UnknownUncertainty.
        Defaults to None.

    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.

    extra_coords : iterable of `tuple`s, each with three entries
        (`str`, `int`, `astropy.units.quantity` or array-like)
        Gives the name, axis of data, and values of coordinates of a data axis
        not included in the WCS object.

    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every
        attribute before saving it while False tries to save every parameter
        as reference. Note however that it is not always possible to save the
        input as reference.
        Default is False.

    Examples
    --------
    >>> from irispy.sjindcube import read_fits
    >>> from irispy.data import sample
    >>> sji = read_fits(sample.SJI_CUBE_1400)
    """

    def __init__(self, data, wcs, uncertainty=None, unit=None, meta=None,
                 extra_coords=None, mask=None, copy=False, missing_axis=None):

        """
        Initialization of Slit Jaw Imager
        """

        # Initialize SJI_NDCube.
        super(SJI_NDCube, self).__init__(data, wcs, uncertainty=uncertainty,
                                         mask=mask, meta=meta, unit=unit,
                                         extra_coords=extra_coords, copy=copy,
                                         missing_axis=missing_axis)

    def __repr__(self):
        return (
            """SJI_NDCube {dtype!s}
    ---------
    Observatory:\t {obs}
    Instrument:\t {inst}
    Detector:\t {det}
    Measurement:\t {meas}
    Wavelength:\t {wave}
    Obs. Start:\t {date_start:{tmf}}
    Obs. End:\t {date_end:{tmf}}
    Num. of Frames:\t {frame_num}
    IRIS Obs. id:\t {obs_id}
    IRIS Obs. Description:\t {obs_desc}
    Dimensions:\t {dim}
    Scale:\t\t {scale}
    """.format(dtype=self.__class__.__name__,
               obs=self.observatory, inst=self.instrument,
               det=self.detector, meas=self.measurement,
               wave=self.wavelength, date_start=self.date[0],
               date_end=self.date[-1], frame_num=len(self),
               dim=u.Quantity(self.dimensions), scale=u.Quantity(self.scale),
               obs_id=self.iris_obs_id, obs_desc=self.iris_obs_description,
               tmf=TIME_FORMAT) + self.data.__repr__())


def read_fits(fitsfile):
    """
    Read IRIS level 2 SJI FITS from an OBS into an SJI_NDCube instance

    Parameters
    ----------
    fitsfile : `str`
        The fits file path

    Returns
    -------
    result: 'irispy.sjindcube.SJI_NDCube'

    """

    # Open a fits file
    my_file = pyfits.open(fitsfile)
    # Derive WCS, data and mask for NDCube from fits file.
    wcs = WCS(my_file[0].header)
    data = my_file[0].data
    data_nan_masked = my_file[0].data
    data_nan_masked[data == BAD_PIXEL_VALUE] = np.nan
    mask = data_nan_masked == BAD_PIXEL_VALUE
    # Derive unit and readout noise from detector.
    exposure_times = my_file[1].data[:, my_file[1].header["EXPTIMES"]]
    unit = iris_tools.DN_UNIT["SJI"]
    readout_noise = iris_tools.READOUT_NOISE["SJI"]
    # Derive uncertainty of data for NDCube from fits file.
    uncertainty = u.Quantity(np.sqrt((data_nan_masked*unit).to(u.photon).value
                                     + readout_noise.to(u.photon).value**2),
                             unit=u.photon).to(unit).value
    # Derive extra coordinates for NDCube from fits file.
    times = np.array([parse_time(my_file[0].header["STARTOBS"])
                      + timedelta(seconds=s)
                      for s in my_file[1].data[:, my_file[1].header["TIME"]]])
    extra_coords = [('TIME', 0, times), ("EXPOSURE TIME", 0, exposure_times)]
    # Extraction of meta for NDCube from fits file.
    number_of_images = my_file[0].data.shape[0]
    try:
        date_obs = parse_time(my_file[0].header["DATE_OBS"])
    except ValueError:
        date_obs = None
    try:
        date_end = parse_time(my_file[0].header["DATE_END"])
    except ValueError:
        date_end = None
    #for i in range(number_of_images):
        meta = {'WAVElNTH': my_file[0].header.get('TWAVE1'),
                'DETECTOR': my_file[0].header.get('INSTRUME'),
                'WAVEUNIT': "Angstrom",
                'OBSRVTRY': my_file[0].header.get('TELESCOP'),
                'DATE_OBS': date_obs,
                'DATE_END': date_end,}
    #    for item in my_file[1].header[7:]:
    #        meta[item] = my_file[1].data[i, my_file[1].header[item]]
    #        if item.count('EXPTIMES'):
    #            meta['EXPTIME'] = my_file[1].data[i, my_file[1].header[item]]


    my_file.close()

    return SJI_NDCube(data_nan_masked, wcs, uncertainty=uncertainty, unit=unit,
                      mask=mask, extra_coords=extra_coords)
    #return data_nan_masked
    #return extra_coords
    #return NDCube(data, wcs)
    #return mask
    #return meta
