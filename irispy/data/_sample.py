# -*- coding: utf-8 -*-
"""IrisPy Sample and test files"""

import os.path
import socket
import warnings
from zipfile import ZipFile
from shutil import move

from astropy.utils.data import download_file

from sunpy.util.net import url_exists
from sunpy.util.config import get_and_create_sample_dir
from sunpy import config

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"

# Creating the directory for sample files to be downloaded
sampledata_dir = get_and_create_sample_dir()

# urls to search for the data
_base_urls = (
    'http://data.sunpy.org/sample-data/',
    'http://hesperia.gsfc.nasa.gov/~schriste/sunpy-sample-data/')

_files = {"SJI_CUBE_2832": ("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits", ""),
          "SJI_CUBE_1330": ("iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits", ""),
          "SJI_CUBE_2796": ("iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits", ""),
          "SJI_CUBE_1400": ("iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits", ""),
          "RASTER": ("iris_l2_20170502_052551_3893010094_raster_t000_r00000.fits", ".tar.gz")
         }

sample_files = {}
for key in _files:
    sample_files[key] = os.path.abspath(os.path.join(sampledata_dir, _files[key][0]))


def download_data(progress=True, overwrite=True, timeout=None):
    """
    Download the sample/test data.
    Parameters
    ----------
    progress: `bool`
        Show a progress bar during download
    overwrite: `bool`
        If exist overwrites the downloaded sample data.
    timeout: `float`
        The timeout in seconds. If `None` the default timeout is used from
        `astropy.utils.data.Conf.remote_timeout`.
    Returns
    -------
    None
    """

    print(sampledata_dir)
    number_of_files_fetched = 0
    print("Downloading sample files to {}".format(sampledata_dir))
    for handle, file_name in _files.items():
        if not overwrite:
            if os.path.isfile(os.path.join(sampledata_dir,
                                           file_name[0])):
                number_of_files_fetched += 1
                continue

        for base_url in _base_urls:
            full_file_name = file_name[0] + file_name[1]
            try:
                exists = url_exists(os.path.join(base_url, full_file_name))
                if exists:
                    f = download_file(os.path.join(base_url, full_file_name))
                    real_name, ext = os.path.splitext(full_file_name)

                    if file_name[1] == '.zip':
                        print("Unpacking: {}".format(real_name))
                        with ZipFile(f, 'r') as zip_file:
                            zip_file.extract(real_name, sampledata_dir)
                        os.remove(f)
                    else:
                        # move files to the data directory
                        move(f, os.path.join(sampledata_dir, file_name[0]))
                    # increment the number of files obtained to check later
                    number_of_files_fetched += 1
                    break
            except (socket.error, socket.timeout) as e:
                warnings.warn("Download failed with error {}. \n"
                              "Retrying with different mirror.".format(e))

    if number_of_files_fetched < len(list(_files.keys())):
        raise URLError("Could not download all samples files."
                       "Problem with accessing sample data servers.")