# -*- coding: utf-8 -*-
"""IRISPy sample data files"""
from __future__ import absolute_import

import sys
from sunpy.data.sample import get_sample_file

_base_urls = (
    'http://data.sunpy.org/sample-data/irispy/',
    'https://github.com/sunpy/sample-data/raw/master/irispy/'
)

_sample_files = {"SJI_CUBE_2832": "iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits",
                 "SJI_CUBE_1330": "iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits",
                 "SJI_CUBE_2796": "iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits",
                 "SJI_CUBE_1400": "iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits",
                 "RASTER": "iris_l2_20170502_052551_3893010094_raster.fits.tar.zip"
}

file_list = []
file_dict = {}
for _key in _sample_files:
    f = get_sample_file(_sample_files[_key], _base_urls)
    setattr(sys.modules[__name__], _key, f)
    file_list.append(f)
    file_dict.update({_key: f})

__all__ = list(_sample_files.keys()) + ['file_dict', 'file_list']