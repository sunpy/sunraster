# -*- coding: utf-8 -*-
# """Tests for functions in sji.py"""
import pytest
import astropy.units as u
from irispy.new_sji import SJICube, read_iris_sji_level2_fits

#iris_dir = '/Users/schriste/Developer/repositories/sample-data/irispy/'
iris_dir = '/mn/stornext/u3/baptistp/iris/level2/2017/05/02/20170502_052551_3893010094/'
f1330 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits'
f1400 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits'
f2796 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits'
f2832 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits'
all_sji_files = [f1330, f1400, f2796, f2832]



@pytest.mark.parametrize("sji_file", all_sji_files)
def test_open_cubes(sji_file):
    """test that SJICube open lvl2 sji files"""
    assert isinstance(read_iris_sji_level2_fits(sji_file), SJICube)



@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_observatory_meta(sji_file):
    assert read_iris_sji_level2_fits(sji_file).meta['TELESCOP'] == 'IRIS'



@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_detector_meta(sji_file):
    assert read_iris_sji_level2_fits(sji_file).meta['INSTRUME'] == 'SJI'




@pytest.mark.parametrize("sji_file, measurement", [(all_sji_files[0], 1330),
                                                   (all_sji_files[1], 1400),
                                                   (all_sji_files[2], 2796),
                                                   (all_sji_files[3], 2832)])
def test_cubes_measurement_value_is_correct(sji_file, measurement):
    assert read_iris_sji_level2_fits(sji_file).meta['TWAVE1'] == measurement



@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_indexing_to_SJICube(sji_file):
    assert isinstance(read_iris_sji_level2_fits(sji_file)[0:4], SJICube)
