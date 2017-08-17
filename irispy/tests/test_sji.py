# -*- coding: utf-8 -*-
# """Tests for functions in sji.py"""
import pytest
import astropy.units as u
from irispy.sji import SJICube, SJIMap

iris_dir = '/Users/schriste/Developer/repositories/sample-data/irispy/'
f1330 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits'
f1400 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits'
f2796 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits'
f2832 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits'
all_sji_files = [f1330, f1400, f2796, f2832]


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_open_cubes(sji_file):
    """test that SJICube open lvl2 sji files"""
    assert isinstance(SJICube(sji_file), SJICube)


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_observatory_meta(sji_file):
    assert SJICube(sji_file).observatory == 'IRIS'


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_detector_meta(sji_file):
    assert SJICube(sji_file).detector == 'SJI'


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_measurement_is_quantity(sji_file):
    assert isinstance(SJICube(sji_file).measurement, u.Quantity)


@pytest.mark.skip
@pytest.mark.parametrize("sji_file, measurement", [(all_sji_files[0], 1330),
                                                   (all_sji_files[1], 1400),
                                                   (all_sji_files[2], 2796),
                                                   (all_sji_files[3], 2832)])
def test_cubes_measurement_value_is_correct(sji_file, measurement):
    assert SJICube(sji_file).measurement.value == measurement


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_length(sji_file):
    assert len(SJICube(sji_file)) == SJICube(sji_file).data.shape[0]


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_length(sji_file):
    assert len(SJICube(sji_file)) == SJICube(sji_file).data.shape[0]


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_indexing_to_map(sji_file):
    assert isinstance(SJICube(sji_file)[0], SJIMap)


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_cubes_indexing_to_SJICube(sji_file):
    assert isinstance(SJICube(sji_file)[0:4], SJICube)


@pytest.mark.skip
@pytest.mark.parametrize("sji_file", all_sji_files)
def test_making_SJIMap(sji_file):
    sji_cube = SJICube(sji_file)
    header = sji_cube.meta(0)
    data = sji_cube.data[0,:,:]
    assert isinstance(SJIMap(header, data), SJICube)


        #def singlefile():
#    mc=sji.SJI_fits_to_cube(fname[0],0,10)
#    mc_db = sji.dustbuster(mc)
#    print('Single file input Passed')

#def listfits():
#    mc=sji.SJI_fits_to_cube(fname)
#    print('Multiple file input Passed')
    # Bad Pixel Map (BPM) before dustbuster
#    m1 = ma.masked_inside(mc[0].data, -199, .1)

#    mc_db=sji.dustbuster(mc)
    # Fail condition: If mask finds unchanged values in BPM will return 'numpy.array'
    # Pass condition: If mask finds no unchanged values in BPM will return 'numpy.bool'
#    fail_count=0
#    for i, map in enumerate(mc_db):
#        m = ma.masked_inside(map.data, -199, .1)
#        if type(m.mask)==type(m1.mask):
#            print('DUSTBUSTER Failed: ' +str(i))
#            fail_count+=1
#    if fail_count==0:
#        print('DUSTBUSTER Passed')


#test1=singlefile()
#test2=listfits()
