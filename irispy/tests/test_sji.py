# -*- coding: utf-8 -*-
# """Tests for functions in sji.py"""
import pytest
import numpy.ma as ma
from irispy.sji import SJICube, SJIMap

iris_dir = '/Users/schriste/Developer/repositories/sample-data/irispy/'
f1330 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits'
f1400 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits'
f2796 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits'
f2832 = iris_dir + 'iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits'
all_sji_files = [f1330, f1400, f2796, f2832]

fname = [iris_dir + 'iris_l2_20170228_053417_3600106076_SJI_1330_t000.fits',
              iris_dir + 'iris_l2_20170228_082338_3600106076_SJI_1330_t000.fits',
              iris_dir + 'iris_l2_20170228_104517_3600106076_SJI_1330_t000.fits']


@pytest.mark.parametrize("sji_file", all_sji_files)
def test_open_cubes(sji_file):
    assert isinstance(SJICube(sji_file), SJICube)


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

