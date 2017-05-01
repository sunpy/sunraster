import sji
import matplotlib.pyplot as plt
import numpy.ma as ma

iris_dir = '/Users/shelbe/sunpy/data/sample_data/'
fname = [iris_dir + 'iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits',
              iris_dir + 'iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits',
              iris_dir + 'iris_l2_20130801_074720_4040000014_SJI_1400_t000.fits']

def singlefile():
    mc=sji.SJI_fits_to_cube(fname[0])
    mc.peek()
    plt.show()
    print('Single file input Passed')

def listfits():
    mc=sji.SJI_fits_to_cube(fname)
    print('Multiple file input Passed')
    # Bad Pixel Map (BPM) before dustbuster
    m1 = ma.masked_inside(mc[0].data, -199, .1)

    mc_db=sji.dustbuster(mc)
    # Fail condition: If mask finds unchanged values in BPM will return 'numpy.array'
    # Pass condition: If mask finds no unchanged values in BPM will return 'numpy.bool'
    fail_count=0
    for i, map in enumerate(mc_db):
        m = ma.masked_inside(map.data, -199, .1)
        if type(m.mask)==type(m1.mask):
            print('DUSTBUSTER Failed: ' +str(i))
            fail_count+=1
    if fail_count==0:
        print('DUSTBUSTER Passed')


test1=singlefile()
test2=listfits()

