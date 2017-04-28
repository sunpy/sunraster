import datetime
import sunpy.cm as cm
import sunpy.instr.iris
import sunpy.map
import numpy as np
import numpy.ma as ma
import matplotlib.colors as colors
from astropy.io import fits
import sunpy.io
import sunpy.time


'''
This module provides movie tools for level 2 IRIS SJI fits file
'''

__all__ = ['SJI_fits_to_cube','SJI_to_cube', 'dustbuster']

def SJI_fits_to_cube(filelist, start=0, stop=None, skip=None):
    """
    Read SJI files and return a MapCube. Inherits
    sunpy.instr.iris.SJI_to_cube to stitch multiple sji fits
    files stitched into one mapcube.
    Also sets colormap and normalization
    Assumes hdu index of 0.


    Parameters
    ----------
    filelist: `str` or `list`
        File to read, if single file string detected, will
        create a list of length 1.

    start: `int`
        Temporal axis index to create MapCube from

    stop: `int`
        Temporal index to stop MapCube at

    skip: `int`
        Temporal index to skip over MapCube at


    Returns
    -------
    iris_cube: `sunpy.map.MapCube`
        A map cube of the SJI sequence
    """

    if type(filelist) == str:
        filelist = [filelist]
    iris_cube = sunpy.map.MapCube()
    for fname in filelist[0:]:
        newmap = sunpy.instr.iris.SJI_to_cube(fname)
        for frame in newmap[start:stop:skip]:
            if frame.mean() < frame.max():
                cmap = cm.get_cmap(frame._get_cmap_name())
                vmax = frame.mean()+3.*frame.std()
                frame.plot_settings['cmap'] = cmap
                frame.plot_settings['norm'] = colors.LogNorm(1, vmax)
                #  todo: iris_intscale
                iris_cube.maps.append(frame)

    #  todo: pointing correction(rot_hpc)


    return iris_cube

def SJI_to_cube(filename, start=0, stop=None, hdu=0):
    """
    Read a SJI file and return a MapCube
    .. warning::
        This function is a very early beta and is not stable. Further work is
        on going to improve SunPy IRIS support.
    Parameters
    ----------
    filename: string
        File to read
    start: int
        Temporal axis index to create MapCube from
    stop: int
        Temporal index to stop MapCube at
    hdu: int
        Choose hdu index
    Returns
    -------
    iris_cube: sunpy.map.MapCube
        A map cube of the SJI sequence
    """

    hdus = sunpy.io.read_file(filename)
    # Get the time delta
    time_range = sunpy.time.TimeRange(hdus[hdu][1]['STARTOBS'],
                                      hdus[hdu][1]['ENDOBS'])
    splits = time_range.split(hdus[hdu][0].shape[0])

    if not stop:
        stop = len(splits)

    headers = [hdus[hdu][1]]*(stop-start)
    datas = hdus[hdu][0][start:stop]

    # Make the cube:
    iris_cube = sunpy.map.Map(list(zip(datas, headers)), cube=True)
    # Set the date/time

    for i, m in enumerate(iris_cube):
        m.meta['DATE-OBS'] = splits[i].center.isoformat()

    return iris_cube

def dustbuster(mc):
    """
    Read SJI fits files and return Inpaint-corrected fits files.
    Image inpainting involves filling in part of an image or video
    using information from the surrounding area.

    Parameters
    ----------
    mc: `sunpy.map.MapCube`
        Mapcube to read


    Returns
    -------
    mc: `sunpy.map.MapCube`
        Inpaint-corrected Mapcube
    """
    image_result = []
    ndx = len(mc)
    for i, map in enumerate(mc):
        image_orig = map.data
        nx = map.meta.get('NRASTERP')
        firstpos = range(ndx)[0::nx]
        #  Create mask with values < 1, excluding frame (-200)
        m = ma.masked_inside(image_orig,-199,.1)

        if nx <= 50:  # sparse/coarse raster
            skip = 1
            secpos = [-1]
            thirdpos = [-1]
        elif nx > 50:  # dense raster
            skip = 5
            secpos = range(ndx)[1::nx]
            thirdpos = range(ndx)[2::nx]

        if (i in firstpos) or (i in secpos) or (i in thirdpos):
            image_inpaint = mc[i + skip].data.copy()  # grab next frame
        else:
            image_inpaint = mc[i - skip].data.copy()  # grab prev frame

        # Inpaint mask onto image
        image_orig[m.mask] = image_inpaint[m.mask]

        map.data = image_orig

    return mc