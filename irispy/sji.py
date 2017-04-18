import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.cm as cm
import sunpy.instr.iris
import sunpy.map
import numpy as np
import matplotlib.colors as colors
from astropy.io import fits
'''
This module provides movie tools for level 2 IRIS SJI fits file
'''


def sji_fits_to_cube(filelist, start=0, stop=None, skip=None, grid=False):
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

    #  Option to overlay sun coordinates
    if grid:
        iris_cube[0].draw_grid(linestyle='--', color='orange')

    return iris_cube


def save_to_mp4(ani, outputfile, fps=60):
    '''
    Set up ffmpeg writer
    Save an animation of the MapCube
    Example:
    mc = sunpy.map.Map(files, cube=True)
    ani = mc.plot()
    save_to_mp4(ani, "myfile.mp4")

    Parameters
    ----------
    ani:`matplotlib.animation.FuncAnimation`
        Input animation
    Outputfile: `str`
        Path name of output file (.mp4)
    fps: `int`
        Desired frames per seconds (default=60)

    Return
    -----------
    NoneType
        mp4 to outputfile path
    '''

    plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
    writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='SunPy'),
                                    bitrate=64000)
    if not writer.isAvailable():
        path = input('Set FFMpeg path: (ex. /usr/local/bin/ffmpeg)')
        plt.rcParams['animation.ffmpeg_path'] = path
    writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='SunPy'),
                                    bitrate=64000)
    if not writer.isAvailable():
        print('Cannot find available FFMpeg Writer')
        return

    ani.save(outputfile, writer=writer)
    print('Successfully saved to ' + outputfile)


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
               Optional, save inpaint-corrected fits file
           """
    image_result = []
    ndx = len(mc)
    for i, map in enumerate(mc):
        image_orig = map.data
        nx = map.meta.get('NRASTERP')
        firstpos = range(ndx)[0::nx]
        #  Create mask with values < 1, excluding frame (-200)
        mask = np.zeros(image_orig.shape)
        mask[np.where(image_orig < 1)] = 1
        mask[np.where(image_orig == -200)] = 0
        image_fix = image_orig.copy()
        if nx <= 50:  # sparse/coarse raster
            skip = 1
            secpos = [-1]
            thirdpos = [-1]
        elif nx > 50:  # dense raster
            skip = 3
            secpos = range(ndx)[1::nx]
            thirdpos = range(ndx)[2::nx]

        if (i in firstpos) or (i in secpos) or (i in thirdpos):
            image_inpaint = mc[i + skip].data.copy()  # grab next frame
        else:
            image_inpaint = mc[i - skip].data.copy()  # grab prev frame

        # Inpaint mask onto image
        for layer in range(image_fix.shape[-1]):
            image_fix[np.where(mask)] = image_inpaint[np.where(mask)]

        image_result.append(image_fix)  # add corrected image to new list
        map.data = image_fix

    return mc