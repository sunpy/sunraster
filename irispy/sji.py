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
    filelist: string or list
        File to read, if single file string detected, will
        create a list of length 1.

    start: step
        Temporal axis index to create MapCube from

    stop: step
        Temporal index to stop MapCube at

    skip: step
        Temporal index to skip over MapCube at


    Returns
    -------
    iris_cube: sunpy.map.MapCube
        A map cube of the SJI sequence
    """

    if type(filelist) == str:
        filelist = [filelist]
    iris_cube = sunpy.map.MapCube()
    counter = 1
    setup = sunpy.instr.iris.SJI_to_cube(filelist[0])
    xcen = setup[0].meta.get('xcen')
    ycen = setup[0].meta.get('ycen')
    time = datetime.datetime.strptime(setup[0].meta.get('startobs'),
                                      '%Y-%m-%d' + 'T' + '%H:%M:%S.%f')
    print(xcen, ycen, time)
    for fname in filelist[0:]:
        newmap = sunpy.instr.iris.SJI_to_cube(fname)
        for frame in newmap[start:stop:skip]:
            if frame.mean() < frame.max():
                cmap = cm.get_cmap(frame._get_cmap_name())
                vmax = frame.mean()+3.*frame.std()
                frame.plot_settings['cmap'] = cmap
                frame.plot_settings['norm'] = colors.LogNorm(1, vmax)
                #  todo: iris_intscale
                resultmapcube.maps.append(frame)

        print(str(counter)+' of '+str(len(filelist)))
        counter += 1
    #  todo: pointing correction(rot_hpc)

    #  Option to overlay sun coordinates
    if grid:
        iris_cube[0].draw_grid(linestyle='--', color='orange')

    return iris_cube


def save_to_mp4(mc, outputfile, fps=60):
    '''
    Set up ffmpeg writer
    Configures plot axes
    Save an animation of the MapCube

    Parameters
    ----------
    mc: sunpy.map.MapCube
        Input mapcube
    Outputfile: str
        Path name of output file (.mp4)
    fps: int
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
    plt.rcParams['axes.facecolor'] = 'black'
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    ax = plt.axes(zorder=0)

    ani = mc.plot(ax)
    ani.save(outputfile, writer=writer)
    print('Successfully saved to ' + outputfile)


def dustbuster(filelist):
    """
       Read SJI fits files and return Inpaint-corrected fits files.
       **Warning** Will clobber input file.

       Parameters
       ----------
       filelist: string or list
           File to read, if single file string detected, will
           create a list of length 1.


       Returns
       -------
       NoneType
           Inpaint-corrected fits to clobber input file
       """

    if type(filelist) == str:  # convert single file to list
        file_list = [filelist]
    nfile = 0
    for file in filelist:
        nfile += 1
        image_header = fits.getheader(file)
        image_orig = fits.getdata(file)
        nx = image_header.get('NRASTERP')
        print(image_header)
        # TODO mapcube input

        img_shape = image_orig.shape
        ndx = img_shape[0]

        image_result = []
        counter = 0
        firstpos = range(ndx)[0::nx]

        for i in range(0, ndx):
            #  Create mask with values < 1, excluding frame (-200)
            mask = np.zeros(image_orig[i].shape)
            mask[np.where(image_orig[i] < 1)] = 1
            mask[np.where(image_orig[i] == -200)] = 0

            image_fix = image_orig[i].copy()
            if nx <= 50:     # sparse/coarse raster
                skip = 1
            elif nx > 50:    # dense raster
                skip = 3
                secpos = range(ndx)[1::nx]
                thirdpos = range(ndx)[2::nx]

            if (i in firstpos) or (i in secpos) or (i in thirdpos):
                image_inpaint = image_orig[i + skip].copy()  # grab next frame
            else:
                image_inpaint = image_orig[i - skip].copy()  # grab prev frame

            # Inpaint mask onto image
            for layer in range(image_fix.shape[-1]):
                image_fix[np.where(mask)] = image_inpaint[np.where(mask)]

            image_result.append(image_fix)  # add corrected image to new list
            counter += 1
            print(str(counter) + ' of ' + str(ndx))

        #  overwrite old file
        print(len(image_result))
        outputfile = file
        fits.writeto(outputfile, image_result, header=image_header,
                     output_verify='fix', clobber=True)
        print("Saved to " + outputfile)
