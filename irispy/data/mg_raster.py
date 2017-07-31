import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.cm as cm
import irispy
import sunpy.map
import numpy as np
import matplotlib.colors as colors
import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import visualization
import sunpy.physics.differential_rotation as diff_rot

from astropy.io import fits
#pyplot setup
plt.rcParams['animation.ffmpeg_path'] = '/Users/shelbe/anaconda/bin/ffmpeg'
plt.rcParams.update({'font.size': 4,'axes.facecolor':'black'})
writer = animation.FFMpegWriter(fps=10, metadata=dict(artist='SunPy'), bitrate=100000)

#remote/local directories
iris_level2 = '/links/kale/iris/data/level2/2017/'
iris_dir = '/Volumes/HDrive/Users/Shelbe/IRIS/Raster'

#grab all .fits in directory
import os
raster_list=[]
for file in os.listdir(iris_dir):
    if file.endswith(".fits"):
        raster_list.append(os.path.join(iris_dir, file))

print(len(raster_list))

def raster(file_list):

    # Jitter is the manual adjustments made to align the image TODO: auto-alignment
    jitter = [(  20, 0), ( 16, -1), ( 19, 0), ( 15, -1), (17, -1), ( 15,  9), ( 12,  8), ( 15,  8), ( 12,  8), ( 16,  8),
              (-1, -3), (-4, -3), (0, -3), (-3, -3), (0, -1), (-2, -1), ( 0, -2), (-9, -6), (-14, -5), (-11, -6),
              (-15, -4), (-13, -5), (-15, -6), (-14, -10), (-18, -11), (-17, -11), (-18, -9), (-16, -10), (-19, -9), (-16, -9),
              (-16, -11), (-17, -7), (-17, -9), (-16, -9), (-16, -11), (-18, -7), (2, 1), (0, 2), (-1, 0), (-2, 3),
              (-2, 3), (-2, 3), (-2, 2), (-2, -1), (-2, -1), (-2, -1), (3, -6), (1, 1), (-1, -3), (0, -2),
              (0, -3), (-1, -2), (1, -2), (0, -4), (1, 0), (1, -6), (1, -6), (0, -7), (2, -4), (-1, -8),
              (1, -6), (-1, -6), (1, -1), (0, 1), (0, -1), (0, 1), (2, -1), (-2, -2), (1, 1), (1, 1),
              (-1, 1), (1, 1), (-1, 0), (-1, 0), (8, -13), (3, -14), (7, -14), (7, -14), (4, -14), (6, -14),
              (6, -17), (8, -18), (6, -17), (9, -17), (6, -17), (3, -17), (4, -14), (2, -14), (4, -14), (2, -13),
              (4, -7), (1, -7), (3, -2), (4, -2), (4, -2), (2, -2), (2, -2), (2, -2), (3, -2), (3, -2),
              (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]
    #ar12661 jitter [(-5, 5), (8, 0), (0, 2), (0, -1), (16,-1), (-12,0), (13,0), (-12, 0), (0, -3), (5, 0),
    #          (12, 0), (-10, 0), (12, 0), (-8, -1), (3, 0), (0, -1), (2,3), (0, 0), (5, 2), (0, 0),
    #          (2, 0), (0, 0), (0, 0), (0, 5), (0, 0), (0, 0), (0, 0), (0, 0)]

    #Colormaps downloaded from sunpy.cm
    cmap = cm.get_cmap(name='irissji2832')
    cmap2 = cm.get_cmap(name='irissji1330')
    cmap3 = cm.get_cmap(name='irissji1600')
    cmap4 = cm.get_cmap(name='irissji1400')
    cmap5 = cm.get_cmap(name='irissji2796')
    cmap6 = cm.get_cmap(name='sohoeit171')

    #Setup plot and output file
    fig, ax = plt.subplots(2, 3, sharex=True, sharey=True)
    plt.subplots_adjust(bottom=0.07, left=0.07, top=.98, right=.98, wspace=0, hspace=0)
    writer.setup(fig, '/Users/shelbe/Documents/IRIS/raster_test.mp4', dpi=200)
    frame = 35

    j = 0

    for n,i in enumerate(file_list[0:]):
        #Select file, get header info
        print( 'Getting FITS data from ' + i)
        main_header = fits.getheader(i, ext=0)
        exptime = main_header.get('exptime')

        #Find 2832A and Mg extensions by searching "TDESC" keywords
        Mg_header = None
        for h in range(1,9):
            tdesc='TDESC' + str(h)
            if main_header.get(tdesc) != None:
                if '2832' in main_header.get(tdesc):
                    wing_data = fits.getdata(i, ext=h)
                    wing_header = fits.getheader(i, ext=h)
                    wing_header['naxis1']=wing_header['naxis3']
                    wing_header['crpix1'] = wing_header['crpix3']
                    wing_header['crval1'] = wing_header['crval3']
                    wing_header['cdelt1'] = wing_header['cdelt3']
                    wing_header['ctype1'] = wing_header['ctype3']
                    wing_header['cunit1'] = wing_header['cunit3']

                if 'Mg' in main_header.get(tdesc):
                    Mg_data = fits.getdata(i, ext=h)
                    Mg_header = fits.getheader(i, ext=h)
                    print(tdesc, main_header.get(tdesc),Mg_data[0].shape)


        #Check if Mg header load failed
        if Mg_header == None:
            print('Header not found')
            return

        #Initialize header data
        naxis1 = Mg_header.get('NAXIS1')
        naxis2 = Mg_header.get('NAXIS2')
        naxis3 = Mg_header.get('NAXIS3')
        crpix1 = Mg_header.get('CRPIX1')
        crpix2 = Mg_header.get('CRPIX2')
        crpix3 = Mg_header.get('CRPIX3')
        crval1 = Mg_header.get('CRVAL1')
        crval2 = Mg_header.get('CRVAL2')
        crval3 = Mg_header.get('CRVAL3')
        cdelt1 = Mg_header.get('CDELT1')
        cdelt2 = Mg_header.get('CDELT2')
        cdelt3 = Mg_header.get('CDELT3')
        ctype3 = Mg_header.get('CTYPE3')
        cunit3 = Mg_header.get('CUNIT3')
        date_start = datetime.datetime.strptime(main_header.get('STARTOBS'), '%Y-%m-%d' + 'T' + '%H:%M:%S.%f')
        date_end = datetime.datetime.strptime(main_header.get('ENDOBS'), '%Y-%m-%d' + 'T' + '%H:%M:%S.%f')
        dt = main_header.get('RASNRPT') #Number of total raster
        nt = main_header.get('RASRPT')  #Order of raster
        time_delta = (date_end - date_start) / dt
        timestamp = date_start+(time_delta*(nt-.5))
        print(str(nt)+' of '+str(dt))

        # Transpose z-x axes (should be fixed in IRISSG)
        Mg_header['naxis1'] = naxis3
        Mg_header['crpix1'] = crpix3
        Mg_header['crval1'] = crval3
        Mg_header['cdelt1'] = cdelt3
        Mg_header['ctype1'] = ctype3
        Mg_header['cunit1'] = cunit3
        Mg_header['crval2'] = crval2

        #x and y axis setup
        xcen = (crval3 + jitter[j][0]) * u.arcsec   # jitter[j] = (xshift, yshift)
        ycen = (crval2 + jitter[j][1]) * u.arcsec   # shift center to align with prev frame
        xmin = (float(xcen / u.arcsec - frame))
        xmax = (float(xcen / u.arcsec + frame))
        ymin = (float(ycen / u.arcsec - frame))
        ymax = (float(ycen / u.arcsec + frame))

        #Select wavelengths (angstroms)
        wavmin = crval1 - crpix1 * cdelt1
        title = []
        wavelength0 = 2832.0
        title.append('$' + str(wavelength0) + '\AA$')

        wave_idx2 = []
        wavelength=np.arange(0,naxis1)*(cdelt1) + wavmin
        wave_idx = [2795.75, 2798.75, 2796.15, 2796.35, 2796.55]
        #Find closest slices to the wave_idx list
        for i in wave_idx:
            wave=np.where(wavelength >= i)[0] #search results >=
            wave_idx2.append(wave[0])         #save first result
            title.append(str('$'+ str(round(wavelength[wave[0]],1)) + '\AA$'))

        #Transpose z-x axes (should be fixed in IRISSG)
        mg_wing_image = Mg_data.T[wave_idx2[0]]/exptime
        mg_triplet_image = Mg_data.T[wave_idx2[1]]/exptime
        mg_k2v_image = Mg_data.T[wave_idx2[2]]/exptime
        mg_k3_image = Mg_data.T[wave_idx2[3]]/exptime
        mg_k2r_image = Mg_data.T[wave_idx2[4]]/exptime
        wing_image = wing_data.T[15]/exptime


        #Initialize Maps, TODO: Move to IRISSG map object
        mg_wing =sunpy.map.Map(mg_wing_image, Mg_header)
        mg_triplet = sunpy.map.Map(mg_triplet_image, Mg_header)
        mg_k2v = sunpy.map.Map(mg_k2v_image, Mg_header)
        mg_k3 = sunpy.map.Map(mg_k3_image, Mg_header)
        mg_k2r = sunpy.map.Map(mg_k2r_image, Mg_header)
        wing = sunpy.map.Map(wing_image,wing_header)

        #print('Mean: '+ str(mg_k3.data.mean()))

        if (j <= 78 and mg_k3.data.mean()>80) or((j > 78 and mg_k3.data.mean()>10)): #Filter out bad frames
            #Plot Settings
            mg_wing.plot_settings['cmap'] = cmap6
            mg_wing.plot_settings['norm'] = colors.PowerNorm(.3, 0, 80)
            mg_wing.plot(axes=ax[0, 1])
            mg_triplet.plot_settings['cmap'] = cmap2
            mg_triplet.plot_settings['norm'] = colors.PowerNorm(.6, 0, 40)
            mg_triplet.plot(axes=ax[0, 2])
            mg_k2v.plot_settings['cmap'] = cmap4
            mg_k2v.plot_settings['norm'] = colors.PowerNorm(.9, 0, 500)
            mg_k2v.plot(axes=ax[1, 0])
            mg_k3.plot_settings['cmap'] = cmap3
            mg_k3.plot_settings['norm'] = colors.PowerNorm(.9, 0, 500)
            mg_k3.plot(axes=ax[1, 1])
            mg_k2r.plot_settings['cmap'] = cmap5
            mg_k2r.plot_settings['norm'] = colors.PowerNorm(.9, 0, 500)
            mg_k2r.plot(axes=ax[1, 2])
            wing.plot_settings['cmap'] = cmap
            wing.plot_settings['norm'] = colors.PowerNorm(.9, 0, 220)
            wing.plot(axes=ax[0, 0])

            #Label Axes
            ax[0, 0].set_title(title[0], visible=False)
            ax[0, 1].set_title(title[1], visible=False)
            ax[0, 2].set_title(title[2], visible=False)
            ax[1, 0].set_title(title[3], visible=False)
            ax[1, 1].set_title(title[4], visible=False)
            ax[1, 2].set_title(title[5], visible=False)
            ax[0, 0].set_xlabel('X [Arcsec]', visible=False)
            ax[0, 0].set_ylabel('Y [Arcsec]')
            ax[0, 0].set_xlim(xmin, xmax)
            ax[0, 0].set_ylim(ymin, ymax)
            plt.setp(ax[0, 1].get_yticklabels(), visible=False)
            plt.setp(ax[0, 2].get_yticklabels(), visible=False)
            plt.setp(ax[1, 1].get_yticklabels(), visible=False)
            plt.setp(ax[1, 2].get_yticklabels(), visible=False)

            ax[0, 1].set_xlabel('X [Arcsec]', visible=False)
            ax[0, 1].set_ylabel('Y [Arcsec]', visible=False)

            ax[0, 2].set_xlabel('X [Arcsec]')
            ax[0, 2].set_ylabel('Y [Arcsec]', visible=False)

            ax[1, 0].set_xlabel('X [Arcsec]', visible=False)
            ax[1, 0].set_ylabel('Y [Arcsec]')

            ax[1, 1].set_xlabel('X [Arcsec]')
            ax[1, 1].set_ylabel('Y [Arcsec]', visible=False)

            ax[1, 2].set_xlabel('X [Arcsec]', visible=False)
            ax[1, 2].set_ylabel('Y [Arcsec]', visible=False)

            ax[1, 0].annotate('IRIS' + ' ' + timestamp.strftime('%Y/%m/%d %H:%M:%S'),
                              xy=(xmin + 1, ymin + 1), color='white', fontsize=5, zorder=1)
            ax[0, 0].annotate(title[0], xy=(xmin + 1, ymin + 2), color='black', fontsize=6, zorder=1)
            ax[0, 1].annotate(title[1], xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[0, 2].annotate('Mg II Triplet ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 0].annotate('Mg II k$_{2v}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 1].annotate('Mg II k$_{3}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 2].annotate('Mg II k$_{2r}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)

            #Save frame to output file
            writer.grab_frame()

            #Clear Plots
            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[0, 2].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()
            ax[1, 2].cla()

        #New obs
        if nt == dt:
            j += 1






test = raster(raster_list)