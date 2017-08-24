import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.style as style
import sunpy.cm as cm
import sunpy.map
import numpy as np
import matplotlib.colors as colors
import astropy.units as u

from astropy.io import fits
plt.rcParams['animation.ffmpeg_path'] = '/Users/shelbe/anaconda/bin/ffmpeg'
plt.rcParams.update({'font.size': 4})
#style.use( 'dark_background')
writer = animation.FFMpegWriter(fps=45, metadata=dict(artist='SunPy'), bitrate=100000)


iris_level2 = '/links/kale/iris/data/level2/2017/'
iris_dir = '/Volumes/HDrive/Users/Shelbe/IRIS/Raster'
import os
raster_list=[]
for file in os.listdir(iris_dir):
    if file.endswith(".fits"):
        raster_list.append(os.path.join(iris_dir, file))

print(len(raster_list))

def raster(file_list):
    fig, ax = plt.subplots(2, 3, sharex=True, sharey=True)
    plt.subplots_adjust(bottom=0.07, left=0.07, top=.98, right=.98, wspace=0, hspace=0)
    writer.setup(fig, '/Users/shelbe/Documents/IRIS/raster_test.mp4', dpi=120)
    fig.set_facecolor('black')

    j = 0

    for n,i in enumerate(file_list[0:]):
        print( 'Getting FITS data from ' + i)

        main_header = fits.getheader(i, ext=0)
        exptime = main_header.get('exptime')

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



        if Mg_header == None:
            print('Header not found')
            return

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


        #Transpose for Plot x-axis label
        Mg_header['naxis1'] = naxis3
        Mg_header['crpix1'] = crpix3
        Mg_header['crval1'] = crval3
        Mg_header['cdelt1'] = cdelt3
        Mg_header['ctype1'] = ctype3
        Mg_header['cunit1'] = cunit3
        Mg_header['crval2'] = crval2

        date_start = datetime.datetime.strptime(main_header.get('STARTOBS'), '%Y-%m-%d' + 'T' + '%H:%M:%S.%f')
        date_end = datetime.datetime.strptime(main_header.get('ENDOBS'), '%Y-%m-%d' + 'T' + '%H:%M:%S.%f')
        dt = main_header.get('RASNRPT')
        nt = main_header.get('RASRPT')
        time_delta = (date_end - date_start) / dt
        timestamp = date_start+(time_delta*(nt-.5))

        framesize = 40
        xcen = (crval3) * u.arcsec
        ycen = (crval2 ) * u.arcsec
        xmin = (float(xcen / u.arcsec - framesize))
        xmax = (float(xcen / u.arcsec + framesize))
        ymin = (float(ycen / u.arcsec - framesize))
        ymax = (float(ycen / u.arcsec + framesize))


        wavmin = crval1 - crpix1 * cdelt1
        title = []
        wavelength = round(2832.0, 1)

        title.append('$' + str(wavelength) + '\AA$')

        wave_idx2 = []
        wavelength=np.arange(0,naxis1)*(cdelt1) + wavmin
        wave_idx = [2795.75, 2798.75, 2796.15, 2796.35, 2796.55]

        #Find closest wavelength slice to wave_idx
        for i in wave_idx:
            wave=np.where(wavelength >= i)[0]
            print(wavelength[wave[0:5]])
            wave_idx2.append(wave[0])

            title.append(str('$'+ str(round(wavelength[wave[0]],1)) + '\AA$'))

        #Transpose image data to plot solar x v solar y
        mg_wing_image = (Mg_data.T[wave_idx2[0]]/exptime)
        mg_triplet_image = (Mg_data.T[wave_idx2[1]]/exptime)
        mg_k2v_image = (Mg_data.T[wave_idx2[2]]/exptime)
        mg_k3_image = (Mg_data.T[wave_idx2[3]]/exptime)
        mg_k2r_image = (Mg_data.T[wave_idx2[4]]/exptime)
        wing_image = wing_data.T[15]/exptime

        #Initialize Sunpy maps
        mg_wing =sunpy.map.Map(mg_wing_image, Mg_header)
        mg_triplet = sunpy.map.Map(mg_triplet_image, Mg_header)
        mg_k2v = sunpy.map.Map(mg_k2v_image, Mg_header)
        mg_k3 = sunpy.map.Map(mg_k3_image, Mg_header)
        mg_k2r = sunpy.map.Map(mg_k2r_image, Mg_header)
        wing = sunpy.map.Map(wing_image,wing_header)

        if (j <= 80 and mg_k3.data.mean()>65)\
                or(j > 80 and mg_k3.data.mean()>25)\
                or(j>95 and mg_k3.data.mean()>10):  #Filter bad frames
            #Plot Image Data
            mg_wing.plot(axes=ax[0, 1],norm=colors.PowerNorm(.3, 1, 80),cmap=cm.get_cmap(name='sohoeit171'))
            mg_triplet.plot(axes=ax[0, 2], norm=colors.PowerNorm(.5, 1, 40),cmap=cm.get_cmap(name='irissji1330'))
            mg_k2v.plot(axes=ax[1, 0], norm=colors.PowerNorm(.9, 5, 500),cmap=cm.get_cmap(name='irissji1400'))
            mg_k3.plot(axes=ax[1, 1], norm=colors.PowerNorm(.9, 5, 500), cmap=cm.get_cmap(name='irissji1600'))
            mg_k2r.plot(axes=ax[1, 2], norm=colors.PowerNorm(.9, 5, 500), cmap=cm.get_cmap(name='irissji2796'))
            wing.plot(axes=ax[0, 0], norm=colors.PowerNorm(.9, 1, 220),cmap=cm.get_cmap(name='irissji2832'))
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

            ax[1, 0].annotate('IRIS' + ' ' + timestamp.strftime('%Y/%m/%d %H:%M:%S') ,
                              xy=(xmin + 1, ymin + 1), color='white', fontsize=5, zorder=1)
            ax[0, 0].annotate(title[0], xy=(xmin + 1, ymin + 2), color='black', fontsize=6, zorder=1)
            ax[0, 1].annotate(title[1], xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[0, 2].annotate('Mg II Triplet ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 0].annotate('Mg II k$_{2v}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 1].annotate('Mg II k$_{3}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            ax[1, 2].annotate('Mg II k$_{2r}$ ', xy=(xmin + 1, ymin + 5), color='white', fontsize=6, zorder=1)
            #Save to outputfile
            writer.grab_frame()
            #Clear Plots
            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[0, 2].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()
            ax[1, 2].cla()

        if nt == dt:
            j += 1






test = raster(raster_list)