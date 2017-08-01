import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import vso_query
from irispy.sji import SJICube,SJIMap, dustbuster
import sunpycube
import sunpy.map
from astropy.coordinates import SkyCoord as SC

import matplotlib.colors as colors
import astropy.units as u
import time
import sunpy.physics.differential_rotation as diff_rot

#pyplot setup
plt.rcParams['animation.ffmpeg_path'] = '/Users/shelbe/anaconda/bin/ffmpeg'
writer = animation.FFMpegWriter(fps=60, metadata=dict(artist='SunPy'), bitrate=100000)
plt.rcParams.update({'font.size': 7})





#local directories
local='/Users/shelbe/Documents/IRIS/data/AR12641/'
iris_dir = '/Volumes/G-DRIVE with Thunderbolt/Users/Shelbe/IRIS/SJI/'
aia_dir = local + 'AIA/'

#grab all .fits in directory
import os
sji_flist=[]
for file in os.listdir(local):
    if file.endswith(".fits"):
        sji_flist.append(os.path.join(local, file))

aia_flist=[]
for file in os.listdir(aia_dir):
    if file.endswith(".fits"):
        aia_flist.append(os.path.join(aia_dir, file))

dimensions = u.Quantity([250,250],u.pixel) #resample size if using full disk aia image
aia = sunpy.map.Map(aia_flist)             #Initialize aia maps
for n,map in enumerate(aia[0:2]):           #Resample to 250 x 250 pix
    aia[n]=map.resample(dimensions)

#Setup plot and output file
fig = plt.figure(dpi=200)
ax0 = fig.add_axes([0.07, 0.05, .92,  .96], zorder=0)
ax0.get_xaxis().set_tick_params(direction='out', width=1)
ax0.get_yaxis().set_tick_params(direction='out', width=1)
ax1 = fig.add_axes([.815, .779, .2, .2], zorder=1)

writer.setup(fig, "/Users/shelbe/Documents/IRIS/writer_test2.mp4", dpi=200)
start = time.perf_counter()
ndx=1
# Jitter is the manual adjustments made to align the image TODO: auto-alignment
jitter = [( 0,  0), (-1,  0), (-1, -1), ( 1,  1), ( 0,  0), (-7, -1), ( 3,  1), ( 4,  0), (-2,  1), ( 1,  0),
          (-5,  0), ( 1,  0), ( 4,  0), ( 1,  0), (-3,  0), ( 4,  0), (-1,  1), (-3,  0), ( 3,  0), (-3,  0),
          ( 1,  0), ( 1,  1), (-2,  0), ( 2,  1), ( 1,  1), (-4,  0), ( 1,  0), ( 1, -2), (-2,  0), ( 3,  1),
          ( 1, -1), ( 2,  0), (-1,  1), ( 2,  0), (-4, -1), ( 2,  0), (-5,  0), ( 8,  1), (-4,  2), ( 0, -1),
          (-2, -1), ( 6,  1), (-2,  0), ( 2, -1), (-6,  1), ( 4, -1), ( 2,  0), (-6,  0), ( 0,  1), ( 0,  1),
          ( 2,  0), ( 2, -1), (-6,  1), ( 1,  0), ( 5,  0), (-6,  1), ( 2,  1), (-3,  1), ( 6,  0), (-6,  1),
          ( 0, -1), ( 6,  0), (-7,  1), ( 9,  0), ( 1,  0), ( 1,  0), ( 0,  0), (-10, 1), ( 8, -1), (-8,  0),
          ( 0,  1), ( 7,  0), (-1, -1), (-6,  2), ( 0,  0), ( 6,  0), (-1,  0), (-6,  3), ( 0,  0), ( 6,  0),
          ( 1,  0), ( 1,  0), ( 1,  0), (-8,  0), ( 5,  0), ( 1,  0), ( 0,  0), ( 1,  0), (-3,  0), ( 1,  -2),
          ( 0,  -2), (-3,  -2), ( 1,  -2), (-7,  -1), ( 8, -5), (-9,  -1), (-14, 1), ( -2, -6), (-11, -4), (-13,  -1),
          (-2, -23), (-2, -23), (-2, -3),(-2, -3), (-2, -3), (0, 0)]
#jitter61 = [(0, 0), (3, 0), (1, 2), (6, 0), (0, 2), (-5, 1), (0, 2), (9, 2), (0, 0), (-4, 0),
#            (2, 0), (0, 0), (1, 2), (0, 0), (0, -1), (0, 2), (0, 0), (0, 0), (-2, 0), (-18, 0),
#            (-6, 0), (-7, 0), (0, 0), (0, 0)]
#[(0, 0), (-4, 0), (0, 0), (-1, 0), (0, 0), (-3, 0), (-2, -8), (-2, -8), (-4, -8),

for i,file in enumerate(sji_flist[0:2]):

    mc = SJICube(file)
    mc = dustbuster(mc)

    print('Getting FITS data from ' + file)


    crval1 = mc[0].meta.get('CRVAL1')
    crval2 = mc[0].meta.get('CRVAL2')
    print(crval1,crval2)
    cdelt1 = mc[0].meta.get('CDELT1')
    cdelt2 = mc[0].meta.get('CDELT2')
    fovx = mc[0].meta.get('FOVX')
    fovy = mc[0].meta.get('FOVY')
    dsun = mc[0].meta.get('DSUN_OBS')

    wave = int(mc[0].meta.get('TWAVE1'))
    title = mc[0].meta.get('TELESCOP') + ' ' + mc[0].meta.get('INSTRUME') + ' $' + str(wave) + r'\AA$'

    nx = int(mc[0].meta.get('NRASTERP')) # number of raster positions


    #AIA query for data, ask Shelbe for vso_query.aiaVSO
    #run = raw_input("download aia images?  y/n: ")
    #if run == 'y':
    #    vso_query.aiaVSO(date_start, time_delta, 5, '1600', aia_dir + '/AIA' + str(ndx) + '/{file}')





    ndx+=1

    if i == 0:
        #Set Starting Coordinates (arcsec)
        xcen = (mc[0].meta['crval1'])* u.arcsec
        ycen = (mc[0].meta['crval2'])* u.arcsec
        xlength = .45*fovx * u.arcsec #Crop map, TODO:Sunpy submap
        ylength = .51*fovy * u.arcsec
        x0 = xcen - .5 * xlength
        y0 = ycen - .5 * ylength
        bl = SC(x0, y0, frame= 'helioprojective')
        tr = SC(x0 + xlength, y0 + ylength, frame='helioprojective')


        #Plot Boundaries

        xmin = float(x0 / u.arcsec)
        xmax = float((x0 + xlength) / u.arcsec)
        ymin = float(y0 / u.arcsec)
        ymax = float((y0 + ylength) / u.arcsec)

        #Rotate Helio Projective Coordinates from start-end time (time elapsed during obs)
        st = datetime.datetime.strptime(mc[0].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        et = datetime.datetime.strptime(mc[-1].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        (xcen, ycen) = diff_rot.rot_hpc(xcen, ycen, st, et, frame_time='synodic', rot_type='snodgrass')



    if xmax < 990 and i>0:
        #Get new start time
        st = datetime.datetime.strptime(mc[0].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        # Rotate Helio Projective Coordinates from end-start time (time elapsed between obs)
        (xcen, ycen) = diff_rot.rot_hpc(xcen, ycen, et, st, frame_time='synodic', rot_type='snodgrass')

        xcen = xcen + jitter[i][0] * u.arcsec     # jitter[j] = (xshift, yshift)
        ycen = ycen + jitter[i][1] * u.arcsec     # shift center to align with prev frame

        #New Overlay and Plot Boundaries
        x0 = xcen - .5 * xlength
        y0 = ycen - .5 * ylength
        bottom_left = (x0, y0)*u.arcsec
        square = [bottom_left, xlength, ylength]

        xmin = float(x0 / u.arcsec)
        xmax = float((x0 + xlength) / u.arcsec)
        ymin = float(y0 / u.arcsec)
        ymax = float((y0 + ylength) / u.arcsec)

        #Get new end time
        et = datetime.datetime.strptime(mc[-1].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        # Rotate Helio Projective Coordinates from start-end time (time elapsed during obs)
        (xcen, ycen) = diff_rot.rot_hpc(xcen, ycen, st, et, frame_time='synodic', rot_type='snodgrass')

        #(Rotating during obs helps minimizes jitter)




    for j, sji in enumerate(mc):
        plot=False

        if (sji.min() < sji.max()): #Filter out bad frames
            #Plot Settings
            aia[i].plot_settings['norm'] = colors.LogNorm(10, 2000)
            print('Mean: ',sji.mean())
            #sji.plot(axes=ax0)
            submap=sji.submap(bl,tr)
            submap.plot_settings['norm'] = colors.LogNorm(1, 400)
            submap.plot(axes=ax0)
            #Label axes
            ax0.set_title(title, visible=False)
            ax0.set_xlabel('X [Arcsec]')
            ax0.set_ylabel('Y [Arcsec]')
            #ax0.set_xlim(xmin, xmax)
            #ax0.set_ylim(ymin, ymax)
            ax0.set_facecolor('black')
            note = ax0.annotate(title + ' ' + sji.meta['DATE-OBS'][:-7], xy=(xmin, ymin),
                                color='white')

            aia[i].plot(axes=ax1)
            aia[i].draw_rectangle(bl, xlength, ylength, color='black')
            ax1.set_title('', visible=False)
            ax1.set_xlabel('', visible=False)
            ax1.set_ylabel('', visible=False)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_facecolor('black')

            #Save plot to output file
            writer.grab_frame()

            #Clear Plots
            ax0.cla()
            ax1.cla()


end = time.perf_counter() - start

print(end) #Print time ran










