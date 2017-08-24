import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import vso_query
from irispy.sji import SJICube, intscale, dustbuster
import sunpy.map
from astropy.coordinates import SkyCoord as SC

import matplotlib.colors as colors
import astropy.units as u
import time
import sunpy.physics.differential_rotation as diff_rot
### sji_animator is a procedural script that combines IRIS SJI and AIA data sets
### Created by Shelbe Timothy, 2016

#pyplot setup
plt.rcParams['animation.ffmpeg_path'] = '/Users/shelbe/anaconda/bin/ffmpeg'
writer = animation.FFMpegWriter(fps=60, metadata=dict(artist='SunPy'), bitrate=100000)
plt.rcParams.update({'font.size': 7})


#local directories
local='/Users/shelbe/Documents/IRIS/data/CH/'
iris_dir = '/Volumes/G-DRIVE with Thunderbolt/Users/Shelbe/IRIS/SJI/'
aia_dir = local + 'AIA/'
outputfile = '/Users/shelbe/Documents/IRIS/writer_test.mp4'

#grab all .fits in directory
import os
sji_list=[]
for file in os.listdir(local):
    if file.endswith(".fits"):
        sji_list.append(os.path.join(local, file))

aia_list=[]
for file in os.listdir(aia_dir):
    if file.endswith(".fits"):
        aia_list.append(os.path.join(aia_dir, file))

dimensions = u.Quantity([250,250],u.pixel) #resample size if using full disk aia image
aia = sunpy.map.Map(aia_list)              #Initialize aia maps
for n,map in enumerate(aia[0:]):           #Resample to 250 x 250 pix
    aia[n]=map.resample(dimensions)

#Setup plot
fig = plt.figure(dpi=200)
ax0 = fig.add_axes([0.07, 0.05, .92,  .96], zorder=0)
ax0.get_xaxis().set_tick_params(direction='out', width=1)
ax0.get_yaxis().set_tick_params(direction='out', width=1)
ax1 = fig.add_axes([.815, .779, .2, .2], zorder=1)

writer.setup(fig, outputfile, dpi=100)
start = time.perf_counter()

#Build plots
for i,file in enumerate(sji_list[0:15]):

    mc = SJICube(file)
    intscale(mc)
    dustbuster(mc)

    print('Getting FITS data from ' + file)

    crval1 = mc[0].meta.get('CRVAL1')
    crval2 = mc[0].meta.get('CRVAL2')
    cdelt1 = mc[0].meta.get('CDELT1')
    cdelt2 = mc[0].meta.get('CDELT2')
    fovx = mc[0].meta.get('FOVX')
    fovy = mc[0].meta.get('FOVY')
    dsun = mc[0].meta.get('DSUN_OBS')

    wave = int(mc[0].meta.get('TWAVE1'))
    title = mc[0].meta.get('TELESCOP') + ' ' + mc[0].meta.get('INSTRUME') + ' $' + str(wave) + r'\AA$'

    nx = int(mc[0].meta.get('NRASTERP')) # number of raster positions

    if i == 0:
        #Set Starting Coordinates (arcsec)

        st = datetime.datetime.strptime(mc[0].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        et = datetime.datetime.strptime(mc[-1].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')
        cen = SC((mc[0].meta['crval1'])* u.arcsec, (mc[0].meta['crval2'])* u.arcsec, frame='helioprojective', obstime=st)

        xlength = .45*fovx    #Cropping dimensions
        ylength = .8*fovy
        x0 = (cen.data._lon.arcsec - .5 * xlength)
        y0 = (cen.data._lat.arcsec - .5 * ylength)
        bl = SC(x0 * u.arcsec, y0 * u.arcsec, frame='helioprojective', obstime=st) #Bottom Left and Top Right coordinates
        tr = SC((x0  + xlength)* u.arcsec, (y0 + ylength)* u.arcsec, frame='helioprojective', obstime=st)


    if tr.data._lon.arcsec < 990 and i>0:
        #Get new start time
        st = datetime.datetime.strptime(mc[0].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')

        # Rotate Helio Projective Coordinates from end-start time (time elapsed between obs)
        cen = diff_rot.solar_rotate_coordinate(cen, st, frame_time='synodic', rot_type='snodgrass')

        #New Submap Boundaries
        x0 = (cen.data._lon.arcsec - .5 * xlength)
        y0 = (cen.data._lat.arcsec - .5 * ylength)
        print(x0,y0)
        bl = SC(x0* u.arcsec, y0* u.arcsec, frame='helioprojective', obstime=st)  # Bottom Left and Top Right coordinates
        tr = SC((x0 + xlength)* u.arcsec, (y0 + ylength)* u.arcsec, frame='helioprojective', obstime=st)

        #Get new end time
        et = datetime.datetime.strptime(mc[-1].meta['DATE-OBS'], '%Y-%m-%d' + ' ' + '%H:%M:%S.%f')

        #(Note: Rotating during obs helps minimizes jitter)

    for j, sji in enumerate(mc):

        if (sji.mean() < 5000) and j%10==5: #Filter out bad frames

            #Create Cropped Submap

            submap=sji.submap(bl,tr)

            #Plot Data
            submap.plot(axes=ax0, norm = colors.PowerNorm(.2, 1, 100))
            aia[i].plot(axes=ax1, norm=colors.LogNorm(10, 2000))
            print(i,cen)
            #Draw GUIs
            note = ax0.annotate(title + ' ' + sji.meta['DATE-OBS'][:-7],
                                xy=(bl.data._lon.arcsec, bl.data._lat.arcsec), color='white')
            aia[i].draw_rectangle(bl, xlength*u.arcsec, ylength*u.arcsec, color='black')

            #Label axes
            ax0.set_title(title, visible=False)
            ax0.set_xlabel('X [Arcsec]')
            ax0.set_ylabel('Y [Arcsec]')
            ax0.set_xlim(x0, x0 + xlength)
            ax0.set_ylim(y0, y0 + ylength)
            ax0.set_facecolor('black')
            ax1.set_title('', visible=False)
            ax1.set_xlabel('', visible=False)
            ax1.set_ylabel('', visible=False)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_facecolor('black')

            #Save Plot to output file
            writer.grab_frame()

            #Clear Plots
            ax0.cla()
            ax1.cla()
    # Rotate Helio Projective Coordinates from start-end time (time elapsed during obs)
    #cen = diff_rot.solar_rotate_coordinate(cen, et, frame_time='synodic', rot_type='snodgrass')

end = time.perf_counter() - start

print(end) #Print time ran










