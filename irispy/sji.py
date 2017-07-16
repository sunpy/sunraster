'''
This module provides movie tools for level 2 IRIS SJI fits file
'''

from copy import deepcopy
from datetime import timedelta

import numpy as np
import numpy.ma as ma
import matplotlib.colors as colors
import matplotlib.animation
from pandas import DataFrame
from astropy.table import Table
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import visualization
from astropy.wcs import WCS
import sunpy.io
import sunpy.time
import sunpy.cm as cm
import sunpy.instr.iris
import sunpy.map
from sunpy.map import GenericMap
from sunpy.map.map_factory import Map
from sunpy.visualization.mapcubeanimator import MapCubeAnimator
from sunpy.visualization import wcsaxes_compat
from sunpy.time import parse_time
from sunpy.lightcurve import LightCurve

from irispy import iris_tools

__all__ = ['SJI_fits_to_cube','SJI_to_cube', 'dustbuster', 'SJICube', 'SJIMap']

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

# the following value is only appropriate for byte scaled images
BAD_PIXEL_VALUE = -200
# the following value is only appropriate for unscaled images
BAD_PIXEL_VALUE_UNSCALED = -32768


class SJIMap(GenericMap):
    """
    A 2D IRIS Slit Jaw Imager Map.

    The Interface Region Imaging Spectrograph (IRIS) small explorer spacecraft
    provides simultaneous spectra and images of the photosphere, chromosphere,
    transition region, and corona with 0.33 to 0.4 arcsec spatial resolution,
    2-second temporal resolution and 1 km/s velocity resolution over a
    field-of- view of up to 175 arcsec by 175 arcsec.  IRIS consists of a 19-cm
    UV telescope that feeds a slit-based dual-bandpass imaging spectrograph.

    Slit-jaw images in four different passbands (C ii 1330, Si iv 1400,
    Mg ii k 2796 and Mg ii wing 2830  A) can be taken simultaneously with
    spectral rasters that sample regions up to 130 arcsec by 175 arcsec at a
    variety of spatial samplings (from 0.33 arcsec and up).
    IRIS is sensitive to emission from plasma at temperatures between
    5000 K and 10 MK.

    IRIS was launched into a Sun-synchronous orbit on 27 June 2013.

    References
    ----------
    * `IRIS Mission Page <http://iris.lmsal.com>`_
    * `IRIS Analysis Guide <https://iris.lmsal.com/itn26/itn26.pdf>`_
    * `IRIS Instrument Paper <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0196&file_type=pdf>`_
    * `IRIS FITS Header keywords <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0077&file_type=pdf>`_
    """

    def __init__(self, data, header, **kwargs):
        GenericMap.__init__(self, data, header, **kwargs)
        if header.get('lvl_num') == 2:
            self.meta['wavelnth'] = header.get('twave1')
            self.meta['detector'] = header.get('instrume')
            self.meta['waveunit'] = "Angstrom"
        if header.get('lvl_num') == 1:
            self.meta['wavelnth'] = int(header.get('img_path').split('_')[1])
            self.meta['waveunit'] = "Angstrom"

        self.meta['detector'] = "SJI"
        self.meta['waveunit'] = "Angstrom"
        palette = cm.get_cmap('irissji' + str(int(self.meta['wavelnth'])))
        palette.set_bad('black')
        self.plot_settings['cmap'] = palette
        self.plot_settings['norm'] = ImageNormalize(stretch=visualization.AsinhStretch(0.1))


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an IRIS SJI image"""
        tele = header.get('TELESCOP', '').startswith('IRIS')
        obs = header.get('INSTRUME', '').startswith('SJI')
        level = header.get('lvl_num') == 1
        return tele and obs

    def draw_slit(self, axes=None, **kwargs):
        """Draws the slit location over the SJI observation.

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
        Axes to plot limb on or None to use current axes.

        Returns
        -------
        lines: list
            A list of `matplotlib.axvline` objects that have been plotted.

        Notes
        -----
        keyword arguments are passed onto matplotlib.pyplot.plot
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        axes.axvline(self.meta['SLTPX1IX'])


class SJICube(object):
    """
    SJICube

    A series of SJI Images

    Parameters
    ----------

    Attributes
    ----------

    Examples
    --------
    >>> from irispy.sji import SJICube
    >>> from irispy.data import sample
    >>> sji = SJICube(sample.SJI_CUBE_1400)   # doctest: +SKIP

    """
    #pylint: disable=W0613,E1101
    def __init__(self, input):
        """Creates a new instance"""
        if isinstance(input, str):
            fits = pyfits.open(input, memmap=True, do_not_scale_image_data=True)
            # TODO find a new masking value for unscaled data
            #self.data = np.ma.masked_less_equal(fits[0].data, 0)
            self.data = fits[0].data
            self.mask = np.ma.masked_equal(fits[0].data, BAD_PIXEL_VALUE_UNSCALED).mask
            reference_header = deepcopy(fits[0].header)
            table_header = deepcopy(fits[1].header)
            # fix reference header
            if reference_header.get('lvl_num') == 2:
                reference_header['wavelnth'] = reference_header.get('twave1')
                reference_header['detector'] = reference_header.get('instrume')
                reference_header['waveunit'] = "Angstrom"
                reference_header['obsrvtry'] = reference_header.get('telescop')
            # check consistency
            if reference_header['NAXIS3'] != self.data.shape[0]:
                raise ValueError("Something is not right with this file!")

            number_of_images = self.data.shape[0]
            metas = []
            dts = fits[1].data[:, fits[1].header['TIME']]
            file_wcs = WCS(fits[0].header)
            # Caution!! This has not been confirmed for non-zero roll
            # angles.
            self.slit_center_sji_indices_x = fits[1].data[:, fits[1].header['SLTPX1IX']]
            self.slit_center_sji_indices_y = fits[1].data[:, fits[1].header['SLTPX2IX']]
            slit_center_positions = file_wcs.celestial.all_pix2world(
                self.slit_center_sji_indices_x, self.slit_center_sji_indices_y,
                iris_tools.WCS_ORIGIN)
            self.slit_center_position_x = u.Quantity(slit_center_positions[0], 'deg').to("arcsec")
            self.slit_center_position_y = u.Quantity(slit_center_positions[1], 'deg').to("arcsec")

            # append info in second hdu to each header
            for i in range(number_of_images):
                metas.append(deepcopy(reference_header))
                metas[i]['DATE_OBS'] = str(parse_time(reference_header['STARTOBS']) + timedelta(seconds=dts[i]))
                # copy over the individual header fields
                for item in fits[1].header[7:]:
                    metas[i][item] = fits[1].data[i, fits[1].header[item]]
                    if item.count('EXPTIMES'):
                        metas[i]['EXPTIME'] = fits[1].data[i, fits[1].header[item]]

            self._meta = metas
        elif len(input) > 1:
            self.data = input[0]
            self._meta = input[1]
            number_of_images = self.data.shape[0]

        norm = []
        cmap = []
        for i in range(number_of_images):
            norm.append(deepcopy(self[i].plot_settings['norm']))
            cmap.append(deepcopy(self[i].plot_settings['cmap']))

        self.plot_settings = {'norm': norm, 'cmap': cmap}
        self.ref_index = 0

    def _get_map(self, index):
        return SJIMap(self.data[index, :, :], self._meta[index], mask=self.mask[index, :, :])

    def __getitem__(self, key):
        """Overriding indexing operation.  If the key results in a single map,
        then a map object is returned.  This allows functions like enumerate to
        work.  Otherwise, a mapcube is returned."""
        if isinstance(key, int):
            return self._get_map(key)
        elif isinstance(key, slice):
            return SJICube((self.data[key, :, :], self._meta[key]))

    def __len__(self):
        """Return the number of maps in a mapcube."""
        return self.data.shape[0]

    def __repr__(self):
        return (
"""SunPy {dtype!s}
---------
Observatory:\t {obs}
Instrument:\t {inst}
Detector:\t {det}
Measurement:\t {meas}
Wavelength:\t {wave}
Obs. Start:\t {date_start:{tmf}}
Obs. End:\t {date_end:{tmf}}
Num. of Frames:\t {frame_num}
IRIS Obs. id:\t {obs_id}
IRIS Obs. Description:\t {obs_desc}
Dimensions:\t {dim}
Scale:\t\t {scale}
""".format(dtype=self.__class__.__name__,
           obs=self.observatory, inst=self.instrument, det=self.detector,
           meas=self.measurement, wave=self.wavelength, date_start=self.date[0],
           date_end=self.date[-1], frame_num=len(self),
           dim=u.Quantity(self.dimensions),
           scale=u.Quantity(self.scale), obs_id=self.iris_obs_id, obs_desc=self.iris_obs_description,
           tmf=TIME_FORMAT) + self.data.__repr__())

    # Sorting methods
    @classmethod
    def _sort_by_date(cls):
        return lambda m: m.date # maps.sort(key=attrgetter('date'))

    def _derotate(self):
        """Derotates the layers in the MapCube"""
        pass

    @property
    def iris_obs_id(self):
        """IRIS Observation ID"""
        return self._meta[self.ref_index].get('OBSID', "")

    @property
    def iris_obs_description(self):
        """IRIS Observation Description"""
        return self._meta[self.ref_index].get('OBS_DESC', "")

    @property
    def instrument(self):
        """Instrument name"""
        return self._meta[self.ref_index].get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image"""
        return u.Quantity(self._meta[self.ref_index].get('wavelnth', 0), self._meta[self.ref_index].get('waveunit', ""))

    @property
    def wavelength(self):
        """wavelength of the observation"""
        return u.Quantity(self._meta[self.ref_index].get('wavelnth', 0), self._meta[self.ref_index].get('waveunit', ""))

    @property
    def observatory(self):
        """Observatory or Telescope name"""
        return self._meta[self.ref_index].get('obsrvtry', self._meta[self.ref_index].get('telescop', "")).replace("_", " ")

    @property
    def detector(self):
        """Detector name"""
        return self._meta[self.ref_index].get('detector', "")

    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return self._meta[0].get('NAXIS1'), self._meta[0].get('NAXIS2')

    @property
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @property
    def date(self):
        """Observation time"""
        return [parse_time(m.get('date_obs')) for m in self._meta]

    @property
    def scale(self):
        """
        Image scale along the x and y axes in units/pixel (i.e. cdelt1,
        cdelt2)
        """
        #TODO: Fix this if only CDi_j matrix is provided
        return self._get_map(self.ref_index).scale

    @property
    def spatial_units(self):
        """
        Image coordinate units along the x and y axes (i.e. cunit1,
        cunit2).
        """
        return self._get_map(self.ref_index).units

    @property
    def exposure_time(self):
        """Exposure time of each frame in seconds."""
        return u.Quantity([m.get('exptime') for m in self._meta], 's')

    @property
    def rotation_matrix(self):
        return self._get_map(self.ref_index).rotation_matrix

    def meta(self, key):
        """The Meta."""
        result = [meta[key] for meta in self._meta]
        # check to see if they are all the same if so just return one value
        if len(set(result)) == 1:
            result = set(result)
        return result

    def submap(self, range_a, range_b, range_c=None):
        """Returns a submap of the map with the specified range.
        """
        new_maps = []
        for i in range(0, len(self)):
            new_maps.append(self._get_map(i).submap(range_a, range_b))
        data = np.zeros([len(self), new_maps[0].data.shape[0], new_maps[0].data.shape[1]])
        data = np.ma.masked_less_equal(data, 0)
        _meta = []
        for i in range(0, len(self)):
            data[i,:,:] = new_maps[i].data
            _meta.append(deepcopy(new_maps[i].meta))
        return SJICube((data, _meta))

    def lightcurve(self, location_a, location_b, range_c=None):
        """Given a pixel index return a lightcurve."""
        return LightCurve(DataFrame({"{0},{1}".format(location_a, location_b):self.data[:, location_a,location_b]},
                                    index=self.date))

    @u.quantity_input(dimensions=u.pixel, offset=u.pixel)
    def superpixel(self, dimensions, offset=(0, 0)*u.pixel, func=np.sum):
        """Returns a new map consisting of superpixels formed from the
        original data.  Useful for increasing signal to noise ratio in images.
        """
        new_maps = []
        print(offset)
        for i in range(0, len(self)):
            new_maps.append(self._get_map(i).superpixel(dimensions, offset=offset, func=func))
        return SJICube(new_maps)

    @u.quantity_input(dimensions=u.pixel)
    def resample(self, dimensions, method='linear'):
        """Returns a new Map that has been resampled up or down"""
        new_maps = []
        for i in range(0, len(self)):
            new_maps.append(self._get_map(i).resample(dimensions, method=method))
        return SJICube(new_maps)

    def apply_function(self, function, *function_args, **function_kwargs):
        """
        Apply a function that operates on the full 3-d data in the mapcube and
        return a single 2-d map based on that function.
        :param function: a function that takes a 3-d numpy array as its first
        argument.
        :param function_args: function arguments
        :param function_kwargs: function keywords
        :return: `sunpy.map.Map`
            A map that stores the result of applying the function to the 3-d
            data of the mapcube.
        """
        if "mapcube_index" in function_kwargs:
            mapcube_index = function_kwargs.pop("mapcube_index")
        else:
            mapcube_index = 0
        return SJIMap(function(self.data, *function_args, **function_kwargs), self._meta[mapcube_index])

    def std(self):
        """
        Calculate the standard deviation of the data array.
        """
        return SJIMap(np.std(self.data, axis=0), self._meta[self.ref_index])

    def mean(self):
        """
        Calculate the mean of the data array.
        """
        return SJIMap(np.mean(self.data, axis=0), self._meta[self.ref_index])

    def min(self):
        """
        Calculate the minimum value of the data array.
        """
        return SJIMap(np.min(self.data, axis=0), self._meta[self.ref_index])

    def max(self):
        """
        Calculate the maximum value of the data array.
        """
        return SJIMap(np.max(self.data, axis=0), self._meta[self.ref_index])

    def plot(self, axes=None, resample=None, annotate=True, interval=200,
             plot_function=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapCube

        Parameters
        ----------
        gamma: float
            Gamma value to use for the color map

        axes: mpl axes
            axes to plot the animation on, if none uses current axes

        resample: list or False
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to
            plot at 1/4 resolution.
            [Note: this will only work where the map arrays are the same size]

        annotate: bool
            Annotate the figure with scale and titles

        interval: int
            Animation interval in ms

        plot_function : function
            A function to be called as each map is plotted. Any variables
            returned from the function will have their ``remove()`` method called
            at the start of the next frame so that they are removed from the plot.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> import matplotlib.animation as animation
        >>> from sunpy.map import Map

        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.plot(colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map at 1/2 original resolution

        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.plot(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Save an animation of the MapCube

        >>> cube = Map(res, cube=True)   # doctest: +SKIP

        >>> ani = cube.plot()   # doctest: +SKIP

        >>> Writer = animation.writers['ffmpeg']   # doctest: +SKIP
        >>> writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)   # doctest: +SKIP

        >>> ani.save('mapcube_animation.mp4', writer=writer)   # doctest: +SKIP

        Save an animation with the limb at each time step

        >>> def myplot(fig, ax, sunpy_map):
        ...    p = sunpy_map.draw_limb()
        ...    return p
        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.peek(plot_function=myplot)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP
        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self._get_map(self.ref_index).wcs)
        fig = axes.get_figure()

        if not plot_function:
            plot_function = lambda fig, ax, smap: []
        removes = []

        # Normal plot
        def annotate_frame(i):
            axes.set_title("{s.name}".format(s=self[i]))

            # x-axis label
            if self[0].coordinate_system.x == 'HG':
                xlabel = 'Longitude [{lon}'.format(lon=self[i].spatial_units.x)
            else:
                xlabel = 'X-position [{xpos}]'.format(xpos=self[i].spatial_units.x)

            # y-axis label
            if self[0].coordinate_system.y == 'HG':
                ylabel = 'Latitude [{lat}]'.format(lat=self[i].spatial_units.y)
            else:
                ylabel = 'Y-position [{ypos}]'.format(ypos=self[i].spatial_units.y)

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        if resample:
            # This assumes that the maps are homogeneous!
            # TODO: Update this!
            resample = np.array(len(self)-1) * np.array(resample)
            ani_data = [self._get_map(j).resample(resample) for j in range(0, len(self))]
        else:
            ani_data = [self._get_map(j) for j in range(0, len(self))]

        im = ani_data[0].plot(axes=axes, **kwargs)
        im.set_cmap(self.plot_settings['cmap'][0])
        im.set_norm(self.plot_settings['norm'][0])

        def updatefig(i, im, annotate, ani_data, removes):
            while removes:
                removes.pop(0).remove()

            im.set_array(ani_data[i].data)
            im.set_cmap(self.plot_settings['cmap'][i])
            norm = deepcopy(self.plot_settings['norm'][i])
            # The following explicit call is for bugged versions of Astropy's ImageNormalize
            # norm.autoscale_None(ani_data[i].data)
            im.set_norm(norm)

            if wcsaxes_compat.is_wcsaxes(axes):
                im.axes.reset_wcs(self[i].wcs)
                wcsaxes_compat.default_wcs_grid(axes)
            else:
                im.set_extent(np.concatenate((self[i].xrange.value,
                                              self[i].yrange.value)))

            if annotate:
                annotate_frame(i)
            removes += list(plot_function(fig, axes, self[i]))

        ani = matplotlib.animation.FuncAnimation(fig, updatefig,
                                                 frames=list(range(0, len(self))),
                                                 fargs=[im, annotate, ani_data, removes],
                                                 interval=interval,
                                                 blit=False)

        return ani

    def plot_in_notebook(self):
        """Provides ipython widgets to plot the SJI data."""
        # TODO: write function
        import ipywidgets as widgets
        pass

    def to_html5_video(self):
        """Output video of the observation"""
        # TODO: write function
        pass

    def peek(self, resample=None, **kwargs):
        """
        A animation plotting routine that animates each element in the
        MapCube

        Parameters
        ----------
        fig: mpl.figure
            Figure to use to create the explorer

        resample: list or False
            Draws the map at a lower resolution to increase the speed of
            animation. Specify a list as a fraction i.e. [0.25, 0.25] to
            plot at 1/4 resolution.
            [Note: this will only work where the map arrays are the same size]

        annotate: bool
            Annotate the figure with scale and titles

        interval: int
            Animation interval in ms

        colorbar: bool
            Plot colorbar

        plot_function : function
            A function to call to overplot extra items on the map plot.
            For more information see `sunpy.visualization.MapCubeAnimator`.

        Returns
        -------
        mapcubeanim : `sunpy.visualization.MapCubeAnimator`
            A mapcube animator instance.

        See Also
        --------
        sunpy.visualization.mapcubeanimator.MapCubeAnimator

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sunpy.map import Map

        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.peek(colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map at 1/2 original resolution

        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.peek(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Plot the map with the limb at each time step

        >>> def myplot(fig, ax, sunpy_map):
        ...    p = sunpy_map.draw_limb()
        ...    return p
        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.peek(plot_function=myplot)   # doctest: +SKIP
        >>> plt.show()   # doctest: +SKIP

        Decide you want an animation:

        >>> cube = Map(files, cube=True)   # doctest: +SKIP
        >>> ani = cube.peek(resample=[0.5, 0.5], colorbar=True)   # doctest: +SKIP
        >>> mplani = ani.get_animation()   # doctest: +SKIP
        """

        if resample:
            if self.all_maps_same_shape():
                resample = np.array(len(self) - 1) * np.array(resample)
                for i in range(0, len(self)):
                    self._get_map(i).resample(resample)
            else:
                raise ValueError('Maps do not all have the same shape.')

        return MapCubeAnimator(self, **kwargs)

    def running_difference(self, offset=1, use_offset_for_meta='ahead'):
        """
        Calculate the running difference of the mapcube and return a new mapcube

        Parameters
        ----------

        offset : int
           Calculate the running difference between map 'i + offset' and image 'i'.
        use_offset_for_meta : {'ahead', 'behind'}
           Which meta header to use in layer 'i' in the returned mapcube, either
           from map 'i + offset' (when set to 'ahead') and image 'i' (when set to
           'behind').

        Returns
        -------
        sunpy.map.MapCube
           A mapcube containing the running difference of the input mapcube.

        """
        if offset < 1:
            raise ValueError('The value of the offset keyword must be greater than or equal to 1.')

        # Create a list containing the data for the new map object
        new_mc = []
        for i in range(0, len(self) - offset):
            new_data = self.data[:, :, i + offset] - self.data[:, :, i]
            if use_offset_for_meta == 'ahead':
                new_meta = self._meta[i + offset]
            elif use_offset_for_meta == 'behind':
                new_meta = self._meta[i]
            else:
                raise ValueError('The value of the keyword "use_offset_for_meta" has not been recognized.')
            new_mc.append(Map(new_data, new_meta))

        # Create the new mapcube and return
        return SJICube(new_mc)

    def base_difference(self, base=0, fraction=False):
        """
        Calculate the base difference of a mapcube.

        Parameters
        ----------
        base : int, sunpy.map.Map
           If base is an integer, this is understood as an index to the input
           mapcube.  Differences are calculated relative to the map at index
           'base'.  If base is a sunpy map, then differences are calculated
           relative to that map

        fraction : boolean
            If False, then absolute changes relative to the base map are
            returned.  If True, then fractional changes relative to the base map
            are returned

        Returns
        -------
        sunpy.map.MapCube
           A mapcube containing base difference of the input mapcube.

        """

        if not(isinstance(base, GenericMap)):
            base_data = self.data[:, :, base]
        else:
            base_data = base.data

        if base_data.shape != self.data.shape:
            raise ValueError('Base map does not have the same shape as the maps in the input mapcube.')

        # Fractional changes or absolute changes
        if fraction:
            relative = base_data
        else:
            relative = 1.0

        # Create a list containing the data for the new map object
        new_mc = []
        for i in range(0, len(self)):
            new_data = (self.data[:, :, i] - base_data) / relative
            new_mc.append(Map(new_data, self._meta[i]))

        # Create the new mapcube and return
        return SJICube(new_mc)


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
    for fname in filelist:
        newmap = SJI_to_cube(fname)
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
