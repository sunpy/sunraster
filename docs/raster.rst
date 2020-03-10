.. _data_classes:

============
Data Classes
============

The sunraster package provides a number of data classes for manipulating and
visualization spectrograph data. These include `~sunraster.Raster` and 
`~sunraster.RasterSequence`.

.. _raster:

Raster
------

The fundamental data class of the sunraster package is `~sunraster.Raster`. 
It is designed to handle data from a single raster scan or sit-and-stare. 
`~sunraster.Raster` stores its data as a single 3-D array whose 
transformations between pixel and real world coordinates are described by 
a single astropy WCS (World Coordinate System) object. 
For data that cannot be described by a single set of WCS transformations, 
see the :ref:`rastersequence` section below. 
`~sunraster.Raster`'s expected axis types are time/slit step, position along slit, 
and spectral. 
No specific ordering is required. 

`~sunraster.Raster` is subclassed from `ndcube.NDCube` and so inherits the 
same attributes for data, wcs, extra_coords, uncertainty, mask, meta, and unit. It also 
inherits much of the same slicing, WCS, and visualization API and functionality.
In addition, `~sunraster.Raster` provides some additional convenience 
properties relevant to its specific data type. Below we will explain how to 
initialize a `~sunraster.Raster` and utilize its functionalities.

Initialization
^^^^^^^^^^^^^^
To initialize the a basic `~sunraster.Raster` object, all you need is a
3-D `numpy.ndarray` containing the data, and an `astropy.wcs.WCS` object
describing the transformation from array-element space to real world
coordinates.  Let's create an array of data representing a single raster scan 
across as region of the Sun.
Let the shape of the array be (3, 4, 5) and let every value be 1.::

  >>> import numpy as np
  >>> data = np.ones((3, 4, 5))

Let's say the first axis represents the slit steps, the second represent the 
pixels along the slit, and the third represent the spectral axis.

Now let's create an `astropy.wcs.WCS` object describing the
translation from the array element coordinates (which we refer to as pixel 
coordinates) to real world coordinates.
Let the first (slit step) axis have real world coordinates 
of helioprojective longitude and latitude.
Note that although we sometimes think of the x-dimension as longitude and the
y-dimension as latitude, latitude and longitude are in fact coupled dimensions.
This means that except in a small number of edge cases, moving along the slit
in y-direction will cause both the latitude and longitude to change, even if
only slightly.
This is important to understand when interacting with the WCS object,
and hence the `~sunraster.Raster` class.
The second axis is the position along the slit and therefore must also 
have real world coordinates of helioprojective longitude and latitude.
Finally let the third (spectral) axis a real world coordinate of wavelength.
Although a WCS object can often be easily created by feeding a FITS header into
the `astropy.wcs.WCS` initializer, here we will create one manually to be 
explicit.
We can do this with the following code.
Note that due to (confusing) convention, the order of the axes in the
WCS object is reversed relative to the data array.::

  >>> import astropy.wcs
  >>> wcs_input_dict = {
  ... 'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 5,
  ... 'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 4,
  ... 'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 3}
  >>> input_wcs = astropy.wcs.WCS(wcs_input_dict)

Now that we have a data array and a corresponding WCS object, we can
create an `~sunraster.Raster` instance simply by doing::

  >>> from sunraster import Raster
  >>> my_raster = Raster(data, input_wcs)

The data array is stored in the ``my_raster.data`` attribute while the
WCS object is stored in the ``my_raster.wcs`` attribute.  However, when
manipulating/slicing the data is it better to slice the object as a
whole.  (See section on :ref:`raster_slicing`.)  So the ``.data`` attribute
should only be used to extract/access specific values in the data.

Thanks to the fact that `~sunraster.Raster` is subclassed from
`~ndcube.NDCube`, you can also supply additional data to the instance. 
These include: metadata (`dict` or dict-like) located at `Raster.meta`;
a data mask (boolean `numpy.ndarray`) located at `Raster.mask` marking, for
example, reliable and unreliable pixels; 
an uncertainty array (`numpy.ndarray`) located at `Raster.uncertainty` describing the
uncertainty of each data array value; 
and a unit (`astropy.units.Unit` or unit `str`). 
For example::

  >>> mask = np.zeros_like(my_raster.data, dtype=bool)
  >>> meta = {"Description": "This is example Raster metadata."}
  >>> my_raster = Raster(data, input_wcs, uncertainty=np.sqrt(data),
  ...                    mask=mask, meta=meta, unit=None)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

N.B. The above warning is due to the fact that
`astropy.nddata.uncertainty` is recommended to have an
``uncertainty_type`` attribute giving a string describing the type of
uncertainty.  However, this is not required.

Dimensions
^^^^^^^^^^

`~sunraster.Raster` has useful properties for inspecting its data shape and
axis types, `~sunraster.Raster.dimensions` and
`~sunraster.Raster.world_axis_physical_types`::

  >>> my_raster.dimensions
  <Quantity [3., 4., 5.] pix>
  >>> my_raster.world_axis_physical_types
  ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'em.wl')

`~sunraster.Raster.dimensions` returns an `~astropy.units.Quantity` of
pixel units giving the length of each dimension in the
`~sunraster.Raster` while `~sunraster.Raster.world_axis_physical_types`
returns an iterable of strings denoting the type of physical property
represented by each axis.  The axis names are in accordance with the
International Virtual Observatory Alliance (IVOA)
`UCD1+ controlled vocabulary <http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`_.
Here the shape and axis types are given in data order, not WCS order.

Slicing
^^^^^^^

Users can slice a `~sunraster.Raster` instance using the standard slicing notation.
For example::

  >>> my_raster_roi = my_raster[1:3, 0:3, 2:4]

This enables users to access sub-regions of their data while simultaneously slicing
not only the other array attributes (e.g. uncertainty, mask, etc.) but
also the WCS object.  This ensures that even though the data array has
changed size and shape, each array element will still correspond to
the same real world coordinates as they did before.
Therefore, slicing in this way is preferred to accessing the `Raster.data` array 
directly unless the user wants to directly access the data values without 
supporting metadata such as WCS information, uncertainties etc.


Slicing can reduce the dimension of an `~sunraster.Raster`, e.g.::

  >>> my_2d_raster = my_raster[1:3, 0:3, 2:4]

In addition to slicing by index, `~sunraster.Raster` supports a basic
version of slicing/indexing by real world coordinates via the
`~sunraster.Raster.crop_by_coords` method.  This takes a list of
`astropy.units.Quantity` instances representing the minimum real world
coordinates of the region of interest in each dimension.  The
order of the coordinates must be the same as the order of the data
axes. 
A second iterable of `~astropy.units.Quantity` must also be provided 
which gives the maximum real world coordinates of the region of interest 
in each data axis::

  >>> import astropy.units as u
  >>> my_raster_roi = my_raster.crop_by_coords([0.7*u.deg, 1.3e-5*u.deg, 1.04e-9*u.m],
  ...                                          [1.3*u.deg, 1.000013*u.deg, 1.12e-9*u.m])

This method does not rebin or interpolate the data if the region of interest
does not perfectly map onto the array's "pixel" grid.  Instead
it translates from real world to pixel coordinates and rounds to the
nearest integer pixel before indexing/slicing the `~sunraster.Raster`
instance. Therefore it should be noted that slightly different inputs to
this method can result in the same output.

Extra Coordinates
^^^^^^^^^^^^^^^^^

`~sunraster.Raster` allows users to attach additional array-based real world 
coordinates that aren't described by the WCS object via its 
`~sunraster.Raster.extra_coords`.
This can be particularly useful for rastering data because the real world 
coordinates of the slit step axis are actually a convolution of spatial and 
temporal.
Therefore, if the WCS object only supplies (lat, lon) for the x-axis, the 
timestamp of each exposure can be attached as an array of times in a 
`astropy.time.Time` object.
`~sunraster.Raster.extra_coords` is not restricted to time.
The user can supply any additional coordinate as an array or 
`astropy.units.Quantity` or array-like.

Extra coordinates can be supplied during the initiation of a `~sunraster.Raster` 
instance as an iterable of tuples of the form (`str`, `int`,
`~astropy.units.Quantity` or array-like). 
The 0th entry gives the name of the coordinate, the 1st entry gives the data
axis to which the extra coordinate corresponds, and the 2nd entry
gives the value of that coordinate at each pixel along the axis. 
Note that the coordinate array must be the same length as its corresponding 
data axis.
So to add timestamps along the 0th axis of ``my_raster`` we do::

  >>> from datetime import datetime, timedelta
  >>> # Define our timestamps.  Must be same length as data axis.
  >>> axis_length = int(my_raster.dimensions[0].value)
  >>> timestamps = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
  ...                    for i in range(axis_length)], format='datetime', scale='utc')
  >>> extra_coords_input = [("time", 0, timestamps)]
  >>> # Generate Raster as above, except now set extra_coords kwarg.
  >>> my_raster = Raster(data, input_wcs, uncertainty=np.sqrt(data),
  ...                    mask=mask, meta=meta, unit=None,
  ...                    extra_coords=extra_coords_input)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

The `~sunraster.Raster.extra_coords` property returns a dictionary where each key
is a coordinate name entered by the user.  The value of each key is
itself another dictionary with keys ``'axis'`` and ``'value'`` giving the
corresponding data axis number and coordinate value at each pixel as
supplied by the user::

  >>> my_raster.extra_coords # doctest: +SKIP
  {'time': {'axis': 0, 'value': <Time object: scale='utc' format='datetime' value=[datetime.datetime(2000, 1, 1, 0, 0) datetime.datetime(2000, 1, 1, 0, 1) datetime.datetime(2000, 1, 1, 0, 2)]>}}

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data, 
`~sunraster.Raster` inherits a simple-to-use, yet powerful plotting method from 
`ndcube.NDCube`, `~sunraster.Raster.plot`.
It is intended to be a useful quicklook tool and not a
replacement for high quality plots or animations, e.g. for
publications.  The plot method can be called very simply, like so::

  >>> my_raster.plot() # doctest: +SKIP

The type of visualization returned depends on the dimensionality of
the data within the `~sunraster.Raster` object.  For 1-D data a line plot
is produced, similar to `matplotlib.pyplot.plot`.  For 2-D data, an
image is produced similar to that of `matplotlib.pyplot.imshow`.
While for a >2-D data, a
`sunpy.visualization.imageanimator.ImageAnimatorWCS` object is
returned.  This displays a 2-D image with sliders for each additional
dimension which allow the user to animate through the different values
of each dimension and see the effect in the 2-D image.

No args are required.  The necessary information to generate the plot
is derived from the data and metadata in the `~sunraster.Raster`
itself. Setting the x and y ranges of the plot can be done simply by
indexing the `~sunraster.Raster` object itself to the desired region of
interest and then calling the plot method, e.g.::

  >>> my_raster[0, 1:, :].plot() # doctest: +SKIP

In addition, some optional kwargs can be used to customize the
plot.  The ``axis_ranges`` kwarg can be used to set the axes ticklabels.  See the
`~sunpy.visualization.imageanimator.ImageAnimatorWCS` documentation for
more detail.  However, if this is not set, the axis ticklabels are
automatically derived in real world coordinates from the WCS object
within the `~sunraster.Raster`.

By default the final two data dimensions are used for the plot
axes in 2-D or greater visualizations, but this can be set by the user
using the ``images_axes`` kwarg::

  >>> my_raster.plot(image_axes=[0,1]) # doctest: +SKIP

where the first entry in the list gives the index of the data index to
go on the x-axis, and the second entry gives the index of the data
axis to go on the y-axis.

In addition, the units of the axes or the data can be set by the
``unit_x_axis``, ``unit_y_axis``, unit kwargs.  However, if not set,
these are derived from the `~sunraster.Raster` wcs and unit attributes.


