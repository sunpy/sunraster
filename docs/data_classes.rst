.. _data_classes:

============
Data Classes
============

.. _raster:

Raster
------

The fundamental data class of the ``sunraster`` package is `~sunraster.Raster`. 
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
same attributes for data, wcs, extra_coords, uncertainty, mask, meta, and unit.
It also inherits much of the same slicing, WCS, and visualization API and functionality.
In addition, `~sunraster.Raster` provides some additional convenience 
properties relevant to its specific data type.

Initialization
^^^^^^^^^^^^^^
To initialize a basic `~sunraster.Raster` object, all you need is a
3-D `numpy.ndarray` containing the data, and an `astropy.wcs.WCS` object
describing the transformation from array-element (or pixel) space to real
world coordinates.
Let's create an array of data representing a single raster scan
across as region of the Sun.
Let the array shape be (3, 4, 5) and let every value be 1.
Let the first axis represents the slit steps, the second represent the
pixels along the slit, and the third represent the spectral axis.
Although a WCS object can often be easily created by feeding a FITS header into
the `astropy.wcs.WCS` initializer, here we will create one manually to be
explicit.
Note that due to (confusing) convention, the order of the axes in the
WCS object is reversed relative to the data array

.. code-block:: python

  >>> import numpy as np
  >>> data = np.ones((3, 4, 5))
  >>> import astropy.wcs
  >>> wcs_input_dict = {
  ...     'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 5,
  ...     'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 4,
  ...     'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 3}
  >>> input_wcs = astropy.wcs.WCS(wcs_input_dict)

Here we have assigned the first (slit steps) and second (position along slit)
axes to have real world coordinates of helioprojective longitude and latitude.
Although we often think of the x-dimension as longitude and the
y-dimension as latitude, latitude and longitude are in fact coupled dimensions.
This means that except in a small number of edge cases, moving along the slit
in y-direction will cause both the latitude and longitude to change, even if
only slightly.
This is important to understand when interacting with the WCS object,
and hence the `~sunraster.Raster` class.
The 3rd (spectral) axis we have assigned coordinates of wavelength.

Now that we have a data array and a corresponding WCS object, we can
create an `~sunraster.Raster` instance simply by doing:

.. code-block:: python

  >>> from sunraster import Raster
  >>> my_raster = Raster(data, input_wcs)

The data array is stored in the ``my_raster.data`` attribute while the
WCS object is stored in the ``my_raster.wcs`` attribute.  However, when
manipulating/slicing the data is it better to slice the object as a
whole as all relevant data and metadata is sliced simultaneously.
(See section on :ref:`raster_slicing`.)

Thanks to the fact that `~sunraster.Raster` is subclassed from
`~ndcube.NDCube`, you can also supply additional data to the instance. 
These include: metadata (`dict` or dict-like) located at `Raster.meta`;
a data mask (boolean `numpy.ndarray`) located at `Raster.mask` marking, for
example, reliable and unreliable pixels; 
an uncertainty array (`numpy.ndarray`) located at `Raster.uncertainty`
describing the uncertainty of each data array value;
and a unit (`astropy.units.Unit` or unit `str`).
For example:

.. code-block:: python

  >>> mask = np.zeros_like(my_raster.data, dtype=bool)
  >>> meta = {"Description": "This is example Raster metadata."}
  >>> my_raster = Raster(data, input_wcs, uncertainty=np.sqrt(data),
  ...                    mask=mask, meta=meta, unit=None)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

N.B. The above warning is due to the fact that
`astropy.nddata.uncertainty` is recommended to have an
``uncertainty_type`` attribute giving a string describing the type of
uncertainty.
However, this is not currently required by ``sunraster``.

Coordinates
^^^^^^^^^^^

As seen above, coordinate information of a `~sunraster.Raster` instance is
stored in the WCS object.
The coordinate values for each axis and pixel can be accessed via the
`~sunraster.Raster.axis_world_coords method` inherited from `ndcube.NDCube`.
To learn how to use this method, see the
`NDCube.axis_world_coords documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#coordinate-transformations>`_

However, `sunraster.Raster` can also store additional array-based real world
coordinates that aren't described by the WCS object.
These are stored in the `sunraster.Raster.extra_coords` property, also inherited
from `~ndcube.NDCube`.

Extra Coordinates
*****************

`sunraster.Raster.extra_coords` is particularly useful for rastering data
because the real world coordinates of the slit step axis are actually a
convolution of spatial and temporal.
Therefore, if the WCS object only supplies (lat, lon) for the x-axis, the
timestamp of each exposure can be attached as an array of times, e.g. as an
`astropy.time.Time` object.
`~sunraster.Raster.extra_coords` is not restricted to timestamps.
The user can supply any additional coordinate as an `astropy.units.Quantity`
or other array-like.
Metadata that has a relationship with an axis but isn't strictly a coordinate
can also be stored, e.g. exposure time of each image.
(See :ref:`exposure_time_correction` for more on `~sunraster.Raster`'s
handling of exposure times.)

Extra coordinates can be supplied during the initiation of a `~sunraster.Raster`
instance as an iterable of tuples of the form
(`str`, `int`, `~astropy.units.Quantity` or array-like).
The 0th entry gives the name of the coordinate, the 1st entry gives the data
axis to which the extra coordinate corresponds, and the 2nd entry
gives the value of that coordinate at each pixel along the axis.
Note that the coordinate array must be the same length as its corresponding
data axis.
So to add timestamps along the 0th axis of ``my_raster`` we do:

.. code-block:: python

  >>> from datetime import datetime, timedelta
  >>> from astropy.time import Time
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
supplied by the user.

.. code-block:: python

  >>> my_raster.extra_coords # doctest: +SKIP
  {'time': {'axis': 0, 'value': <Time object: scale='utc' format='datetime' value=[datetime.datetime(2000, 1, 1, 0, 0) datetime.datetime(2000, 1, 1, 0, 1) datetime.datetime(2000, 1, 1, 0, 2)]>}}

Coordinate Properties
*********************

For convenience, `~sunraster.Raster` provides shortcuts to the the four primary
coordinates that define raster data.  These are `sunraster.Raster.lon`,
`sunraster.Raster.lat`, `sunraster.Raster.spectral`, and `sunraster.Raster.time`.
These return `~astropy.units.Quantity` objects (or an `astropy.time.Time` object
in the case of `sunraster.Raster.time`) giving the relevant coordinate value of
each pixel.
Note that both `sunraster.Raster.lon` and `sunraster.Raster.lat` return 2-D
`~astropy.units.Quantity` objects because longitude and latitude are couple dimensions.
These properties inspect the WCS and extra coords objects and locate where and
how the relevant coordinate information is stored.
This is possible only if the coordinate name is in the set of those supported
by `sunraster`.
This see these supported names, see ``sunraster.raster.SUPPORTED_LONGITUDE_NAMES``,
``sunraster.raster.SUPPORTED_LATITUDE_NAMES``,
``sunraster.raster.SUPPORTED_SPECTRAL_NAMES``, and
``sunraster.raster.SUPPORTED_TIME_NAMES``.
If the coordinate name cannot be found, these methods raise an error.

In addition to the four primary coordinates, `~sunraster.Raster` also provides a
convenience for the exposure time, `sunraster.Raster.exposure_time`.
The supported exposure time coordinate names can be found under
``sunraster.raster.SUPPORTED_EXPOSURE_NAMES``.


Dimensions
^^^^^^^^^^

The `~sunraster.Raster.dimensions` and
`~sunraster.Raster.world_axis_physical_types` methods on `~sunraster.Raster`
enable users to inspect the shape and axis types of the
`~sunraster.Raster` instance.

.. code-block:: python

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

.. _raster_slicing:

Slicing
^^^^^^^

Users can slice a `~sunraster.Raster` instance using the standard slicing notation.
For example

.. code-block:: python

  >>> my_raster_roi = my_raster[1:3, 0:3, 2:4]

This enables users to access sub-regions of their data while simultaneously
slicing all relevent attributes including uncertainty, mask, wcs, extra_coords, etc.
This ensures that even though the data array has changed size and shape,
each array element will still correspond to the same coordinates, uncertainty
and mask value as before.
Slicing in this way is therefore recommended over accessing the `sunraster.Raster.data`
array directly unless the user wants to access the data values without
its supporting metadata.

Just as with arrays, slicing can reduce the dimensionality of a
`~sunraster.Raster` instance.

.. code-block:: python

  >>> my_2d_raster = my_raster[1:3, 0:3, 2:4]
  >>> my_2d_raster.dimensions
  <Quantity [2., 3., 2.] pix>

In addition to slicing by index, `~sunraster.Raster` supports a basic version of
slicing by real world coordinates via the `~sunraster.Raster.crop_by_coords` method.
This takes a list of `astropy.units.Quantity` instances representing the minimum
real world coordinates of the region of interest in each dimension.
The order of the coordinates must be the same as the order of the data axes.
A second iterable of `~astropy.units.Quantity` must also be provided 
which gives the maximum real world coordinates of the region of interest 
in each data axis.

.. code-block:: python

  >>> import astropy.units as u
  >>> my_raster_roi = my_raster.crop_by_coords([0.7*u.deg, 1.3e-5*u.deg, 1.04e-9*u.m],
  ...                                          [1.3*u.deg, 1.000013*u.deg, 1.12e-9*u.m])

This method does not rebin or interpolate the data if the region of interest
does not perfectly map onto the array's "pixel" grid.  Instead
it returns the smallest rectangular region in pixel space that includes the
entire region defined by the supplied real world coordinates
Therefore it should be noted that slightly different inputs to
this method can result in the same output.

.. _raster_plotting:

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data, 
`~sunraster.Raster` inherits a simple-to-use, yet powerful plotting method from 
`ndcube.NDCube`, `~sunraster.Raster.plot`.
It is intended to be a useful quicklook tool and not a
replacement for high quality plots or animations, e.g. for
publications.  The plot method can be called very simply.

.. code-block:: python

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
interest and then calling the plot method, e.g.:

.. code-block:: python

  >>> my_raster[0, 1:, :].plot() # doctest: +SKIP

In addition, some optional kwargs can be used to customize the
plot.  The ``axis_ranges`` kwarg can be used to set the axes ticklabels.  See the
`~sunpy.visualization.imageanimator.ImageAnimatorWCS` documentation for
more detail.  However, if this is not set, the axis ticklabels are
automatically derived in real world coordinates from the WCS object
within the `~sunraster.Raster`.

By default the final two data dimensions are used as the plot axes.
But this can be set by the user using the ``images_axes`` kwarg:

.. code-block:: python

  >>> my_raster.plot(image_axes=[0, 1]) # doctest: +SKIP

where the first entry in the list gives the index of the data index to
go on the x-axis, and the second entry gives the index of the data
axis to go on the y-axis.

In addition, the units of the axes or the data can be set by the
``unit_x_axis``, ``unit_y_axis``, unit kwargs.  However, if not set,
these are derived from the `~sunraster.Raster` wcs and unit attributes.

.. _raster_exposure_time_correction:

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

An important step in analyzing raster observations is normalizing the data to
the exposure time.
This is important both for converting between instrumental and physical units,
e.g. DN to energy, and when comparing spectral features like line intensity
between exposures and other instruments.

`~sunraster.Raster` provides a simple API for performing this correction:
`~sunraster.Raster.apply_exposure_time_correction`.
It requires that the exposure time is stored the WCS or as a `~astropy.units.Quantity`
in the extra_coords property.
Let's recreate our raster object, but this time with exposure times of 0.5 seconds
and a data unit of counts.

.. code-block:: python

  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> extra_coords_input = [("time", 0, timestamps), ("exposure time", 0, exposure_times)]
  >>> my_raster = Raster(data, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                    meta=meta, unit=u.ct, extra_coords=extra_coords_input)

To apply the exposure time correction simply do:

.. code-block:: python

  >>> # Check the data unit and average data value before applying correction.
  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0
  >>> my_raster = my_raster.apply_exposure_time_correction() # Apply exposure time correction.
  >>> # Check data unit and average data value again.
  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s 2.0

Notice that the average data value has been doubled and the data unit is now counts per second.
This method alters not only the data, but also the uncertainty if any is supplied.
`~sunraster.Raster.apply_exposure_time_correction` does not apply the scaling blindly,
but first checks whether there is a per second (1/s) component in the data unit.
If there is it assumed that correction has already been performed and raises an error.
This helps users more easily keep track of whether they have applied the correction.
However, if for some reason there is a per second component that doesn't refer to the
exposure time and the user still wants to apply the correction, they can set
the ``force`` kwarg to override the check.

.. code-block:: python

  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s 2.0
  >>> my_raster = my_raster.apply_exposure_time_correction(force=True)
  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s / s 4.0

Should users like to undo the correction, they can set the ``undo`` kwarg.

.. code-block:: python

  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s / s 4.0
  >>> my_raster = my_raster.apply_exposure_time_correction(undo=True)
  >>> my_raster = my_raster.apply_exposure_time_correction(undo=True) # Undo correction twice.
  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0

As before, `~sunraster.Raster.apply_exposure_time_correction` only undoes the
correction if there is a time component in the unit.
Again as before, users can override this check by setting the ``force`` kwarg.

.. code-block:: python

  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0
  my_raster = my_raster.apply_exposure_time_correction(undo=True, force=True)
  >>> print(my_raster.unit, my_raster.data.mean())
  ct s 0.5

.. _rastersequence:

RasterSequence
--------------

The `~sunraster.RasterSequence` class is designed to handle multiple raster scans,
where each raster scan is described by a `~sunraster.Raster` object.
It inherits from `ndcube.NDCubeSequence` and whose API is particularly useful
for interpreting the additional dimensionality.
(See section on :ref:`sequence_dimensions`.)

Initialization
^^^^^^^^^^^^^^

To initialize a `~sunraster.RasterSequence`, we first need multiple raster scans
stored in `~sunraster.Raster` instances.
Let's create some using what we learned in the :ref:`raster` section and include
timestamps and exposure times as extra coords.

.. code-block:: python

  >>> import numpy as np
  >>> import astropy.wcs
  >>> import astropy.units as u
  >>> from datetime import datetime, timedelta
  >>> from astropy.time import Time
  >>> from sunraster import Raster

  >>> # Define primary data array and WCS object.
  >>> data = np.ones((3, 4, 5))
  >>> wcs_input_dict = {
  ...     'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 5,
  ...     'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 4,
  ...     'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 3}
  >>> input_wcs = astropy.wcs.WCS(wcs_input_dict)
  >>> mask = np.zeros_like(data.shape[0], dtype=bool)
  >>> meta = {"Description": "This is example Raster metadata."}

  >>> # Define exposure times.
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> axis_length = int(data.shape[0])

  >>> # Create 1st raster
  >>> timestamps0 = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
  ...                     for i in range(axis_length)], format='datetime', scale='utc')
  >>> extra_coords_input0 = [("time", 0, timestamps0), ("exposure time", 0, exposure_times)]
  >>> raster0 = Raster(data, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input0)

  >>> # Create 2nd raster
  >>> timestamps1 = Time([timestamps0[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input1 = [("time", 0, timestamps1), ("exposure time", 0, exposure_times)]
  >>> raster1 = Raster(data*2, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input1)

  >>> # Create 3rd raster
  >>> timestamps2 = Time([timestamps1[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input2 = [("time", 0, timestamps2), ("exposure time", 0, exposure_times)]
  >>> raster2 = Raster(data*0.5, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input2)

If we choose, we can define some sequence-level metadata in addition to any
metadata attached to the individual raster scans:

.. code-block:: python

  >>> seq_meta = {"description": "This is a RasterSequence."}

The last thing we need to do is to identity the ``slit_step_axis`` of each
`~sunraster.Raster`.
This required to correctly handle to convolution of x-position and time.
While `~sunraster.Raster` does not require the physical axes to be in any
particular order, the current implementation of `~sunraster.RasterSequence`
does require that the ``slit_step_axis`` be in the same for each `~sunraster.Raster`.
In our example, the 0th axis of each `~sunraster.Raster` corresponds to time/slit step.
So ``slit_step_axis = 0``.
We can now define our `~sunraster.RasterSequence` by doing:

.. code-block:: python

  >>> from sunraster import RasterSequence
  >>> my_sequence = RasterSequence([raster0, raster1, raster2], meta=seq_meta, slit_step_axis=0)

Dimensions
^^^^^^^^^^

Because the x-axis of raster data corresponds to both space and time, a
`~sunraster.RasterSequence` can be thought of as either 4-D and 3-D.
In the 4-D case, the dimensions represent ``scan number``, ``slit step``,
``position along slit`` and ``spectral axis``.
In the 3-D case the dimensions are ``time``, ``position along slit`` and
``spectral axis``.
Both are perfectly valid representations and do not change the underlying data.
Instead it affects the data's relationship with the real world coordinates and
the way in which users may want to index the data.
Moreover, users may want to switch back anf forth between the different
representations depending on the specific task within their workflow.
To faciliate this,  `~sunraster.RasterSequence` has two ways in which to inspect
the lengths and physical axis types of the data.

The `sunraster.RasterSequence.dimensions`, analagous to `sunraster.Raster.dimensions`,
shows the lengths of the dimensions in the 4-D case:

.. code-block:: python

  >>> my_sequence.dimensions
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)

A tuple of `astropy.units.Quantity` instances with pixel units is returned
giving the length of each axis.
This is in constrast to the single `~astropy.units.Quantity` returned by
`~sunraster.Raster`.
The reason is that `~sunraster.RasterSequence` supports sub-cubes of different
lengths along the ``slit_step_axis``.
In that case, the length of the corresponding quantity in the dimensions tuple
will be equal to the number of raster scans and give the length of each sub-cube
along the ``slit_step_axis``.

Users can view the dimensions in the 3-D case by using the
`sunraster.RasterSequence.cube_like_dimensions`.

.. code-block:: python

  >>> my_sequence.cube_like_dimensions
  <Quantity [9., 4., 5.] pix>

A single `~astropy.units.Quantity` in pixel units is returned giving the
lengths of each dimension.
Because in this representation, the slit step and scan number axes are combined,
there is no need to return a separate `~astropy.units.Quantity` for each dimension.

To view the physical axis types in the 4-D, use
`sunraster.RasterSequence.world_axis_physical_types`.

.. code-block:: python

  >>> my_sequence.world_axis_physical_types
  ('meta.obs.sequence', 'custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'em.wl')

This returns a tuple of the same `IVOA UCD1+ controlled words
<http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>` that are
used by `sunraster.Raster.world_axis_physical_types`.
The sequence axis is given the label ``'meta.obs.sequence'`` as it is the
IVOA UCD1+ controlled word that best describes it.

To view the physical axis types in the 3-D representation, users can employ the
`sunraster.Raster.world_axis_physical_types` method.

.. code-block:: python

  >>> my_sequence.cube_like_world_axis_physical_types  # doctest: +SKIP

Coordinates
^^^^^^^^^^^

Coordinate Properties
*********************

Just like `~sunraster.Raster`, `~sunraster.RasterSequence` provides convenience
properties to retrieve the real world coordinate values for each pixel along
each axis, namely `sunraster.Raster.lon`, `sunraster.Raster.lat`,
`sunraster.Raster.spectral`, `sunraster.Raster.time` and `sunraster.Raster.exposure_time`.
Since their is no guarantee that each `~sunraster.Raster`'s WCS transformations
are consistent between scans, `sunraster.Raster.lon` and `sunraster.Raster.lat`
return 3-D `~astropy.units.Quantity` instances and `sunraster.Raster.spectral` returns a
2-D `~astropy.units.Quantity` where the additional dimension represent the
coordinates for different raster scans.
Meanwhile, since time is sequential across raster scans, both
`sunraster.Raster.time` and `sunraster.Raster.exposure_time` return 1-D
`~astropy.time.Time` and `~astropy.units.Quantity` instance, respectively,
each of the same length as the 0th element of the output of
`~sunraster.RasterSequence..cube_like_dimensions`.


Time Axis Extra Coordinates
***************************

Still to write.

Scan Number Extra Coordinates
*****************************

Still to write.

Slicing
^^^^^^^

Still to write.

Plotting
^^^^^^^^

Still to write.
