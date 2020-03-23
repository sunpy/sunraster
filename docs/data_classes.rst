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
(For data that is described by multiple WCS objects, see the
:ref:`rastersequence` section below.)
`~sunraster.Raster`'s expected axis types are time/slit step, position along slit, 
and spectral.
No specific ordering is required.

`~sunraster.Raster` is subclassed from `ndcube.NDCube` and so inherits the 
same attributes for data, wcs, extra_coords, uncertainty, mask, meta, and unit.
It also inherits much of the same slicing, coordinate transformation and
visualization API.
`~sunraster.Raster` also provides some additional convenience
properties relevant to its specific data type.

Initialization
^^^^^^^^^^^^^^
To initialize a basic `~sunraster.Raster` object, all you need is a
3-D `numpy.ndarray` containing the data, and an `astropy.wcs.WCS` object
describing the transformation from array-element (or pixel) space to real
world coordinates.
Let's create an array of data representing a single raster scan
across a region of the Sun.
Let the array shape be (3, 4, 5) and let every value be 1.
Let the first axis represents the slit steps, the second represent the
pixels along the slit, and the third represent the spectral axis.
Although a WCS object can often be easily created by feeding a FITS header into
the `astropy.wcs.WCS` initializer, here we will create one manually to be
explicit.
Note that due to (confusing) convention, the order of the axes in the
WCS object is reversed relative to the data array.

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

WCS Coordinates
***************
As seen above, coordinate information of a `~sunraster.Raster` instance is
stored in the WCS object.
The coordinate values for each axis and pixel can be accessed via the
`~sunraster.Raster.axis_world_coords`, `~sunraster.Raster.pixel_to_world` and
`~sunraster.Raster.world_to_pixel` methods inherited from `ndcube.NDCube`.
To learn how to use the coordinate transformation methods, see the
`NDCube coordinate transformations documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#coordinate-transformations>`_

Extra Coordinates
*****************

`~sunraster.Raster` can also store array-based real world coordinates that
aren't described by the WCS object.
These are stored in the `sunraster.Raster.extra_coords` property, also inherited
from `~ndcube.NDCube`.
`~sunraster.Raster.extra_coords` is particularly useful for rastering data
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
handling of exposure times and for an example of initializing and using
`~sunraster.Raster.extra_coords`.)
To learn how to attach extra coordinates to a `~sunraster.Raster` instance
and how to access them once attached, see the `NDCube extra coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_

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
This is possible only if the coordinate name is supported by `sunraster`.
This see these supported names, see ``sunraster.raster.SUPPORTED_LONGITUDE_NAMES``,
``sunraster.raster.SUPPORTED_LATITUDE_NAMES``,
``sunraster.raster.SUPPORTED_SPECTRAL_NAMES``, and
``sunraster.raster.SUPPORTED_TIME_NAMES``.
If the coordinate name cannot be found, these properties raise an error.

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
represented by the axes.  The axis names are in accordance with the
International Virtual Observatory Alliance (IVOA)
`UCD1+ controlled vocabulary <http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`_.
Here the shape is given in numpy order, not WCS order.

.. _raster_slicing:

Slicing
^^^^^^^

`~sunraster.Raster` inherits a powerful and simple slicing API from `~ndcube.NDCube`.
It enables users to access sub-regions of their data while simultaneously
slicing all relevent attributes including uncertainty, mask, wcs, extra_coords, etc.
Slicing in pixel space is achieved via the standard Python slicing API while a
separate API is provided for cropping a `~sunraster.Raster` instance by real
world coordinates.
See the
`NDCube slicing documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#slicing>`_
to learn more.

.. _raster_plotting:

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data, 
`~sunraster.Raster` inherits a simple-to-use, yet powerful plotting method from 
`ndcube.NDCube`.
It is intended to be a useful quicklook tool and not a
replacement for high quality plots or animations, e.g. for
publications.  The plot method can be called very simply.

.. code-block:: python

  >>> my_raster.plot() # doctest: +SKIP

For more detail on how this plotting method works and how to customize its
output, see the
`NDCube plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#plotting>`_

.. _raster_exposure_time_correction:

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

An important step in analyzing raster observations is normalizing the data to
the exposure time.
This is important both for converting between instrumental and physical units,
e.g. DN to energy, and comparing spectral features like line intensity
between exposures.

`~sunraster.Raster` provides a simple API for performing this correction:
`~sunraster.Raster.apply_exposure_time_correction`.
It requires that the exposure time is stored the WCS or as a `~astropy.units.Quantity`
in the `~sunraster.Raster.extra_coords` property.
Let's recreate our raster object, but this time with exposure times of
0.5 seconds stored as an extra coordinate and a data unit of counts.
Note that the API for supplying extra coordinates is an iterable of
tuples of the form (`str`, `int`, `~astropy.units.Quantity` or array-like).
The 0th entry gives the name of the coordinate, the 1st entry gives the data
axis to which the extra coordinate corresponds, and the 2nd entry
gives the value of that coordinate at each pixel along the axis.
Also note that the coordinate array must be the same length as its corresponding
data axis.
See the
`NDCube extra coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_
for more.

.. code-block:: python

  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> extra_coords_input = [("exposure time", 0, exposure_times)]
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

The `~sunraster.RasterSequence` class inherits from `ndcube.NDCubeSequence`
and is designed to handle multiple raster scans,
where each raster scan is described by a `~sunraster.Raster` object.

Initialization
^^^^^^^^^^^^^^

To initialize a `~sunraster.RasterSequence`, we first need multiple raster scans
stored in `~sunraster.Raster` instances.
Let's create some using what we learned in the :ref:`raster` section and include
timestamps and exposure times as extra coordinates.

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
This required to correctly handle the convolution of x-position and time.
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
Moreover, users may want to switch back and forth between the different
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
This is because `~sunraster.RasterSequence` supports sub-cubes of different
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
<http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`
used by `sunraster.Raster.world_axis_physical_types`.
The sequence axis is given the label ``'meta.obs.sequence'``.

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
Since there is no guarantee that each `~sunraster.Raster`'s WCS transformations
are consistent between scans, `sunraster.Raster.lon` and `sunraster.Raster.lat`
return 3-D `~astropy.units.Quantity` instances and `sunraster.Raster.spectral`
returns a 2-D `~astropy.units.Quantity` where the additional dimension
represent the coordinates for different raster scans.
Meanwhile, since time is sequential across raster scans, both
`sunraster.Raster.time` and `sunraster.Raster.exposure_time` return 1-D
`~astropy.time.Time` and `~astropy.units.Quantity` instance, respectively,
each of the same length as the 0th element of the output of
`~sunraster.RasterSequence.cube_like_dimensions`.

Time Axis Extra Coordinates
***************************

Still to write.

Scan Number Extra Coordinates
*****************************

Still to write.

Slicing
^^^^^^^

Still to write.

4-D Slicing
***********



3-D Slicing
***********



Plotting
^^^^^^^^

Still to write.
