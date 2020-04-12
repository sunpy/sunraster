.. _data_classes:

============
Data Classes
============

.. _raster:

Raster
------

The fundamental data class of the ``sunraster`` package is `~sunraster.SpectrogramCube`.
It is designed to handle data from a single raster scan or sit-and-stare.
`~sunraster.SpectrogramCube` stores its data as a single 3-D array whose
transformations between pixel and real world coordinates are described by
a single astropy WCS (World Coordinate System) object.
(For data that is described by multiple WCS objects, see the
:ref:`rastersequence` section below.)
`~sunraster.SpectrogramCube`'s expected axis types are time/slit step, position along slit,
and spectral.
No specific ordering is required.

`~sunraster.SpectrogramCube` is subclassed from `ndcube.NDCube` and so inherits the
same attributes for data, wcs, extra_coords, uncertainty, mask, meta, and unit.
It also inherits much of the same slicing, coordinate transformation and
visualization API.
`~sunraster.SpectrogramCube` also provides some additional convenience
properties relevant to its specific data type.

Initialization
^^^^^^^^^^^^^^
To initialize a basic `~sunraster.SpectrogramCube` object, all you need is a
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
and hence the `~sunraster.SpectrogramCube` class.
The 3rd (spectral) axis we have assigned coordinates of wavelength.

Now that we have a data array and a corresponding WCS object, we can
create an `~sunraster.SpectrogramCube` instance simply by doing:

.. code-block:: python

  >>> from sunraster import SpectrogramCube
  >>> my_raster = SpectrogramCube(data, input_wcs)

The data array is stored in the ``my_raster.data`` attribute while the
WCS object is stored in the ``my_raster.wcs`` attribute.  However, when
manipulating/slicing the data is it better to slice the object as a
whole as all relevant data and metadata is sliced simultaneously.
(See section on :ref:`raster_slicing`.)

Thanks to the fact that `~sunraster.SpectrogramCube` is subclassed from
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
  >>> my_raster = SpectrogramCube(data, input_wcs, uncertainty=np.sqrt(data),
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
As seen above, coordinate information of a `~sunraster.SpectrogramCube` instance is
stored in the WCS object.
The coordinate values for each axis and pixel can be accessed via the
`~sunraster.SpectrogramCube.axis_world_coords`, `~sunraster.SpectrogramCube.pixel_to_world` and
`~sunraster.SpectrogramCube.world_to_pixel` methods inherited from `ndcube.NDCube`.
To learn how to use the coordinate transformation methods, see the
`NDCube coordinate transformations documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#coordinate-transformations>`_

Extra Coordinates
*****************

`~sunraster.SpectrogramCube` can also store array-based real world coordinates that
aren't described by the WCS object.
These are stored in the `sunraster.SpectrogramCube.extra_coords` property, also inherited
from `~ndcube.NDCube`.
`~sunraster.SpectrogramCube.extra_coords` is particularly useful for rastering data
because the real world coordinates of the slit step axis are actually a
convolution of spatial and temporal.
Therefore, if the WCS object only supplies (lat, lon) for the x-axis, the
timestamp of each exposure can be attached as an array of times, e.g. as an
`astropy.time.Time` object.
`~sunraster.SpectrogramCube.extra_coords` is not restricted to timestamps.
The user can supply any additional coordinate as an `astropy.units.Quantity`
or other array-like.
Metadata that has a relationship with an axis but isn't strictly a coordinate
can also be stored, e.g. exposure time of each image.
(See :ref:`raster_exposure_time_correction` for more on `~sunraster.SpectrogramCube`'s
handling of exposure times and for an example of initializing and using
`~sunraster.SpectrogramCube.extra_coords`.)
To learn how to attach extra coordinates to a `~sunraster.SpectrogramCube` instance
and how to access them once attached, see the `NDCube extra coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_

Coordinate Properties
*********************

For convenience, `~sunraster.SpectrogramCube` provides shortcuts to the the four primary
coordinates that define raster data.  These are `sunraster.SpectrogramCube.lon`,
`sunraster.SpectrogramCube.lat`, `sunraster.SpectrogramCube.spectral`, and `sunraster.SpectrogramCube.time`.
These return `~astropy.units.Quantity` objects (or an `astropy.time.Time` object
in the case of `sunraster.SpectrogramCube.time`) giving the relevant coordinate value of
each pixel.
Note that both `sunraster.SpectrogramCube.lon` and `sunraster.SpectrogramCube.lat` return 2-D
`~astropy.units.Quantity` objects because longitude and latitude are couple dimensions.
These properties inspect the WCS and extra coords objects and locate where and
how the relevant coordinate information is stored.
This is possible only if the coordinate name is supported by `sunraster`.
This see these supported names, see ``sunraster.SpectrogramCube.SUPPORTED_LONGITUDE_NAMES``,
``sunraster.SpectrogramCube.SUPPORTED_LATITUDE_NAMES``,
``sunraster.SpectrogramCube.SUPPORTED_SPECTRAL_NAMES``, and
``sunraster.SpectrogramCube.SUPPORTED_TIME_NAMES``.
If the coordinate name cannot be found, these properties raise an error.

In addition to the four primary coordinates, `~sunraster.SpectrogramCube` also provides a
convenience for the exposure time, `sunraster.SpectrogramCube.exposure_time`.
The supported exposure time coordinate names can be found under
``sunraster.SpectrogramCube.SUPPORTED_EXPOSURE_NAMES``.

Dimensions
^^^^^^^^^^

The `~sunraster.SpectrogramCube.dimensions` and
`~sunraster.SpectrogramCube.world_axis_physical_types` methods on `~sunraster.SpectrogramCube`
enable users to inspect the shape and axis types of the
`~sunraster.SpectrogramCube` instance.

.. code-block:: python

  >>> my_raster.dimensions
  <Quantity [3., 4., 5.] pix>
  >>> my_raster.world_axis_physical_types
  ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'em.wl')

`~sunraster.SpectrogramCube.dimensions` returns an `~astropy.units.Quantity` of
pixel units giving the length of each dimension in the
`~sunraster.SpectrogramCube` while `~sunraster.SpectrogramCube.world_axis_physical_types`
returns an iterable of strings denoting the type of physical property
represented by the axes.  The axis names are in accordance with the
International Virtual Observatory Alliance (IVOA)
`UCD1+ controlled vocabulary <http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`_.
Here the shape is given in numpy order, not WCS order.

.. _raster_slicing:

Slicing
^^^^^^^

`~sunraster.SpectrogramCube` inherits a powerful and simple slicing API from `~ndcube.NDCube`.
It enables users to access sub-regions of their data while simultaneously
slicing all relevent attributes including uncertainty, mask, wcs, extra_coords, etc.
Slicing in pixel space is achieved via the standard Python slicing API while a
separate API is provided for cropping a `~sunraster.SpectrogramCube` instance by real
world coordinates.
See the
`NDCube slicing documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#slicing>`_
to learn more.

.. _raster_plotting:

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data,
`~sunraster.SpectrogramCube` inherits a simple-to-use, yet powerful plotting method from
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

`~sunraster.SpectrogramCube` provides a simple API for performing this correction:
`~sunraster.SpectrogramCube.apply_exposure_time_correction`.
It requires that the exposure time is stored the WCS or as a `~astropy.units.Quantity`
in the `~sunraster.SpectrogramCube.extra_coords` property.
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

  >>> import astropy.units as u
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> extra_coords_input = [("exposure time", 0, exposure_times)]
  >>> my_raster = SpectrogramCube(data, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                    meta=meta, unit=u.ct, extra_coords=extra_coords_input)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

To apply the exposure time correction simply do:

.. code-block:: python

  >>> # Check the data unit and average data value before applying correction.
  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0
  >>> my_raster = my_raster.apply_exposure_time_correction() # Apply exposure time correction.
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
  >>> # Check data unit and average data value again.
  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s 2.0

Notice that the average data value has been doubled and the data unit is now counts per second.
This method alters not only the data, but also the uncertainty if any is supplied.
`~sunraster.SpectrogramCube.apply_exposure_time_correction` does not apply the scaling blindly,
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
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s2 4.0

Should users like to undo the correction, they can set the ``undo`` kwarg.

.. code-block:: python

  >>> print(my_raster.unit, my_raster.data.mean())
  ct / s2 4.0
  >>> my_raster = my_raster.apply_exposure_time_correction(undo=True, force=True)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
  >>> my_raster = my_raster.apply_exposure_time_correction(undo=True) # Undo correction twice.
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0

As before, `~sunraster.SpectrogramCube.apply_exposure_time_correction` only undoes the
correction if there is a time component in the unit.
Again as before, users can override this check by setting the ``force`` kwarg.

.. code-block:: python

  >>> print(my_raster.unit, my_raster.data.mean())
  ct 1.0
  >>> my_raster = my_raster.apply_exposure_time_correction(undo=True, force=True)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
  >>> print(my_raster.unit, my_raster.data.mean())
  ct s 0.5

.. _rastersequence:

RasterSequence
--------------

The `~sunraster.RasterSequence` class inherits from `ndcube.NDCubeSequence`
and is designed to handle multiple raster scans,
where each raster scan is described by a `~sunraster.SpectrogramCube` object.

Initialization
^^^^^^^^^^^^^^

To initialize a `~sunraster.RasterSequence`, we first need multiple raster scans
stored in `~sunraster.SpectrogramCube` instances.
Let's create some using what we learned in the :ref:`raster` section and include
timestamps and exposure times as extra coordinates.

.. code-block:: python

  >>> import numpy as np
  >>> import astropy.wcs
  >>> import astropy.units as u
  >>> from datetime import datetime, timedelta
  >>> from astropy.time import Time
  >>> from sunraster import SpectrogramCube

  >>> # Define primary data array and WCS object.
  >>> data = np.ones((3, 4, 5))
  >>> wcs_input_dict = {
  ...     'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 5,
  ...     'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 4,
  ...     'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 3}
  >>> input_wcs = astropy.wcs.WCS(wcs_input_dict)
  >>> mask = np.zeros(data.shape, dtype=bool)
  >>> meta = {"Description": "This is example Raster metadata."}

  >>> # Define exposure times.
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> axis_length = int(data.shape[0])

  >>> # Create 1st raster
  >>> timestamps0 = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
  ...                     for i in range(axis_length)], format='datetime', scale='utc')
  >>> extra_coords_input0 = [("time", 0, timestamps0), ("exposure time", 0, exposure_times)]
  >>> raster0 = SpectrogramCube(data, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input0)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

  >>> # Create 2nd raster
  >>> timestamps1 = Time([timestamps0[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input1 = [("time", 0, timestamps1), ("exposure time", 0, exposure_times)]
  >>> raster1 = SpectrogramCube(data*2, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input1)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

  >>> # Create 3rd raster
  >>> timestamps2 = Time([timestamps1[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input2 = [("time", 0, timestamps2), ("exposure time", 0, exposure_times)]
  >>> raster2 = SpectrogramCube(data*0.5, input_wcs, uncertainty=np.sqrt(data), mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input2)
  INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]

If we choose, we can define some sequence-level metadata in addition to any
metadata attached to the individual raster scans:

.. code-block:: python

  >>> seq_meta = {"description": "This is a RasterSequence."}

The last thing we need to do is to identity the ``slit_step_axis`` of each
`~sunraster.SpectrogramCube`.
This required to correctly handle the convolution of x-position and time.
While `~sunraster.SpectrogramCube` does not require the physical axes to be in any
particular order, the current implementation of `~sunraster.RasterSequence`
does require that the ``slit_step_axis`` be in the same for each `~sunraster.SpectrogramCube`.
In our example, the 0th axis of each `~sunraster.SpectrogramCube` corresponds to time/slit step.
So ``slit_step_axis = 0``.
We can now define our `~sunraster.RasterSequence` by doing:

.. code-block:: python

  >>> from sunraster import RasterSequence
  >>> my_sequence = RasterSequence([raster0, raster1, raster2], meta=seq_meta, common_axis=0)

Dimensions
^^^^^^^^^^

Because the x-axis of raster data corresponds to both space and time, a
`~sunraster.RasterSequence` can be thought of as either 4-D and 3-D.
In the 4-D, or raster, case, the dimensions represent
``scan number``, ``slit step``, ``position along slit`` and ``spectral axis``.
In the 3-D, or sit-and-stare, case the dimensions represent
``time``, ``position along slit`` and ``spectral axis``.
In this guide we use sit-and-stare (SnS) to refer to this 3-D way of representing
the data, regardless of whether the slit is scanning.
Both the raster and sit-and-stare representations are perfectly valid and do
not change the underlying data.
Instead they affect the data's relationship with the real world coordinates and
the way in which users may want to index the data.
Moreover, users may want to switch back and forth between the different
representations depending on the specific task within their workflow.
To faciliate this,  `~sunraster.RasterSequence` has two ways in which to inspect
the lengths and physical axis types of the data.

The `sunraster.RasterSequence.raster_dimensions`, analagous to
`sunraster.SpectrogramCube.dimensions`, shows the lengths of the dimensions in the
4-D case:

.. code-block:: python

  >>> my_sequence.raster_dimensions
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)

A tuple of `astropy.units.Quantity` instances with pixel units is returned
giving the length of each axis.
This is in constrast to the single `~astropy.units.Quantity` returned by
`~sunraster.SpectrogramCube`.
This is because `~sunraster.RasterSequence` supports sub-cubes of different
lengths along the ``slit_step_axis``.
In that case, the length of the ``slit_step_axis`` quantity will be equal to
the number of raster scans and give the number of slit steps in each raster.

Users can view the dimensions in the sit-and-stare representation by using
the `sunraster.RasterSequence.SnS_dimensions` property.

.. code-block:: python

  >>> my_sequence.SnS_dimensions
  <Quantity [9., 4., 5.] pix>

A single `~astropy.units.Quantity` in pixel units is returned giving the
lengths of each dimension.
Because in this representation, the slit step and scan number axes are combined,
there is no need to return a separate `~astropy.units.Quantity` for each dimension.

To view the physical axis types in the raster representation, use
`sunraster.RasterSequence.raster_world_axis_physical_types`.

.. code-block:: python

  >>> my_sequence.raster_world_axis_physical_types
  ('meta.obs.sequence', 'custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'em.wl')

This returns a tuple of the same `IVOA UCD1+ controlled words
<http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`
used by `sunraster.SpectrogramCube.world_axis_physical_types`.
The sequence axis is given the label ``'meta.obs.sequence'``.

To view the physical axis types in the sit-and-stare representation, users
can employ the `sunraster.SpectrogramCube.world_axis_physical_types` method.

.. code-block:: python

  >>> my_sequence.SnS_world_axis_physical_types  # doctest: +SKIP

Coordinates
^^^^^^^^^^^

Coordinate Properties
*********************

Just like `~sunraster.SpectrogramCube`, `~sunraster.RasterSequence` provides convenience
properties to retrieve the real world coordinate values for each pixel along
each axis, namely `sunraster.SpectrogramCube.lon`, `sunraster.SpectrogramCube.lat`,
`sunraster.SpectrogramCube.spectral`, `sunraster.SpectrogramCube.time` and `sunraster.SpectrogramCube.exposure_time`.
Since there is no guarantee that each `~sunraster.SpectrogramCube`'s WCS transformations
are consistent between scans, `sunraster.SpectrogramCube.lon` and `sunraster.SpectrogramCube.lat`
return 3-D `~astropy.units.Quantity` instances and `sunraster.SpectrogramCube.spectral`
returns a 2-D `~astropy.units.Quantity` where the additional dimension
represent the coordinates for different raster scans.
Meanwhile, since time is sequential across raster scans, both
`sunraster.SpectrogramCube.time` and `sunraster.SpectrogramCube.exposure_time` return 1-D
`~astropy.time.Time` and `~astropy.units.Quantity` instance, respectively,
each of the same length as the 0th element of the output of
`~sunraster.RasterSequence.cube_like_dimensions`.

SnS Axis Extra Coordinates
**************************

As well as `sunraster.SpectrogramCube.time` and `sunraster.SpectrogramCube.exposure_time`,
`sunraster.SpectrogramCube.extra_coords` may contain other coordinates that are aligned
with the time/slit step axis.
The `sunraster.RasterSequence.SnS_axis_extra_coords` property enables users
to access these coordinates at the `~sunraster.RasterSequence` level in the
form of an abbreviated ``extra_coords`` dictionary.
Just like `sunraster.SpectrogramCube.time` and `sunraster.SpectrogramCube.exposure_time`,
the coordinates are concatenated so they mimic the sit-and-stare-like dimensionality
returned in the 0th element of `sunraster.RasterSequence.SnS_dimensions`.
`sunraster.RasterSequence.SnS_axis_extra_coords` is equivalent to
`ndcube.NDCubeSequence.common_axis_extra_coords`.
To see examples of how to use this property, see the
`NDCubeSequence Common Axis Extra Coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#common-axis-extra-coordinates>`_.

Raster Axis Extra Coordinates
*****************************

Analgous to `~sunraster.RasterSequence.SnS_axis_extra_coords`, it is also
possible to access the extra coordinates that are not assigned to any
`~sunraster.SpectrogramCube` data axis via the
`~sunraster.RasterSequence.raster_axis_extra_coords` property.
Whereas `~sunraster.RasterSequence.SnS_axis_extra_coords` returns all the
extra coords with an ``'axis'`` value equal to the time/slit step axis,
`~sunraster.RasterSequence.scan_axis_extra_coords` returns all extra coords
with an ``'axis'`` value of None.
Another way of thinking about this is that these coordinates correspond to
the repeat raster scan number axis.
Hence the property’s name.

Slicing
^^^^^^^

`~sunraster.RasterSequence` not only enables users to inspect their data in
the raster and sit-and-stare representations.
It also enables users to slice the data in either representation as well.
This is doen via the `~sunraster.RasterSequence.slice_as_raster` and
`~sunraster.RasterSequence.slice_as_SnS` properties.
As with `~sunraster.SpectrogramCube`, these slicing properties ensure that not only
the data is sliced, but all relevant supporting metadata as well including
incertainties, mask, WCS object, extra_coords, etc.

To slice a `~sunraster.RasterSequence` using the raster representation, do
the following:

.. code-block:: python

  >>> print(my_sequence.raster_dimensions)  # Check dimensionality before slicing.
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)

  >>> my_sequence_roi = my_sequence.slice_as_raster[1:3, 0:2, 1:3, 1:4]

  >>> print(my_sequence_roi.raster_dimensions) # See how slicing has changed dimensionality.
  (<Quantity 2. pix>, <Quantity 2. pix>, <Quantity 2. pix>, <Quantity 3. pix>)
  >>> my_sequence_roi.SnS_dimensions  # Dimensionality can still be represented in SnS form.
  <Quantity [4., 2., 3.] pix>

To slice in the sit-and-stare representation, do the following:

.. code-block:: python

  >>> print(my_sequence.SnS_dimensions)  # Check dimensionality before slicing.
  [9. 4. 5.] pix

  >>> my_sequence_roi = my_sequence.slice_as_SnS[1:7, 1:3, 1:4]

  >>> print(my_sequence_roi.SnS_dimensions)  # See how slicing has changed dimensionality.
  [6. 2. 3.] pix
  >>> print(my_sequence_roi.raster_dimensions)  # Dimensionality can still be represented in raster form.
  (<Quantity 3. pix>, <Quantity [2., 3., 1.] pix>, <Quantity 2. pix>, <Quantity 3. pix>)

Notice that after slicing the data can still be inspected and interpreted in
the raster or sit-and-stare format, irrespective of which slicing
representation was used.
Also notice that the ``my_sequence.slice_as_SnS[1:7, 1:3, 1:4]`` command led
different `~sunraster.SpectrogramCube` objects to have different lengths along the
slit step axis.
And that this can be seen from the fact that the slit step axis entry in the
output of ``my_sequence_roi.raster_dimensions`` has a length greater than 1.

Slicing can reduce the dimensionality of `~sunraster.RasterSequence` instances.
For example, let's slice out the 2nd pixel along the slit.
This reduces the number of dimensions in the raster representation to 3
(``raster scan``, ``slit step``, ``spectral``) and to 2 in the sit-and-stare
representation (``time``, ``spectral``).
However, the raster and sit-and-stare representations are still valid.

.. code-block:: python

  >>> slit_pixel_sequence = my_sequence.slice_as_raster[:, :, 2]
  >>> print(slit_pixel_sequence.raster_dimensions)
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 5. pix>)
  >>> print(slit_pixel_sequence.SnS_dimensions)
  [9. 5.] pix

This demonstrates the the difference between the raster and sit-and-stare
representations is more subtle than simply a 4-D or 3-D dimensionality.
The difference is whether the raster scan and slit step axes are convolved
into a time axis or whether they are represented separately.
And because of this definition, the raster and sit-and-stare representations
are valid and accessible for any dimensionality in which the raster scan and
slit step axes are maintained.

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data,
`~sunraster.RasterSequence` supplies simple-to-use, yet powerful plotting APIs.
They are intended to be a useful quicklook tool and not a
replacement for high quality plots or animations, e.g. for
publications.
As with slicing, there are two plot methods for plotting in each of the
raster and sit-and-stare representations.

To visualize in the raster representation, simply call the following:

.. code-block:: python

  >>> my_sequence.plot_as_raster() # doctest: +SKIP

To visualize in the sit-and-stare representation, do:

.. code-block:: python

  >>> my_sequence.plot_as_SnS() # doctest: +SKIP

These methods produce different types of visualizations including line plots,
2-D images and 1- and 2-D animations.
Which is displayed depends on the dimensionality of the `~sunraster.RasterSequence`
and the inputs of the user.
`~sunraster.RasterSequence.plot_as_raster` and
`~sunraster.RasterSequence.plot_as_SnS` are in fact simply aliases for the
`ndcube.NDCubeSequence.plot` and `ndcube.NDCubeSequence.plot_as_cube` methods,
respectively.
For learn more about how these routines work and the optional inputs that
enable users to customize their output, see the
`NDCubeSequence plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting>`_.

Extracting Data Arrays
^^^^^^^^^^^^^^^^^^^^^^

It is possible that you may have some procedures that are designed to operate on arrays instead of
`~sunraster.RasterSequence` objects.
Therefore it may be useful to extract the data (or other array-like information
such as `uncertainty` or `mask`) in the `~sunraster.RasterSequence` into a single `~numpy.ndarray`.
A succinct way of doing this operation is using python's list comprehension features.

To make a 4-D array from the data arrays of the `~sunraster.SpectrogramCube` within ``my_sequence``,
use `numpy.stack`.

.. code-block:: python

    >>> print(my_sequence.raster_dimensions)  # Print sequence dimensions as a reminder.
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
    >>> data = np.stack([cube.data for cube in my_sequence.data])
    >>> print(data.shape)
    (3, 3, 4, 5)

To define a 3D array where every `~sunraster.SpectrogramCube` data array in the
`~sunraster.RasterSequence` is appended together, we can use `numpy.vstack`.

.. code-block:: python

    >>> data = np.vstack([cube.data for cube in my_sequence.data])
    >>> print(data.shape)
    (9, 4, 5)

To create 3D arrays by slicing `~sunraster.RasterSequence` objects:

.. code-block:: python

    >>> data = np.stack([cube[2].data for cube in my_sequence.data])
    >>> print(data.shape)
    (3, 4, 5)

Raster Collections
------------------

During analysis of slit spectrograph data, it is often desirable to group
different data sets together together.
For example, you may have several `~sunraster.SpectrogramCube` or
`~sunraster.RasterSequence` objects representing observations in different
spectral windows.
Or we may have fit a spectral line in each pixel and extracted a property
such as linewidth, thus collapsing the spectral axis.
In both these cases, the `~sunraster.RasterSequence` objects share a common
origin and coordinate transformations with the original observations
(except in the spectral axis).
However, they do not have a sequential relationship in their common coordinate spaces
and in the latter case the data represents a different physical property to the
original observations.
Therefore, combining them in a `~sunraster.RasterSequence` is not appropriate.

`sunraster`` does not provide a suitable object for this purpose.
However, because `~sunraster.SpectrogramCube` and `~sunraster.RasterSequence` are
instances of `ndcube` classes underneath, users can employ the `ndcube.NDCollection`
class for this purpose.
`~ndcube.NDCollection` is a `dict`-like class that provides additional slicing
capabilities of its constituent data cubes along aligned axes.
To see whether `~ndcube.NDCollection` could be helpful for your research, see
the
`NDCollection documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcollection.html>`_.
