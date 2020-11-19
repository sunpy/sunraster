.. _data_classes:

============
Data Classes
============

.. _spectrogramcube:

SpectrogramCube
---------------

The fundamental data class of the ``sunraster`` package is `~sunraster.SpectrogramCube`.
It is designed to handle data representing one or more spectrograms of solar regions.
`~sunraster.SpectrogramCube` stores its data as an array whose
transformations between pixel and real world coordinates are described by
a single ``astropy`` WCS (World Coordinate System) object.
(For data that is described by multiple WCS objects, see the
:ref:`sequence` :ref:`raster_sequence` sections below.)
`~sunraster.SpectrogramCube` is subclassed from `ndcube.NDCube` and so inherits the
same attributes for data, wcs, extra_coords, uncertainty, mask, meta, and unit.
It also inherits much of the same slicing, coordinate transformation and
visualization API and provides some additional convenience properties relevant to
spectrogram data.

Initialization
^^^^^^^^^^^^^^
To initialize a basic `~sunraster.SpectrogramCube` object, all you need is an
array containing the data and an `astropy.wcs.WCS` object
describing the transformation from array-element (or pixel) space to real
world coordinates.
Let's create a 3D `numpy.ndarray` representing a series of spectrograms.
Let the array shape be (3, 4, 5) and let every value be 1.
Let the first axis represent time (and/or space if the spectrogram slit
is rastering across a solar region).
Let the second represent the position along a dispersing slit,
and the third represent the spectral axis.
Although a WCS object can often be easily created by feeding a FITS header into
the `astropy.wcs.WCS` class, we will create one manually here to be explicit.
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

We have defined the first axis to be spatial (helioprojective longitude and latitude)
which implies that this series of spectrograms represents a raster scan across a
solar region. The second axis (position along slit) also has coordinates of
helioprojective longitude and latitude.
Although we often think of the x-dimension as longitude and the
y-dimension as latitude, latitude and longitude are in fact coupled dimensions.
This means that -- except in a small number of edge cases -- moving along the slit
in y-direction will cause both the latitude AND longitude to change, even if
only slightly. This is important to understand when interacting
with the WCS object, and hence the `~sunraster.SpectrogramCube` class.
The 3rd axis (spectral) has coordinates of wavelength.

Now that we have a data array and a corresponding WCS object, we can
create a `~sunraster.SpectrogramCube` instance simply by doing:

.. code-block:: python

  >>> from sunraster import SpectrogramCube
  >>> my_spectrograms = SpectrogramCube(data, input_wcs)

The data array is stored in the ``my_spectrograms.data`` attribute while the
WCS object is stored in the ``my_spectrograms.wcs`` attribute.  However, when
manipulating/slicing the data is it better to slice the object as a
whole as all relevant data and metadata is sliced simultaneously.
(See section on :ref:`spectrogram_slicing`.)

Thanks to the fact that `~sunraster.SpectrogramCube` is subclassed from
`~ndcube.NDCube`, you can also supply additional data to the instance.
These include: metadata (``dict`` or dict-like) located in `sunraster.SpectrogramCube.meta`;
a data mask (boolean ``numpy.ndarray``) located in ``sunraster.SpectrogramCube.mask``
for marking reliable and unreliable pixels;
a unit (``astropy.units.Unit`` or unit `str`) located at ``sunraster.SpectrogramCube.unit``;
and an uncertainty array (`numpy.ndarray`) located in `~sunraster.SpectrogramCube.uncertainty`
describing the uncertainty of each data array value. It is advised that you use
one of astropy's uncertainty classes to describe your uncertainty.
However, this is not required by `~sunraster.SpectrogramCube`.
A simple array will still work but will cause a warning to be raised.
Here is an example of how to instantiate these attributes.

.. code-block:: python

  >>> import astropy.units as u
  >>> from astropy.nddata import StdDevUncertainty
  >>> uncertainties = StdDevUncertainty(np.sqrt(data))
  >>> # Create a mask where all pixels are unmasked, i.e. all mask values are False.
  >>> mask = np.zeros_like(data, dtype=bool)
  >>> meta = {"Description": "This is example SpectrogramCube metadata."}
  >>> my_spectrograms = SpectrogramCube(data, input_wcs, uncertainty=uncertainties,
  ...                                   mask=mask, meta=meta)

Coordinates
^^^^^^^^^^^

WCS Coordinates
***************

The primary location for coordinate information in a `~sunraster.SpectrogramCube`
instance is its WCS.
The coordinate values for each axis and pixel can be accessed via the
`~sunraster.SpectrogramCube.axis_world_coords`, `~sunraster.SpectrogramCube.pixel_to_world` and
`~sunraster.SpectrogramCube.world_to_pixel` methods inherited from ``ndcube.NDCube``.
To learn how to use these coordinate transformation methods, see the
`NDCube coordinate transformations documentation
<https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#coordinate-transformations>`_

Extra Coordinates
*****************

`~sunraster.SpectrogramCube` can also store array-based real world coordinates
that aren't described by the WCS object.
These can be accessed via the ``sunraster.SpectrogramCube.extra_coords`` property,
also inherited from `~ndcube.NDCube`.
`~sunraster.SpectrogramCube.extra_coords` is particularly useful if
the temporal axis is convolved with space, as is the case for raster scans.
Therefore, if the WCS object only supplies (lat, lon) for the x-axis, the
timestamp of each exposure can be attached separately, e.g. as an
``astropy.time.Time`` object. `~sunraster.SpectrogramCube.extra_coords`
is not restricted to timestamps. The user can supply any additional coordinate
as an ``astropy.units.Quantity`` or other array-like.
Metadata that has a relationship with an axis but isn't strictly a coordinate
can also be stored, e.g. the exposure time of each image.
(See :ref:`cube_exposure_time_correction` for more on `~sunraster.SpectrogramCube`'s
handling of exposure times.)
To learn how to attach extra coordinates to a `~sunraster.SpectrogramCube` instance
and how to access them once attached, see the
`NDCube extra coordinates documentation
<https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_

Coordinate Properties
*********************

For convenience, `~sunraster.SpectrogramCube` provides shortcuts to the
four primary coordinates that define spectrogram data.
These are `sunraster.SpectrogramCube.lon`, `sunraster.SpectrogramCube.lat`,
`sunraster.SpectrogramCube.spectral`, and `sunraster.SpectrogramCube.time`
which return the relevant coordinate values of each pixel.
Note that both `sunraster.SpectrogramCube.lon` and `sunraster.SpectrogramCube.lat`
return 2-D data because longitude and latitude are couple dimensions.
These properties inspect the WCS and extra coords objects and locate where and
how the relevant coordinate information is stored.
This is possible only if the coordinate name is supported by ``sunraster``.
To see these supported names, see
``sunraster.SpectrogramCube.SUPPORTED_LONGITUDE_NAMES``,
``sunraster.spectrogram.SUPPORTED_LATITUDE_NAMES``,
``sunraster.spectrogram.SUPPORTED_SPECTRAL_NAMES``, and
``sunraster.spectrogram.SUPPORTED_TIME_NAMES``.
If the coordinate name cannot be found, these properties will raise an error.
If you think additional coordinate names should be supported,
please let us know by `raising an issue on our GitHub repo. <https://github.com/sunpy/sunraster/issues>`

In addition to the four primary coordinates, there is also a convenience
for the exposure time, ``sunraster.SpectrogramCube.exposure_time``.
The supported exposure time coordinate names can be found under
``sunraster.spectrogram.SUPPORTED_EXPOSURE_NAMES``.

Dimensions
^^^^^^^^^^

The `~sunraster.SpectrogramCube.dimensions` and
`~sunraster.SpectrogramCube.array_axis_physical_types` methods
enable users to inspect the shape and WCS axis types of the
`~sunraster.SpectrogramCube` instance.

.. code-block:: python

  >>> my_spectrograms.dimensions
  <Quantity [3., 4., 5.] pix>
  >>> my_spectrograms.array_axis_physical_types
  [('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('em.wl',)]

`~sunraster.SpectrogramCube.dimensions` returns a `~astropy.units.Quantity`
giving the length of each dimension in pixel units while
`~sunraster.SpectrogramCube.array_axis_physical_types`
returns an list of tuples where each tuple contains the types of physical properties
associated with each array axis. Since more than one physical type be associated with
an array axis because they are dependent, e.g. latitude/longitude, or because of the
rastering naturing of the instrument, e.g. latitude/longitude and time, the length of each
tuple can be greater than one.
The axis names are in accordance with the International Virtual Observatory Alliance (IVOA)
`UCD1+ controlled vocabulary <http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`_.

.. _spectrogram_slicing:

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

.. _spectrogram_plotting:

Plotting
^^^^^^^^

To quickly and easily visualize spectrograms,
`~sunraster.SpectrogramCube` inherits a simple-to-use,
yet powerful plotting method from `~ndcube.NDCube`.
It is intended to be a useful quicklook tool and not a
replacement for high quality plots or animations, e.g. for
publications.  The plot method can be called very simply.

.. code-block:: python

  >>> my_spectrograms.plot() # doctest: +SKIP

This method produces different types of visualizations including line plots,
2-D images and 1- and 2-D animations.
Which is displayed depends on the dimensionality of the `~sunraster.SpectrogramCube`
and the inputs of the user.
For learn more about how to customize plots and animations through the
`~sunraster.SpectrogramCube.plot` method, see the
`NDCubeSequence plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting>`_.

.. _cube_exposure_time_correction:

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

An important step in analyzing any form of photon-based observations is normalizing
the data to the exposure time.
This is important both for converting between instrumental and physical units,
e.g. DN to energy, and comparing spectral features between exposure, e.g. line intensity.

`~sunraster.SpectrogramCube` provides a simple API for performing this correction:
`~sunraster.SpectrogramCube.apply_exposure_time_correction`.
It requires that the exposure time is stored the WCS or as a `~astropy.units.Quantity`
in the `~sunraster.SpectrogramCube.extra_coords` property.
Let's recreate our spectrogram object again, but this time with exposure times of
0.5 seconds stored as an extra coordinate and a data unit of counts.

.. code-block:: python

  >>> import astropy.units as u
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> extra_coords_input = [("exposure time", 0, exposure_times)]
  >>> my_spectrograms = SpectrogramCube(data, input_wcs, uncertainty=uncertainties,
  ...                                   mask=mask, meta=meta, unit=u.ct,
  ...                                   extra_coords=extra_coords_input)

Note that the API for supplying extra coordinates is an iterable of
tuples of the form (``str``, ``int``, `~astropy.units.Quantity` or array-like).
The 0th entry gives the name of the coordinate, the 1st entry gives the data
axis to which the extra coordinate corresponds, and the 2nd entry
gives the value of that coordinate at each pixel along the axis.
Also note that the coordinate array must be the same length as its corresponding
data axis.
See the
`NDCube extra coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates>`_
for more.

Applying the exposure time correction is now simple.

.. code-block:: python

  >>> # First check the data unit and average data value before applying correction.
  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct 1.0
  >>> my_spectrograms = my_spectrograms.apply_exposure_time_correction() # Apply exposure time correction.
  >>> # Confirm effect by checking data unit and average data value again.
  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct / s 2.0

Notice that the average data value has been doubled and
the data unit is now counts per second.
This method alters not only the data, but also the uncertainty if any is supplied.
`~sunraster.SpectrogramCube.apply_exposure_time_correction`
does not apply the scaling blindly, but first checks whether there is
a per second (1/s) component in the data unit.
If there is, it assumes that the correction has already been performed
and raises an error.
This helps users more easily keep track of whether they have applied the correction.
However, if for some reason there is a per second component that
doesn't refer to the exposure time and the user still wants to apply the correction,
they can set the ``force`` kwarg to override the check.

.. code-block:: python

  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct / s 2.0
  >>> my_spectrograms = my_spectrograms.apply_exposure_time_correction(force=True)
  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct / s2 4.0

Should users like to undo the correction, they can set the ``undo`` kwarg.

.. code-block:: python

  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct / s2 4.0
  >>> my_spectrograms = my_spectrograms.apply_exposure_time_correction(undo=True, force=True)
  >>> my_spectrograms = my_spectrograms.apply_exposure_time_correction(undo=True) # Undo correction twice.
  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct 1.0

As before, `~sunraster.SpectrogramCube.apply_exposure_time_correction` only undoes the
correction if there is a time component in the unit.
And again as before, users can override this check by setting the ``force`` kwarg.

.. code-block:: python

  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct 1.0
  >>> my_spectrograms = my_spectrograms.apply_exposure_time_correction(undo=True, force=True)
  >>> print(my_spectrograms.unit, my_spectrograms.data.mean())
  ct s 0.5

.. _sequence:

SpectrogramSequence
-------------------

In some cases, a series of spectrograms may not be describable by a
single set of WCS transformations.
However, it still may make sense to combine them in order along a dimension.
This is the purpose of the `~sunraster.SpectrogramSequence` class.
It stores a sequence of `~sunraster.SpectrogramCube` instances and provides
equivalent or analagous APIs so users can interact with the data as if it were
a single data cube.
`~sunraster.SpectrogramSequence` inherits from `~ndcube.NDCubeSequence` and
so inherits much of the same API.

Initialization
^^^^^^^^^^^^^^

To initialize a `~sunraster.SpectrogramSequence`, we first need spectrograms
stored in multiple `~sunraster.SpectrogramCube` instances.
Let's create some using what we learned in the :ref:`spectrogramcube` section and include
timestamps and exposure times as extra coordinates.

.. code-block:: python

  >>> from datetime import datetime, timedelta
  >>> import numpy as np
  >>> import astropy.wcs
  >>> import astropy.units as u
  >>> from astropy.nddata import StdDevUncertainty
  >>> from astropy.time import Time
  >>> from sunraster import SpectrogramCube

  >>> # Define primary data array and WCS object.
  >>> data = np.ones((3, 4, 5))
  >>> wcs_input_dict = {
  ...     'CTYPE1': 'WAVE    ', 'CUNIT1': 'Angstrom', 'CDELT1': 0.2, 'CRPIX1': 0, 'CRVAL1': 10, 'NAXIS1': 5,
  ...     'CTYPE2': 'HPLT-TAN', 'CUNIT2': 'deg', 'CDELT2': 0.5, 'CRPIX2': 2, 'CRVAL2': 0.5, 'NAXIS2': 4,
  ...     'CTYPE3': 'HPLN-TAN', 'CUNIT3': 'deg', 'CDELT3': 0.4, 'CRPIX3': 2, 'CRVAL3': 1, 'NAXIS3': 3}
  >>> input_wcs = astropy.wcs.WCS(wcs_input_dict)
  >>> # Define a mask with all pixel unmasked, i.e. mask values = False
  >>> mask = np.zeros(data.shape, dtype=bool)
  >>> # Define uncertaines for data, 2*data and data/2.
  >>> uncertainties = StdDevUncertainty(np.sqrt(data))
  >>> uncertainties2 = StdDevUncertainty(np.sqrt(data * 2))
  >>> uncertainties05 = StdDevUncertainty(np.sqrt(data * 0.5))

  >>> # Define exposure times.
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> axis_length = int(data.shape[0])

  >>> # Create 1st cube of spectrograms.
  >>> timestamps0 = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
  ...                     for i in range(axis_length)], format='datetime', scale='utc')
  >>> extra_coords_input0 = [("time", 0, timestamps0), ("exposure time", 0, exposure_times)]
  >>> spectrograms0 = SpectrogramCube(data, input_wcs, uncertainty=uncertainties, mask=mask,
  ...                                 meta=meta, unit=u.ct, extra_coords=extra_coords_input0)

  >>> # Create 2nd cube of spectrograms.
  >>> timestamps1 = Time([timestamps0[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input1 = [("time", 0, timestamps1), ("exposure time", 0, exposure_times)]
  >>> spectrograms1 = SpectrogramCube(data*2, input_wcs, uncertainty=uncertainties2, mask=mask,
  ...                                 meta=meta, unit=u.ct, extra_coords=extra_coords_input1)

  >>> # Create 3rd cube of spectrograms.
  >>> timestamps2 = Time([timestamps1[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input2 = [("time", 0, timestamps2), ("exposure time", 0, exposure_times)]
  >>> spectrograms2 = SpectrogramCube(data*0.5, input_wcs, uncertainty=uncertainties05, mask=mask,
  ...                                 meta=meta, unit=u.ct, extra_coords=extra_coords_input2)

If we choose, we can define some sequence-level metadata in addition to any
metadata attached to the individual raster scans:

.. code-block:: python

  >>> seq_meta = {"description": "This is a SpectrogramSequence."}

To create a `~sunraster.SpectrogramSequence`, simply supply the class with a
list of `~sunraster.SpectrogramCube` instances.

.. code-block:: python

  >>> from sunraster import SpectrogramSequence
  >>> my_sequence = SpectrogramSequence([spectrograms0, spectrograms1, spectrograms2],
  ...                                   meta=seq_meta)

Dimensions
^^^^^^^^^^

In order to inspect the dimensionlity of our sequence and the physical properties
to which the axes correspond, we can use the
`~sunraster.SpectrogramSequence.dimensions` and
`~sunraster.SpectrogramSequence.array_axis_physical_types` properties.

.. code-block:: python

  >>> my_sequence.dimensions
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
  >>> my_sequence.array_axis_physical_types
  [('meta.obs.sequence',),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('em.wl',)]

Note that this is the same API as `~sunraster.SpectrogramCube` except that
`sunraster.SpectrogramSequence.dimensions` returns an iterable of
`~astropy.units.Quantity` objects, one for each axis.
This is because of its inheritance from  `~ndcube.NDCubeSequence`
rather than `~ndcube.NDCube`.
Also note that there are now four dimensions, as the sequence is treated
as though it were an additional data axis.
This can be very helpful if you have a series of 2D spectrograms and
want to use the sequence axis to represent time.
`sunraster.SpectrogramSequence.array_axis_physical_types`
returns a list of tuples of the same `IVOA UCD1+ controlled words
<http://www.ivoa.net/documents/REC/UCD/UCDlist-20070402.html>`
used by `sunraster.SpectrogramCube.array_axis_physical_types`.
The sequence axis is given the label ``'meta.obs.sequence'``.

.. _sequence_coords:

Coordinates
^^^^^^^^^^^

Coordinate Properties
*********************

Just like `~sunraster.SpectrogramCube`, `~sunraster.SpectrogramSequence`
provides convenience properties to retrieve the real world coordinate values
for each pixel along each axis, namely
`sunraster.SpectrogramSequence.lon`, `sunraster.SpectrogramSequence.lat`,
`sunraster.SpectrogramSequence.spectral`, `sunraster.SpectrogramSequence.time` and
`sunraster.SpectrogramSequence.exposure_time`.
Since there is no guarantee that `~sunraster.SpectrogramCube`'s WCS transformations
are consistent between `~sunraster.SpectrogramCube` s, `sunraster.SpectrogramCube.lon`
and `sunraster.SpectrogramCube.lat` return 3-D `~astropy.units.Quantity` instances
and `sunraster.SpectrogramCube.spectral` returns a 2-D `~astropy.units.Quantity`
where the additional dimension represent the coordinates for different
`~sunraster.SpectrogramCube` instances.

.. _sequence_slicing:

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

Analogous to `~sunraster.SpectrogramCube`, `~sunraster.SpectrogramSequence`
also provides a `~sunraster.SpectrogramSequence.apply_exposure_time_coorection`
method. This is simply a wrapper around the `~sunraster.SpectrogramCube` version
that saves users from apply or removing the exposure time correction to each
`~sunraster.SpectrogramCube` manually. To remind yourself how that method works,
see the `~sunraster.SpectrogramCube` :ref:`cube_exposure_time_correction` section.

Slicing
^^^^^^^

`~sunraster.SpectrogramSequence` provides an identical slicing API to
`~sunraster.SpectrogramCube`.
Although recall that a `~sunraster.SpectrogramSequence` has an additional dimension.
As with `~sunraster.SpectrogramCube`, the slicing API manipulates not only the
data, but also all relevant supporting metadata including
uncertainties, mask, WCS object, extra_coords, etc.

To slice a `~sunraster.SpectrogramSequence`, simply do:

.. code-block:: python

  >>> my_sequence_roi = my_sequence[1:3, 0:2, 1:3, 1:4]

We can check the effect of the slicing via the
`~sunraster.SpectrogramSequence.dimensions` property.

.. code-block:: python

  >>> print(my_sequence.dimensions)  # Check dimensionality before slicing.
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
  >>> print(my_sequence_roi.dimensions) # See how slicing has changed dimensionality.
  (<Quantity 2. pix>, <Quantity 2. pix>, <Quantity 2. pix>, <Quantity 3. pix>)

Slicing can reduce the dimensionality of `~sunraster.SpectrogramSequence` instances.
For example, let's slice out the 2nd pixel along the slit.

.. code-block:: python

  >>> my_3d_sequence = my_sequence[:, :, 2]
  >>> print(my_3d_sequence.dimensions)
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 5. pix>)

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data,
`~sunraster.SpectrogramSequence` supplies a simple, yet powerful plotting API.
It is intended as a useful quicklook tool and not a replacement
for high quality plots or animations, e.g. for publications or presentations.

.. code-block:: python

  >>> my_sequence.plot() # doctest: +SKIP

As with `~sunraster.SpectrogramCube`, this method produces different types of
visualizations including line plots, 2-D images and 1- and 2-D animations.
Which is displayed depends on the dimensionality of the `~sunraster.SpectrogramSequence`
and the inputs of the user.
For learn more about how to customize plots and animations through the
`~sunraster.SpectrogramSequence.plot` method, see the
`NDCubeSequence plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting>`_.

.. _raster_sequence:

RasterSequence
--------------

Slit spectrographs are often used to produce rasters.
(In fact, it is from this data product that ``sunraster`` derives its name.)
A raster is produced by scanning the slit in discrete steps perpendicular
to its long axis, recording an exposure at each position.
Thus a spectral image over a region is built up over time despite the
slit spectrograph's necessarily narrow horizontal field of view.
Another motivation can be to perform fast repeat raster scans in order to improve
the chances of catching an event with the slit, e.g. a solar flare.
In a raster, the slit-step axis is convolved with time.
Depending on the type of analysis being performed,
users may want to think of their data as if it were in
raster mode/4D (``scan number``, ``slit step``, ``position along slit``, ``wavelength``)
or sit-and-stare mode/3D (``time``, ``position along slit``, ``spectral``).
In order to access the data in the way they want, scientists may often have two
copies, a 3D version and a 4D version.
However, this means scientists have to keep track of two data structures which is
memory intensive both for the scientist and the computer and increases the chances
mistakes in analysis.

Solving this problem is the purpose of the `~sunraster.RasterSequence` class.
It inherits from `~sunraster.SpectrogramSequence` but enables users to label
one of the axes as the slit-step axis.
This in turn facilitates a new set of APIs which allows users to interact with their data
in sit-and-stare (SnS) or rastering mode seemlessly and interchangeably
without having to reformat their data.

Initialization
^^^^^^^^^^^^^^

A `~sunraster.RasterSequence`, is instantiated just like a `~sunraster.SpectrogramCube`.
Let's first create some `~sunraster.SpectrogramCube` instances where each represents a
single raster scan.
As before, we will add the timestamps and exposure times as extra coordinates.

.. code-block:: python

  >>> import numpy as np
  >>> import astropy.wcs
  >>> import astropy.units as u
  >>> from astropy.nddata import StdDevUncertainty
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
  >>> # Define a mask with all pixel unmasked, i.e. mask values = False
  >>> mask = np.zeros(data.shape, dtype=bool)
  >>> # Define some RasterSequence metadata.
  >>> seq_meta = {"description": "This is a RasterSequence."}

  >>> # Define uncertaines for data, 2*data and data/2.
  >>> uncertainties = StdDevUncertainty(np.sqrt(data))
  >>> uncertainties2 = StdDevUncertainty(np.sqrt(data * 2))
  >>> uncertainties05 = StdDevUncertainty(np.sqrt(data * 0.5))

  >>> # Define exposure times.
  >>> exposure_times = np.ones(data.shape[0])/2 * u.s
  >>> axis_length = int(data.shape[0])

  >>> # Create 1st raster
  >>> timestamps0 = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
  ...                     for i in range(axis_length)], format='datetime', scale='utc')
  >>> extra_coords_input0 = [("time", 0, timestamps0), ("exposure time", 0, exposure_times)]
  >>> raster0 = SpectrogramCube(data, input_wcs, uncertainty=uncertainties, mask=mask,
  ...                           meta=meta, unit=u.ct, extra_coords=extra_coords_input0)

  >>> # Create 2nd raster
  >>> timestamps1 = Time([timestamps0[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input1 = [("time", 0, timestamps1), ("exposure time", 0, exposure_times)]
  >>> raster1 = SpectrogramCube(data*2, input_wcs, uncertainty=uncertainties, mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input1)

  >>> # Create 3rd raster
  >>> timestamps2 = Time([timestamps1[-1].to_datetime() + timedelta(minutes=i)
  ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
  >>> extra_coords_input2 = [("time", 0, timestamps2), ("exposure time", 0, exposure_times)]
  >>> raster2 = SpectrogramCube(data*0.5, input_wcs, uncertainty=uncertainties, mask=mask,
  ...                  meta=meta, unit=u.ct, extra_coords=extra_coords_input2)

The last thing we need to do before creating our `~sunraster.RasterSequence` is
to identity the slit-step of the `~sunraster.SpectrogramCube` s.
In the above ``raster`` instances both the 0th and 1st axes correspond to spatial dimensions.
Therefore let's define the 0th axes as the slit-step.
We will do this by setting the ``common_axis`` argument 0.

.. code-block:: python

  >>> from sunraster import RasterSequence
  >>> my_rasters = RasterSequence([raster0, raster1, raster2], common_axis=0, meta=seq_meta)

Dimensions
^^^^^^^^^^

`~sunraster.RasterSequence` provides a version of the
`~sunraster.SpectrogramSequence.array_axis_physical_axis_types` property for
both raster and SnS representations.

.. code-block:: python

  >>> my_rasters.raster_array_axis_physical_types
  [('meta.obs.sequence',),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('em.wl',)]

  >>> my_rasters.SnS_array_axis_physical_types
  [('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'),
   ('em.wl',)]

In the raster case, ``'meta.obs.sequence'`` represents the raster scan number axis.
For those familiar with `~ndcube.NDCubeSequence`, these are simply aliases for the
`~ndcube.NDCubeSequence.array_axis_physical_axis_types` and
`~ndcube.NDCubeSequence.cube_like_world_axis_physical_axis_types`, respectively.

The length of each axis can also be displayed in either the raster or SnS representation.

.. code-block:: python

  >>> my_rasters.raster_dimensions
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)

`~sunraster.RasterSequence.raster_dimensions` always represents the length of
the scan number axis in the 0th position.
We can therefore see that we have 3 raster scans in our `~sunraster.RasterSequence`.
This means that the slit-step axis is shifted by one.
Since we defined ``common_axis=0`` during instantiation, this means that the
length of the slit-step can be found in the 1st element.
From this we can see that we have 3 slit positions per raster scan.

To see the length of the axes as though the data is in sit-and-stare mode, simply do:

.. code-block:: python

  >>> my_rasters.SnS_dimensions
  <Quantity [9., 4., 5.] pix>

Note that scan number and slit-step axes have been combined into the 0th position.
From this we can see that we have 9 (3x3) spectrograms or times in our
`~sunraster.RasterSequence`.

Coordinates
^^^^^^^^^^^

Coordinate Properties
*********************

`~sunraster.RasterSequence` provides the same convenience
properties as `~sunraster.SpectrogramSequence` to retrieve the real world
coordinate values for each pixel along each axis.
`sunraster.RasterSequence.lon`, `sunraster.RasterSequence.lat`,
and `sunraster.RasterSequence.spectral` return their values in the raster
representation while `sunraster.RasterSequence.time` and
`sunraster.RasterSequence.exposure_time` return their values in the SnS
representation.

SnS Axis Extra Coordinates
**************************

As well as `~sunraster.RasterSequence.time` and
`~sunraster.RasterSequence.exposure_time`, some
`sunraster.SpectrogramCube.extra_coords` may contain other coordinates
that are aligned with the slit step axis.
The `sunraster.RasterSequence.SnS_axis_extra_coords` property enables users
to access these coordinates at the `~sunraster.RasterSequence` level in the
form of an abbreviated ``extra_coords`` dictionary.
Just like `~sunraster.RasterSequence.time` and `sunraster.RasterSequence.exposure_time`,
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
with an ``'axis'`` value of ``None``.
Another way of thinking about an ``extra_coord`` with and axis value of ``None``,
is that these coordinates correspond to the raster scan number axis.
Hence the property’s name.

Slicing
^^^^^^^

`~sunraster.RasterSequence` not only enables users to inspect their data in
the raster and sit-and-stare representations.
It also enables them to slice the data in either representation as well.
This is done via the `~sunraster.RasterSequence.slice_as_raster` and
`~sunraster.RasterSequence.slice_as_SnS` properties.
As with `~sunraster.SpectrogramCube` and `~sunraster.SpectrogramSequence`,
these slicing properties ensure that not only the data is sliced,
but also all relevant supporting metadata including uncertainties, mask,
WCS object, extra_coords, etc.

To slice a `~sunraster.RasterSequence` using the raster representation, do:

.. code-block:: python

  >>> my_rasters_roi = my_rasters.slice_as_raster[1:3, 0:2, 1:3, 1:4]

We can see the result of slicing using the ``dimensions`` properties.

.. code-block:: python

  >>> print(my_rasters.raster_dimensions)  # Check dimensionality before slicing.
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
  >>> print(my_rasters_roi.raster_dimensions) # See how slicing has changed dimensionality.
  (<Quantity 2. pix>, <Quantity 2. pix>, <Quantity 2. pix>, <Quantity 3. pix>)
  >>> my_rasters_roi.SnS_dimensions  # Dimensionality can still be represented in SnS form.
  <Quantity [4., 2., 3.] pix>

To slice in the sit-and-stare representation, do the following:

.. code-block:: python

  >>> my_rasters_roi = my_rasters.slice_as_SnS[1:7, 1:3, 1:4]

Let's check the effect of the slicing once again.

.. code-block:: python

  >>> print(my_rasters.SnS_dimensions)  # Check dimensionality before slicing.
  [9. 4. 5.] pix
  >>> print(my_rasters_roi.SnS_dimensions)  # See how slicing has changed dimensionality.
  [6. 2. 3.] pix
  >>> print(my_rasters_roi.raster_dimensions)  # Dimensionality can still be represented in raster form.
  (<Quantity 3. pix>, <Quantity [2., 3., 1.] pix>, <Quantity 2. pix>, <Quantity 3. pix>)

Notice that after slicing the data can still be inspected and interpreted in
the raster or sit-and-stare format, irrespective of which slicing
representation was used.
Also notice that the ``my_sequence.slice_as_SnS[1:7, 1:3, 1:4]`` command led to
different `~sunraster.SpectrogramCube` objects to have different lengths along the
slit step axis.
This can be seen from the fact that the slit step axis entry in the
output of ``my_sequence_roi.raster_dimensions`` has a length greater than 1.
Each element represents the length of each `~sunraster.SpectrogramCube` in the
`~sunraster.SpectrogramSequence` along that axis.

As with `~sunraster.SpectrogramSequence`, slicing can reduce a
`~sunraster.RasterSequence`'s dimensionality.
As in the :ref:`sequence_slicing` section, let's slice out the 2nd pixel along the slit.
This reduces the number of dimensions in the raster representation to 3
(``raster scan``, ``slit step``, ``spectral``) and to 2 in the sit-and-stare
representation (``time``, ``spectral``).
However, the raster and sit-and-stare representations are still valid.

.. code-block:: python

  >>> slit_pixel_rasters = my_rasters.slice_as_raster[:, :, 2]
  >>> print(slit_pixel_rasters.raster_dimensions)
  (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 5. pix>)
  >>> print(slit_pixel_rasters.SnS_dimensions)
  [9. 5.] pix

This demonstrates that the difference between the raster and sit-and-stare
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

  >>> my_rasters.plot_as_raster() # doctest: +SKIP

To visualize in the sit-and-stare representation, do:

.. code-block:: python

  >>> my_rasters.plot_as_SnS() # doctest: +SKIP

These methods produce different types of visualizations including line plots,
2-D images and 1- and 2-D animations.
Which is displayed depends on the dimensionality of the `~sunraster.RasterSequence`
and the inputs of the user.
`~sunraster.RasterSequence.plot_as_raster` and
`~sunraster.RasterSequence.plot_as_SnS` are in fact simply aliases for the
``ndcube.NDCubeSequence.plot`` and ``ndcube.NDCubeSequence.plot_as_cube`` methods,
respectively.
For learn more about how these routines work and the optional inputs that
enable users to customize their output, see the
`NDCubeSequence plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting>`_.

Extracting Data Arrays
^^^^^^^^^^^^^^^^^^^^^^

It is possible that you may have some procedures that are designed to operate
on arrays instead of `~sunraster.SpectrogramSequence` or `~sunraster.RasterSequence`
objects.
Therefore it may be useful to extract the data (or other array-like information
such as `uncertainty` or `mask`) into a single `~numpy.ndarray`.
A succinct way of doing this operation is using python's list comprehension.

To make a 4-D array from the data arrays in ``my_sequence``, use `numpy.stack`.

.. code-block:: python

    >>> print(my_sequence._dimensions)  # Print sequence dimensions as a reminder.
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
    >>> data = np.stack([cube.data for cube in my_sequence.data])
    >>> print(data.shape)
    (3, 3, 4, 5)

To define a 3D array where the data arrays of each `~sunraster.SpectrogramCube`
in the sequence is concatenated along an axis, use `numpy.vstack`.

.. code-block:: python

    >>> data = np.vstack([cube.data for cube in my_sequence.data])
    >>> print(data.shape)
    (9, 4, 5)

To create 3D arrays by slicing sequences, do:

.. code-block:: python

    >>> data = np.stack([cube[2].data for cube in my_sequence.data])
    >>> print(data.shape)
    (3, 4, 5)

Spectrogram Collections
-----------------------

During analysis of slit spectrograph data, it is often desirable to group
different data sets together.
For example, you may have several `~sunraster.SpectrogramCube` or
`~sunraster.RasterSequence` objects representing observations in different
spectral windows.
Or we may have fit a spectral line in each pixel and extracted a property
such as linewidth, thus collapsing the spectral axis.
In both these cases, the `~sunraster.RasterSequence` objects share a common
origin and set of coordinate transformations with the original observations
(except in the spectral axis in the latter example).
However, they do not have a sequential relationship in their common coordinate spaces
and in the latter case the data represents a different physical property to the
original observations.
Therefore, combining them in a `~sunraster.RasterSequence` is not appropriate.

``sunraster`` does not provide a suitable object for this purpose.
However, because `~sunraster.SpectrogramCube` `~sunraster.SpectrogramSequence`
and `~sunraster.RasterSequence` are instances of ``ndcube`` classes underneath,
users can employ the `ndcube.NDCollection` class for this purpose.
`~ndcube.NDCollection` is a ``dict``-like class that provides additional slicing
capabilities of its constituent data cubes along aligned axes.
To see whether `~ndcube.NDCollection` could be helpful for your research, see
the
`NDCollection documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcollection.html>`_.
