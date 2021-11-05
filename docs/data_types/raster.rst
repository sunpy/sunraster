.. _raster_sequence:

RasterSequence
--------------

Slit spectrographs are often used to produce rasters.
In fact, it is from this data product that ``sunraster`` derives its name.

A raster is produced by scanning the slit in discrete steps perpendicular to its long axis, recording an exposure at each position.
Thus a spectral image over a region is built up over time despite the slit spectrograph's necessarily narrow horizontal field of view.
Another motivation can be to perform fast repeat raster scans in order to improve the chances of catching an event with the slit, e.g., a solar flare.
In a raster, the slit-step axis is convolved with time.

Depending on the type of analysis being performed, users may want to think of their data as if it were in raster mode/4D (``scan number``, ``slit step``, ``position along slit``, ``wavelength``) or sit-and-stare mode/3D (``time``, ``position along slit``, ``spectral``).


In order to access the data in the way they want, scientists may often have two copies, a 3D version and a 4D version.
However, this means scientists have to keep track of two data structures which is memory intensive both for the scientist and the computer and increases the chances mistakes in analysis.

Solving this problem is the purpose of the `~sunraster.RasterSequence` class.
It inherits from `~sunraster.SpectrogramSequence` but enables users to label one of the axes as the slit-step axis.
This in turn facilitates a new set of APIs which allows users to interact with their data in sit-and-stare (sns) or rastering mode seamlessly and interchangeably without having to reformat their data.

Initialization
^^^^^^^^^^^^^^

A `~sunraster.RasterSequence`, is instantiated just like a `~sunraster.SpectrogramCube`.
Let's first create some `~sunraster.SpectrogramCube` instances where each represents a single raster scan.
As before, we will add the timestamps and exposure times as extra coordinates.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.wcs
    >>> import astropy.units as u
    >>> from astropy.nddata import StdDevUncertainty
    >>> from datetime import datetime, timedelta
    >>> from astropy.time import Time
    >>> from sunraster import SpectrogramCube
    >>> from sunraster.meta import Meta

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
    >>> exposure_times = np.ones(data.shape[0])/2 * u.s
    >>> scan_meta = Meta({"exposure time": exposure_times}, axes={"exposure time": 0},
    ...                  data_shape=data.shape)
    >>> seq_meta = {"description": "This is a RasterSequence.", "exposure time" : exposure_times}

    >>> # Define uncertainties for data, 2*data and data/2.
    >>> uncertainties = StdDevUncertainty(np.sqrt(data))
    >>> uncertainties2 = StdDevUncertainty(np.sqrt(data * 2))
    >>> uncertainties05 = StdDevUncertainty(np.sqrt(data * 0.5))

    >>> # Create 1st raster
    >>> axis_length = int(data.shape[0])
    >>> timestamps0 = Time([datetime(2000, 1, 1) + timedelta(minutes=i)
    ...                     for i in range(axis_length)], format='datetime', scale='utc')
    >>> extra_coords_input0 = [("time", 0, timestamps0)]
    >>> raster0 = SpectrogramCube(data, input_wcs, uncertainty=uncertainties, mask=mask,
    ...                           meta=scan_meta, unit=u.ct)
    >>> for extra in extra_coords_input0:
    ...     raster0.extra_coords.add(*extra)
    >>> # Create 2nd raster
    >>> timestamps1 = Time([timestamps0[-1].to_datetime() + timedelta(minutes=i)
    ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
    >>> extra_coords_input1 = [("time", 0, timestamps1)]
    >>> raster1 = SpectrogramCube(data*2, input_wcs, uncertainty=uncertainties, mask=mask,
    ...                  meta=scan_meta, unit=u.ct)
    >>> for extra in extra_coords_input1:
    ...     raster1.extra_coords.add(*extra)
    >>> # Create 3rd raster
    >>> timestamps2 = Time([timestamps1[-1].to_datetime() + timedelta(minutes=i)
    ...                     for i in range(1, axis_length+1)], format='datetime', scale='utc')
    >>> extra_coords_input2 = [("time", 0, timestamps2)]
    >>> raster2 = SpectrogramCube(data*0.5, input_wcs, uncertainty=uncertainties, mask=mask,
    ...                  meta=scan_meta, unit=u.ct)
    >>> for extra in extra_coords_input2:
    ...     raster2.extra_coords.add(*extra)

The last thing we need to do before creating our `~sunraster.RasterSequence` is to identity the slit-step of the `~sunraster.SpectrogramCube`.
In the above ``raster`` instances both the 0th and 1st axes correspond to spatial dimensions.
Therefore let's define the 0th axes as the slit-step.
We will do this by setting the ``common_axis`` argument 0.

.. code-block:: python

    >>> from sunraster import RasterSequence
    >>> my_rasters = RasterSequence([raster0, raster1, raster2], common_axis=0, meta=seq_meta)

Dimensions
^^^^^^^^^^

`~sunraster.RasterSequence` provides a version of the `~sunraster.SpectrogramSequence.array_axis_physical_axis_types` property for both raster and sns representations.

.. code-block:: python

    >>> my_rasters.raster_array_axis_physical_types
    [('meta.obs.sequence',), ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon', 'time'), ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'), ('em.wl',)]

    >>> my_rasters.sns_array_axis_physical_types
    [('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon', 'time'), ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon'), ('em.wl',)]

In the raster case, ``'meta.obs.sequence'`` represents the raster scan number axis.
For those familiar with `~ndcube.NDCubeSequence`, these are simply aliases for the `~ndcube.NDCubeSequence.array_axis_physical_axis_types` and `~ndcube.NDCubeSequence.cube_like_world_axis_physical_axis_types`, respectively.

The length of each axis can also be displayed in either the raster or sns representation.

.. code-block:: python

    >>> my_rasters.raster_dimensions
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)

`~sunraster.RasterSequence.raster_dimensions` always represents the length of the scan number axis in the 0th position.
We can therefore see that we have 3 raster scans in our `~sunraster.RasterSequence`.
This means that the slit-step axis is shifted by one.
Since we defined ``common_axis=0`` during instantiation, this means that the length of the slit-step can be found in the 1st element.
From this we can see that we have 3 slit positions per raster scan.

To see the length of the axes as though the data is in sit-and-stare mode, simply do:

.. code-block:: python

    >>> my_rasters.sns_dimensions
    <Quantity [9., 4., 5.] pix>

Note that scan number and slit-step axes have been combined into the 0th position.
From this we can see that we have 9 (3x3) spectrograms or times in our `~sunraster.RasterSequence`.

Coordinates
^^^^^^^^^^^

Coordinate properties
*********************

`~sunraster.RasterSequence` provides the same convenience properties as `~sunraster.SpectrogramSequence` to retrieve the real world coordinate values for each pixel along each axis.
`sunraster.RasterSequence.celestial`, and `sunraster.RasterSequence.spectral` return their values in the raster representation while `sunraster.RasterSequence.time` and `sunraster.RasterSequence.exposure_time` return their values in the sns representation.

sns axis extra coordinates
**************************

As well as `~sunraster.RasterSequence.time` and `~sunraster.RasterSequence.exposure_time`, some `sunraster.SpectrogramCube.extra_coords` may contain other coordinates that are aligned with the slit step axis.
The `sunraster.RasterSequence.sns_axis_coords` property enables users to access these coordinates at the `~sunraster.RasterSequence` level in the form of an abbreviated ``extra_coords`` dictionary.
Just like `~sunraster.RasterSequence.time` and `sunraster.RasterSequence.exposure_time`, the coordinates are concatenated so they mimic the sit-and-stare-like dimensionality returned in the 0th element of `sunraster.RasterSequence.sns_dimensions`.
`sunraster.RasterSequence.sns_axis_coords` is equivalent to `ndcube.NDCubeSequence.common_axis_extra_coords`.
To see examples of how to use this property, see the `NDCubeSequence Common Axis Extra Coordinates documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#common-axis-extra-coordinates>`__.

Raster axis extra coordinates
*****************************

Analogous to `~sunraster.RasterSequence.sns_axis_coords`, it is also possible to access the coordinates that are not assigned to any `~sunraster.SpectrogramCube` data axis via the `~sunraster.RasterSequence.raster_axis_coords` property.
This property is equivalent to `ndcube.NDCubeSequence.sequence_axis_coords` and can be used to return coordinates along the repeat raster scan axis.

Slicing
^^^^^^^

`~sunraster.RasterSequence` not only enables users to inspect their data in the raster and sit-and-stare representations.
It also enables them to slice the data in either representation as well.
This is done via the `~sunraster.RasterSequence.slice_as_raster` and `~sunraster.RasterSequence.slice_as_sns` properties.
As with `~sunraster.SpectrogramCube` and `~sunraster.SpectrogramSequence`, these slicing properties ensure that not only the data is sliced, but also all relevant supporting metadata including uncertainties, mask, WCS object, extra_coords, etc.

To slice a `~sunraster.RasterSequence` using the raster representation, do:

.. code-block:: python

    >>> my_rasters_roi = my_rasters.slice_as_raster[1:3, 0:2, 1:3, 1:4]

We can see the result of slicing using the ``dimensions`` properties.

.. code-block:: python

    >>> print(my_rasters.raster_dimensions)  # Check dimensionality before slicing.
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
    >>> print(my_rasters_roi.raster_dimensions) # See how slicing has changed dimensionality.
    (<Quantity 2. pix>, <Quantity 2. pix>, <Quantity 2. pix>, <Quantity 3. pix>)
    >>> my_rasters_roi.sns_dimensions  # Dimensionality can still be represented in sns form.
    <Quantity [4., 2., 3.] pix>

To slice in the sit-and-stare representation, do the following:

.. code-block:: python

    >>> my_rasters_roi = my_rasters.slice_as_sns[1:7, 1:3, 1:4]

Let's check the effect of the slicing once again.

.. code-block:: python

    >>> print(my_rasters.sns_dimensions)  # Check dimensionality before slicing.
    [9. 4. 5.] pix
    >>> print(my_rasters_roi.sns_dimensions)  # See how slicing has changed dimensionality.
    [6. 2. 3.] pix
    >>> print(my_rasters_roi.raster_dimensions)  # Dimensionality can still be represented in raster form.
    (<Quantity 3. pix>, <Quantity [2., 3., 1.] pix>, <Quantity 2. pix>, <Quantity 3. pix>)

Notice that after slicing the data can still be inspected and interpreted in the raster or sit-and-stare format, irrespective of which slicing representation was used.
Also notice that the ``my_sequence.slice_as_sns[1:7, 1:3, 1:4]`` command led to different `~sunraster.SpectrogramCube` objects to have different lengths along the slit step axis.
This can be seen from the fact that the slit step axis entry in the output of ``my_sequence_roi.raster_dimensions`` has a length greater than 1.
Each element represents the length of each `~sunraster.SpectrogramCube` in the `~sunraster.SpectrogramSequence` along that axis.

As with `~sunraster.SpectrogramSequence`, slicing can reduce a `~sunraster.RasterSequence` dimensionality.
As in the :ref:`sequence_slicing` section, let's slice out the 2nd pixel along the slit.
This reduces the number of dimensions in the raster representation to 3 (``raster scan``, ``slit step``, ``spectral``) and to 2 in the sit-and-stare representation (``time``, ``spectral``).
However, the raster and sit-and-stare representations are still valid.

.. code-block:: python

    >>> slit_pixel_rasters = my_rasters.slice_as_raster[:, :, 2]
    >>> print(slit_pixel_rasters.raster_dimensions)
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 5. pix>)
    >>> print(slit_pixel_rasters.sns_dimensions)
    [9. 5.] pix

This demonstrates that the difference between the raster and sit-and-stare representations is more subtle than simply a 4-D or 3-D dimensionality.
The difference is whether the raster scan and slit step axes are convolved into a time axis or whether they are represented separately.
And because of this definition, the raster and sit-and-stare representations are valid and accessible for any dimensionality in which the raster scan and slit step axes are maintained.

Plotting
^^^^^^^^

To quickly and easily visualize slit spectrograph data, `~sunraster.RasterSequence` supplies simple-to-use, yet powerful plotting APIs.
They are intended to be a useful quicklook tool and not a replacement for high quality plots or animations, e.g. for publications.
As with slicing, there are two plot methods for plotting in each of the raster and sit-and-stare representations.

To visualize in the raster representation, simply call the following:

.. code-block:: python

    >>> my_rasters.plot_as_raster() # doctest: +SKIP

To visualize in the sit-and-stare representation, do:

.. code-block:: python

    >>> my_rasters.plot_as_sns() # doctest: +SKIP

These methods produce different types of visualizations including line plots, 2-D images and 1- and 2-D animations.
Which is displayed depends on the dimensionality of the `~sunraster.RasterSequence` and the inputs of the user.
`~sunraster.RasterSequence.plot_as_raster` and `~sunraster.RasterSequence.plot_as_sns` are in fact simply aliases for the ``ndcube.NDCubeSequence.plot`` and ``ndcube.NDCubeSequence.plot_as_cube`` methods, respectively.
For learn more about how these routines work and the optional inputs that enable users to customize their output, see the `NDCubeSequence plotting documentation <https://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting>`__.

Extracting Data Arrays
^^^^^^^^^^^^^^^^^^^^^^

It is possible that you may have some procedures that are designed to operate on arrays instead of `~sunraster.SpectrogramSequence` or `~sunraster.RasterSequence` objects.
Therefore it may be useful to extract the data (or other array-like information such as `uncertainty` or `mask`) into a single `~numpy.ndarray`.
A succinct way of doing this operation is using python's list comprehension.

To make a 4-D array from the data arrays in ``my_rasters``, use `numpy.stack`.

.. code-block:: python

    >>> print(my_rasters._dimensions)  # Print sequence dimensions as a reminder.
    (<Quantity 3. pix>, <Quantity 3. pix>, <Quantity 4. pix>, <Quantity 5. pix>)
    >>> data = np.stack([cube.data for cube in my_rasters.data])
    >>> print(data.shape)
    (3, 3, 4, 5)

To define a 3D array where the data arrays of each `~sunraster.SpectrogramCube`
in the sequence is concatenated along an axis, use `numpy.vstack`.

.. code-block:: python

    >>> data = np.vstack([cube.data for cube in my_rasters.data])
    >>> print(data.shape)
    (9, 4, 5)

To create 3D arrays by slicing sequences, do:

.. code-block:: python

    >>> data = np.stack([cube[2].data for cube in my_rasters.data])
    >>> print(data.shape)
    (3, 4, 5)
