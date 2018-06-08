.. _IRISMapCube:

===========
IRISMapCube
===========

This class provides access to IRIS Slit-jaw image (SJI) files.

This class inherits from NDCube_, so many methods of IRISMapCube are linked to their
description in NDCube_ documentation.

Initialisation
--------------

To initialize an IRISMapCube, we need to open an IRIS level 2 SJI FITS file. Such data can be downloaded from the IRIS_. The level 2 FITS files can be loaded with the function ``read_iris_sji_level2_fits``. This function loads a FITS file into an IRISMapCube object.

Assuming we have a FITS file ``my_fits_file``, an IRISMapCube can be initialised in the following way: ::

    >>> from irispy.sji import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file)

If you don't have a lot of RAM memory or if you are loading a huge file, we recommend to
use the ``memmap=True`` keyword. By using it, you will only load the data when needed. However,
some methods that require all data in memory will not be accessible. You can use memmap
by doing: ::

    >>> from irispy.sji import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file, memmap=True)

Once we have our ``IRISMapCube`` object, we can quickly inspect it. For that, we can use the ``__repr__`` property of our object, which shows some metadata. This property is called just by writing the name of our
object in the command line: ::

    >>> my_cube
    IRISMapCube
    -----------
    Observatory:		 IRIS
    Instrument:			 SJI
    Bandpass:			 1330.0
    Obs. Start:			 2018-04-26T23:07:22.780000
    Obs. End:			 2018-04-27T01:39:47.122000
    Instance Start:		 2018-04-26T23:07:22.880000
    Instance End:		 2018-04-27T01:36:40.490000
    Total Frames in Obs.:	 49
    IRIS Obs. id:		 3690015104
    IRIS Obs. Description:	 Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x
    Cube dimensions:		 [  49. 1095. 1018.] pix
    Axis Types:			 ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

Here, we find some information about our object, such as when the observation was made,
how many frames it contains, its dimensions, etc ...

Class Structure
---------------

Our IRISMapCube object, here called ``my_cube``, is now created and inspected, we look into its structure.

Data
^^^^

The ``data`` attribute provides direct access to the data arrays. These values have a unit
(see :ref:`Unit`), an uncertainty (see :ref:`Uncertainty`), and can be scaled or not
(see :ref:`Scaled_Value`). A quick representation of these data can be seen by entering: ::

    >>>my_cube.data

This will show a summary start and and of the data array. It such printouts it is common to see many ``nan`` values. This is because IRIS SJI images are often padded, and these values are masked (see :ref:`Mask`).

World Coordinates System
^^^^^^^^^^^^^^^^^^^^^^^^

The ``wcs`` attribute stores coordinate information (space and time) of our object. We can do a quick inspection of these information by entering: ::

    >>> my_cube.wcs
    WCS Keywords

    Number of WCS axes: 3
    CTYPE : 'HPLN-TAN'  'HPLT-TAN'  'Time'
    CRVAL : 0.11982805555555555  0.03607138888888889  4478.91
    CRPIX : 509.5  548.0  25.0
    PC1_1 PC1_2 PC1_3  : 0.99993622303  0.0112722599879  0.0
    PC2_1 PC2_2 PC2_3  : -0.0112722599879  0.99993622303  0.0
    PC3_1 PC3_2 PC3_3  : 0.0  0.0  1.0
    CDELT : 4.620833333333333e-05  4.620833333333333e-05  186.617
    NAXIS : 1018  1095  49

Here, we can find that our object has 3 axes: a helioprojective longitude, a helioprojective latitude and time. The dimension order is from the FITS file (written in Fortran order), and will be reversed when handling the arrays in Python (which uses C order). Therefore, time is usually the first dimension in the ``data`` array.

The ``IRISMapCube`` object is aware of the WCS matrix as read from the FITS file: ``PC`` values represent the rotation matrix, ``CRVAL`` values represent coordinate values at the pixel positions of ``CRPIX``, and ``CDELT`` represent the pixel size in coordinate units. ``NAXIS`` represents the maximum dimension of each each axis.

.. _Uncertainty:

Uncertainty
^^^^^^^^^^^

We can also found another array inside our object, stored in the ``uncertainty`` attribute.
The uncertainty is calculated as the square root of our object data plus squared reading
noise in photon unit. We can inspect the array by entering: ::

    >>> my_cube.uncertainty

This will return a summary of the uncertainty values. The ``nan`` values are
also coming from the mask (see :ref:`Mask`).

.. _Unit:

Units
^^^^^

Inside the ``unit`` attribute, we can find the data units, typically in data number (DN). These can be changed by applying some methods (e.g. :ref:`Exposure_Time_Correction` method). We can inspect the units by entering: ::

    >>> my_cube.unit
    Unit("DN_IRIS_SJI")

By default, the units are ``Unit("DN_IRIS_SJI")``, which is calculated by dividing the
detector gain by the detector yield in photon units.

Meta
^^^^

The ``meta`` attribute is storing a dictionary with some metadata about our object. We can inspect it by entering: ::

    >>> my_cube.meta
    {'ENDOBS': datetime.datetime(2018, 4, 27, 1, 39, 47, 122000),
     'INSTRUME': 'SJI',
     'NBFRAMES': 49,
     'OBSID': '3690015104',
     'OBS_DESC': 'Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x',
     'STARTOBS': datetime.datetime(2018, 4, 26, 23, 7, 22, 780000),
     'TELESCOP': 'IRIS',
     'TWAVE1': 1330.0}

And we can also select only one key (eg. ``OBSID``) with the line: ::

    >>> my_cube.meta['OBSID']
    '3690015104'

.. _Mask:

Mask
^^^^

The mask attribute is a boolean array with the same shape as the data. When ``True`` (masked values), it represents regions with invalid or missing data. For example,
we can use it to mask the dust particle positions on the data by using the
:ref:`Dust_Particle_Mask` method. We can inspect the mask by entering: ::

    >>> my_cube.mask

By default, the mask is set to include the unexposed pixels of the detector.

Extra Coordinates
^^^^^^^^^^^^^^^^^

As our ``IRISMapCube`` object inherits from NDCube_, this attribute is explained in the
NDCube.Extra_Coordinates_ section in the NDCube_ documentation. We can access this
dictionary with: ::

    >>> my_cube.extra_coords

To select only one key (eg. ``TIME``), we can do: ::

    >>> my_cube.extra_coords['TIME']
    {'axis': 0, 'value': array([datetime.datetime(2018, 4, 26, 23, 7, 22, 880000), ... ], dtype=object)}

We can see that this is an other dictionary, so we can select the first value of the
``TIME`` by doing: ::

    >>> my_cube.extra_coords['TIME']['value'][0]
    datetime.datetime(2018, 4, 26, 23, 7, 22, 880000)

Missing axes
^^^^^^^^^^^^

This ``NDCube`` attribute is explained in the NDCube.Missing_Axes_ section in the
NDCube_ documentation. We can inspect it by entering: ::

    >>> my_cube.missing_axis
    [False, False, False]

.. _Scaled_Value:

Scaled values
^^^^^^^^^^^^^

This boolean attribute is used to check if the data values are scaled or not. Scaling happens when the data are read from the FITS file and converted to data number (DN), and unscaled data are read directly from the FITS file without any conversion (from 16-bit integer to 32-bit float). The default value is ``True``, and it is ``False`` when using memmap (see above) during the creation of the object. We can inspect it by entering: ::

    >>> my_cube.scaled
    True

Dimensions
----------

As ``IRISMapCube`` is inherited from NDCube_, we can use the two properties of NDCube_
which allow us to get the data shape and the axis types of our ``IRISMapCube`` object. These
properties are described in the NDCube.Dimensions_ section.

Cropping and Indexing
---------------------

One of the most powerful capabilities of ``IRISMapCube``, coming from NDCube_, is the slicing
ability. There are two ways to slice: using array-like index or by coordinates. These are described in the NDCube.Slicing_ section.

Manipulation the Data
---------------------

We can manipulate an ``IRISMapCube`` object with the methods listed below.

.. _Exposure_Time_Correction:

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

We can apply the exposure time correction to data and to the uncertainty and
this method adjusts the unit: ::

    >>> my_cube.apply_exposure_time_correction()

We also can undo the exposure time correction by doing: ::

    >>> my_cube.apply_exposure_time_correction(undo=True)

The correction is only applied (undone) if the object's unit doesn't (does) already
include inverse time. This can be overridden so that correction is applied (undone)
regardless of unit by setting ``force=True``. Use one of the two lines above to apply
(undone) by using the force kwarg: ::

    >>> my_cube.apply_exposure_time_correction(force=True)
    >>> my_cube.apply_exposure_time_correction(undo=True, force=True)

.. _Dust_Particle_Mask:

Dust particle mask
^^^^^^^^^^^^^^^^^^

If we take time to look inside the data, we see that some pixels are obscured by dust,
and so do not reflect the emission from the Sun at that location. The ``apply_dust_mask``
method of IRISMapCube can be used to add the locations of the dust pixels to the mask
so that we can easily exclude them from our analysis. We can use this method by doing: ::

    >>> my_cube.apply_dust_mask()

Now, our ``my_cube.mask`` contains the dust particles positions and we can use it to
select only the data we want. If we want to remove the dust particles positions from
our mask, we can call again this method with the ``undo`` kwarg. ::

    >>> my_cube.apply_dust_mask(undo=True)

If we don't remember or we don't know if the dust particles positions are applied or not
in our ``my_cube.mask``, we can check an attribute of our object. ::

    >>> my_cube.dust_masked

If the result is ``True``, the dust particles positions are applied in our ``my_cube.mask``.
If the result is ``False``, the dust particles positions are not applied.

Visualization
-------------

As the IRISMapCube object inherits from NDCube_, we can use the plotting method of NDCube_
which allow us to see the data in plots or animations. This method is described in the
NDCube.Plotting_ section.

.. _NDCube: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html
.. _IRIS: http://iris.lmsal.com/search/
.. _NDCube.Extra_Coordinates: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#extra-coordinates
.. _NDCube.Missing_Axes: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#missing-axes
.. _NDCube.Dimensions: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#dimensions
.. _NDCube.Slicing: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#slicing
.. _NDCube.Plotting: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#plotting
