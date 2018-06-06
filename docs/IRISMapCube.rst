.. _IRISMapCube:

===========
IRISMapCube
===========

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS).

This class inherits from NDCube_, so many methods of IRISMapCube are linked to their
description in NDCube_ documentation.

Initialization
--------------

To initialize an IRISMapCube, we will need to open a level 2 fits file. You can find one
of these files in the IRIS_ website by selecting a file either by clicking on the map or by
setting some information in the left panel. Now that we get a level 2 fits file, we can
load it into the reading method for fits files. This method will read all the fits file and
will take all the information that it needs to create an IRISMapCube object.

Let's assume that we will call our fits file ``my_fits_file`` and IRISMapCube object as
``my_cube``: ::

    >>> from irispy.sji import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file)

If you don't have a lot of RAM memory or if you are loading a huge file, we recommend to
use the memmap kwarg. By using it, you will only load what you need to run but some
methods that requires all the file will not be accessible. You can use the memmap
kwarg by doing: ::

    >>> from irispy.sji import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file, memmap=True)

Now, that we have create our IRISMapCube object, we can check if the creation has been
successfully done. We can use the representation property of our object, where are stored
many information, to do that. This property is called just by writing the name of our
object in the console. ::

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

Here, we can found some information about our object like when the observation was made,
how many frames its contains, its dimensions, etc ...

Class Structure
---------------

Our IRISMapCube object, called ``my_cube`` is now created and checked, we can be interested
in the structure of our object.

Data
^^^^

The ``data`` attribute is storing the value of the Sun emission. These values have a unit
(see :ref:`Unit`), an uncertainty (see :ref:`Uncertainty`), and can be scaled or not
(see :ref:`Scaled_Value`). This attribute is an array that can be visualize by doing: ::

    >>>my_cube.data

The return of this should be an array where we can find a lot of ``nan`` values. This is due
to the mask of our object (see :ref:`Mask`)

World Coordinates System
^^^^^^^^^^^^^^^^^^^^^^^^

This attribute stores information about the axes of our object. We can have a visualization
of these information by doing: ::

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

Here, we can find that our object has 3 axes which are a longitude, a latitude and
the time. We can see that the order of these axes is the inverse of the indexing order.
In fact, we can access to the time axis with the first index of slicing.

The ``PC`` values are representing a rotation matrix that we can convert to a tilt angle
of our object axes.

The ``CRVAL`` values are representing the center of the data in each axis in data unit for
the first frame while ``CRPIX`` represents the same thing in pixel unit.

The ``CDELT`` values are the values that are added per frame in each axis. For example,
the data center in the time axis in the second frame will be ``4478.91 + 186.617 = 4665.527``
in the data unit.

Then, ``NAXIS`` represents the maximum value of the data in each axis in pixel unit.

.. _Uncertainty:

Uncertainty
^^^^^^^^^^^

We can also found an other array inside our object, stored in the ``uncertainty`` attribute.
The uncertainty is calculated as the square root of our object data plus squared reading
noise in photon unit. We can see the array by doing: ::

    >>> my_cube.uncertainty

The return of this line is an array with the same shape than data. The ``nan`` values are
also coming from the mask (see :ref:`Mask`).

.. _Unit:

Unit
^^^^

Inside the ``unit`` attribute, we can found which unit is set to our data. The unit can
change by using some methods (like the :ref:`Exposure_Time_Correction` method). The unit
can be displayed with this line: ::

    >>> my_cube.unit
    Unit("DN_IRIS_SJI")

By default, the unit is ``Unit("DN_IRIS_SJI")`` which is calculated by dividing the
detector gain by the detector yield in photon unit.

Meta
^^^^

The ``meta`` attribute is storing a dictionary with some information used by the
representation property of our object. We can see this dictionary by doing: ::

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

The mask attribute is also an array with the same shape than data and uncertainty arrays.
This array stores Boolean values than we can use to select a part of our data. For example,
we can use it to mask the dust particle positions on the data by using the
:ref:`Dust_Particle_Mask` method. We can access to the mask with this line: ::

    >>> my_cube.mask

By default, the mask is set to remove the unexposed pixels of the detector. This mask is
used when the object is created, this is why all other arrays can store ``nan`` values,
corresponding to the unexposed pixels.

Extra Coordinates
^^^^^^^^^^^^^^^^^

As our IRISMapCube object inherits from NDCube_, this attribute is explained in the
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

As previously, this attribute is explain  in the NDCube.Missing_Axes_ section in the
NDCube_ documentation. We can see this array with this line: ::

    >>> my_cube.missing_axis
    [False, False, False]

.. _Scaled_Value:

Scaled values
^^^^^^^^^^^^^

This attribute is storing a Boolean that remind us if the data values are scaled or not.
The default value of this attribute is ``True`` but the value can be ``False`` if the
memmap kwarg has been set during the creation of the object. We can check if we are using
the memmap kwarg by doing: ::

    >>> my_cube.scaled
    True

Dimensions
----------

As the IRISMapCube object inherits from NDCube_, we can use the two properties of NDCube_
which allow us to get the data shape and the axis types of our IRISMapCube object. These
properties are described in the NDCube.Dimensions_ section.

Cropping and Indexing
---------------------

One of the most powerful capability of IRISMapCube, coming from NDCube_, is the slicing
process. To slice the cube, we can slice the IRISMapCube with an Array-like Indexing or
we can crop it by the Real World Coordinates. As the IRISMapCube object inherits from
NDCube_, we can use the described processes in the NDCube.Slicing_ section.

Data Manipulation
-----------------

Now we have our IRISMapCube object and we know how access to all the information it contains,
we can manipulate the data with the below presented methods.

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
