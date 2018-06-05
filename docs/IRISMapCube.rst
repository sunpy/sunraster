.. _IRISMapCube:

===========
IRISMapCube
===========

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS).

This class inherits from NDCube_, so we can can apply all the methods derived from this class.

Initialization
--------------

To initialize an IRISMapCube, we will need to open a level 2 fits file. You can find one
of these files in the IRIS_ website by selecting a file either by clicking on the map or by
setting some information in the left panel. Now that we get a level 2 fits file, we can
load it into the reading method for fits files. This method will read all the fits file and
will take all the information that it needs to create an IRISMapCube object.

Let's assume that we will call our fits file ``my_fits_file`` and IRISMapCube object as
``my_cube``: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file)

If you don't have a lot of RAM memory or if you are loading a huge file, we recommend to
use the memmap kwarg. By using it, you will only load what you need to run but some
methods that requires all the file will not be accessible. You can use the memmap
kwarg by doing: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file, memmap=True)

Class Structure
---------------

So now, ``my_cube`` is an IRISMapCube object and we can access to a lot of information like:

- ``my_cube.data`` : In this attribute, we can find the data array.
- ``my_cube.wcs`` : This attribute contains the WCS.
- ``my_cube.mask`` : This attribute can be used to mask some data (eg. the dust particle).
- ``my_cube.uncertainty`` : The data uncertainty is calculated as the square root of data
  plus squared reading noise in photon unit.

Dimensions
----------

As the IRISMapCube object inherits from NDCube, we can use the two properties of NDCube
which allow us to the data shape and the axis types of our IRISMapCube object. ::

  >>> my_cube.dimensions
  <Quantity [3., 4., 5.] pix>
  >>> my_cube.world_axis_physical_types
  ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

Cropping and Indexing
---------------------

One of the most powerful capability of IRISMapCube, coming from NDCube, is the slicing
process. To slice the cube, we can slice the IRISMapCube with an Array-like Indexing or
we can crop it by the Real World Coordinates.

Array-like Indexing
^^^^^^^^^^^^^^^^^^^

As the IRISMapCube object inherits from NDCube, we can use the Array-like Indexing process
described in NDCube.Slicing_ paragraph.

Cropping by Real World Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As the IRISMapCube object inherits from NDCube, we can use the Cropping by Real World
Coordinates process described in NDCube.Slicing_ paragraph.

Data Manipulation
-----------------

Now we have our IRISMapCube object and we know how access to all the information it contains,
we can manipulate the data with the below presented methods.

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

Dust particle mask
^^^^^^^^^^^^^^^^^^

If we take time to look inside the data, we can see that some spots are not relevant with
the data we want to study. These spots are dust particles that can come from space or a
not perfect construction of the satellite. In IRISMapCube, we can use a method that will
modify the mask of our data with the dust particles positions. We can use this method
by doing: ::

    >>> my_cube.apply_dust_mask()

Now, our ``my_cube.mask`` contains the dust particles positions and we can use it to
select only the data we want. If we want to remove the dust particles positions from
our mask, we can call again this method with the ``undo`` kwarg. ::

    >>> my_cube.apply_dust_mask(undo=True)

Visualization
-------------

We can also print a representation of the IRISMapCube object and we can see something
like this: ::

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

Or we can plot the data of our IRISMapCube object in an image or animation by doing that: ::

    >>> my_cube.plot()

We can customize the visualization by using standard matplotlib kwargs relevant to the type of
visualization produces by the plot method. For example, for a 2D image/animation, we can use
``vmin`` and ``vmax`` to set the floor and ceiling of the color map like so: ::

    >>> my_cube.plot(vmin=0, vmax=300)

.. _NDCube: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html
.. _IRIS: http://iris.lmsal.com/search/
.. _NDCube.Dimensions: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#dimensions
.. _NDCube.Slicing: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#slicing
