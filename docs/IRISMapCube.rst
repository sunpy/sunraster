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

As the IRISMapCube object inherits from NDCube_, we can use the two properties of NDCube_
which allow us to the data shape and the axis types of our IRISMapCube object. These
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
modify the mask of our object with the dust particles positions. We can use this method
by doing: ::

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

There is two different ways to visualize our IRISMapCube object. The first one is to use
the representation of our object and the second one is to plot the data of our object.

Representation
^^^^^^^^^^^^^^

In this part, we can have a look of our IRISMapCube object by using its representation
property. ::

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
how many frames its contains, the dimensions of it, etc ...

Plotting
^^^^^^^^

As the IRISMapCube object inherits from NDCube_, we can use the plotting method of NDCube_
which allow us to see the data in plots or animations. This method is described in the
NDCube.Plotting_ section.

.. _NDCube: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html
.. _IRIS: http://iris.lmsal.com/search/
.. _NDCube.Dimensions: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#dimensions
.. _NDCube.Slicing: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#slicing
.. _NDCube.Plotting: http://docs.sunpy.org/projects/ndcube/en/stable/ndcube.html#plotting
