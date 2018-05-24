===========
IRISMapCube
===========

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS).

Initialization
--------------

To initialize an IRISMapCube, we will need to open a level 2 fits file. You can find one
of these files in the IRIS_ website by selecting a file either by clicking on the map or by
setting some information in the left panel. Now that we get a level 2 fits file, we can
load it into the reading method for fits files. This method will read all the fits file and
will take every information that it needs to create an IRISMapCube object.

Let assume that we will call our fits file ``my_fits_file`` and IRISMapCube object as
``my_cube``:

.. code-block:: python

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file)

If you don't have a lot of RAM memory of if you are loading a huge file, we recommend to
use the memmap kwarg. By using it, you will only load what you need to run but some
methods that requires all the file will not be accessible:

.. code-block:: python

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_cube = read_iris_sji_level2_fits(my_fits_file, memmap=True)

So now, ``my_cube`` is an IRISMapCube object and we can access to a lot of information like:

- ``my_cube.data`` : In this attribute, we can find the data array.
- ``my_cube.wcs`` : This attribute contains the WCS.
- ``my_cube.mask`` : This attribute can be used to mask some data (eg. the dust particle).

We can also print a representation of the IRISMapCube object and we can see something like this:

.. code-block:: python

    >>> my_cube
    IRISMapCube
    ---------
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

Or we can print the data of our IRISMapCube object in an image or animation by doing that:

.. code-block:: python

    >>> my_cube.plot()

We can upgrade the visualization by using kwargs to set value limits. Here, ``vmin`` is
the minimum value and every lower value will be considered as the ``vmin`` value.
This is the same for ``vmax`` that corresponds to the maximum value:

   .. code-block:: python

    >>> my_cube.plot(vmin=0, vmax=300)

Now that we have created our IRISMapCube object, we can use one of the followed methods
to correct the data.

Exposure Time Correction
------------------------

We can apply the exposure time correction to data and to the uncertainty and
this method adjusts the unit:

.. code-block:: python

    >>> my_cube.apply_exposure_time_correction()

We also can undo the exposure time correction by doing:

.. code-block:: python

    >>> my_cube.apply_exposure_time_correction(undo=True)

The correction is only applied (undone) if the object's unit doesn't (does) already
include inverse time. This can be overridden so that correction is applied (undone)
regardless of unit by setting ``force=True``. Use one of the two lines above to apply
(undone) by using the force kwarg:

.. code-block:: python

    >>> my_cube.apply_exposure_time_correction(force=True)
    >>> my_cube.apply_exposure_time_correction(undo=True, force=True)

Dust particle mask
------------------

Waiting the incoming PR

.. _IRIS: http://iris.lmsal.com/search/
