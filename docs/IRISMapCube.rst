===========
IRISMapCube
===========

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS).

Initialization
--------------

To initialize an IRISMapCube, you will need to open a level 2 fits file. You can find one
of these files in the IRIS_ website by selecting a file either by clicking on the map or by
setting some information in the left panel. Now that you get a level 2 fits file, you can
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

Some of these information can be modified by using one of followed methods.

Exposure Time Correction
------------------------

You can apply the exposure time correction to data and to the uncertainty and
this method adjusts the unit:

.. code-block:: python

    >>> my_cube.apply_exposure_time_correction()

You also can undo the exposure time correction by doing:

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
