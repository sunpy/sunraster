===================
IRISMapCubeSequence
===================

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS). This class is a list of IRISMapCube objects, one per time step, from the
same observation.

Initialization
--------------

To initialize an IRISMapCubeSequence, you will need to open some level 2 fits file.
You can find these files in the IRIS_ website by selecting a file either by clicking
on the map or by setting some information in the left panel. Now that you get a level 2
fits file, you can load it into the reading method for fits files. This method will read
all the fits file and will take every information that it needs to create an
IRISMapCubeSequence object.

Let assume that we will build an IRISMapCubeSequence with two fits files that we will
call ``my_fits_file_0`` and ``my_fits_file_1``. You can create an IRISMapCubeSequence
with as much fits files as you want but from the same observation ! Next, we will call
our IRISMapCubeSequence object as ``my_sequence``:

.. code-block:: python

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files)

As for IRISMapCube, if you don't have a lot of RAM memory of if you are loading a huge file,
we recommend to use the memmap kwarg. By using it, you will only load what you need to run
but some methods that requires all the file will not be accessible:

.. code-block:: python

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_cube = read_iris_sji_level2_fits(my_fits_files, memmap=True)

So now, ``my_sequence`` is an IRISMapCubeSequence object and we can access to a lot of
information like:

- ``my_sequence.data`` : In this attribute, we can find the data array of our sequence.

But we also have access to all information stored in the IRISMapCube objects like:

- ``my_sequence[0].data`` : The data attribute of the first IRISMapCube object.
- ``my_sequence[1].wcs`` : The WCS attribute of the second IRISMapCube object.
- ``my_sequence[2].mask`` : The mask attribute of the third IRISMapCube object.

Some of these information can be modified by using one of followed methods.

Exposure Time Correction
------------------------

Waiting the incoming PR

.. _IRIS: http://iris.lmsal.com/search/
