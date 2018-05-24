===================
IRISMapCubeSequence
===================

This class represents Slit Jaw Images from IRIS described by a single World Coordinate
System (WCS). This class is a list of IRISMapCube objects, one per time step, from the
same observation.

Initialization
--------------

To initialize an IRISMapCubeSequence, we will need to open some level 2 fits file.
You can find these files in the IRIS_ website by selecting a file either by clicking
on the map or by setting some information in the left panel. Now that we get some level 2
fits file, we can load it into the reading method for fits files. This method will read
all the fits file and will take every information that it needs to create an
IRISMapCubeSequence object.

Let assume that we will build an IRISMapCubeSequence with two fits files that we will
call ``my_fits_file_0`` and ``my_fits_file_1``. We can create an IRISMapCubeSequence
with as much fits files as we want but from the same observation ! Next, we will call
our IRISMapCubeSequence object as ``my_sequence``: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files)

As for IRISMapCube, if you don't have a lot of RAM memory of if you are loading a huge file,
we recommend to use the memmap kwarg. By using it, you will only load what you need to run
but some methods that requires all the file will not be accessible: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_cube = read_iris_sji_level2_fits(my_fits_files, memmap=True)

So now, ``my_sequence`` is an IRISMapCubeSequence object and we can access to some
information like:

- ``my_sequence.data`` : In this attribute, we can find the data array of our sequence.

As an IRISMapCubeSequence object is a list of IRISMapCube objects, we can do the same things
that we have seen previously in the IRISMapCube documentation. To access to the information
stored in a IRISMapCube object, we can select the first IRISMapCube object by writting
``my_sequence[0]`` or the Nth object with ``my_sequence[N-1]``. As they are IRISMapCube
objects, we can access to their information by doing:

- ``my_sequence[0].data`` : The data attribute of the first IRISMapCube object.
- ``my_sequence[1].wcs`` : The WCS attribute of the second IRISMapCube object.
- ``my_sequence[2].mask`` : The mask attribute of the third IRISMapCube object.

We can also print a representation of the IRISMapCubeSequence object and we can see
something like this: ::

    >>> my_sequence
    IRISMapCubeSequence
    ---------------------
    Observatory:		IRIS
    Instrument:		 SJI

    OBS ID:			3690015104
    OBS Description:	 Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x
    OBS period:		 2018-04-26T23:07:22.780000 -- 2018-04-27T01:39:47.122000

    Sequence period:	 2018-04-26T23:07:22.880000 -- 2018-04-27T01:36:40.490000
    Sequence Shape:		[  98. 1095. 1018.] pix
    Axis Types:		 ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

Or we can print the data of our IRISMapCubeSequence object in an image or animation
by doing that: ::

    >>> my_sequence.plot()

We can upgrade the visualization by using kwargs to set value limits. Here, ``vmin`` is
the minimum value and every lower value will be considered as the ``vmin`` value.
This is the same for ``vmax`` that corresponds to the maximum value: ::

    >>> my_sequence.plot(vmin=0, vmax=300)

Now that we have created our IRISMapCubeSequence object, we can use one of the followed methods
to correct the data.

Exposure Time Correction
------------------------

Waiting the incoming PR

.. _IRIS: http://iris.lmsal.com/search/
