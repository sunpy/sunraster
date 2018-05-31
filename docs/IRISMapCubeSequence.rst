===================
IRISMapCubeSequence
===================

This class represents Slit Jaw Images from IRIS as a list of IRISMapCube, where these
sub-cubes are 2D or 3D described by their own World Coordinate System (WCS).

This class inherits from the NDCubeSequence_, so we can can apply all the methods derived
from this class.

Initialization
--------------

To initialize an IRISMapCubeSequence, we will need to open some level 2 fits file.
You can find these files in the IRIS_ website by selecting a file either by clicking
on the map or by setting some information in the left panel. Now that we get some level 2
fits file, we can load it into the reading method for fits files. This method will read
all the fits file and will take all the information that it needs to create an
IRISMapCubeSequence object.

Let's assume that we will build an IRISMapCubeSequence with two fits files that we will
call ``my_fits_file_0`` and ``my_fits_file_1``. We can create an IRISMapCubeSequence
with as many fits files as we want but from the same observation. Next, we will call
our IRISMapCubeSequence object as ``my_sequence``. ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files)

As for IRISMapCube, if you don't have a lot of RAM memory or if you are loading a huge file,
we recommend to use the memmap kwarg. By using it, you will only load what you need to run
but some methods that requires all the file will not be accessible. You can use the memmap
kwarg by doing: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files, memmap=True)

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
- ``my_sequence[3].uncertainty`` : The uncertainty attribute of the fourth IRISMapCube object.

We can also print a representation of the IRISMapCubeSequence object and we can see
something like this: ::

    >>> my_sequence
    IRISMapCubeSequence
    -------------------
    Observatory:	 IRIS
    Instrument:		 SJI

    OBS ID:		 3690015104
    OBS Description:	 Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x
    OBS period:		 2018-04-26T23:07:22.780000 -- 2018-04-27T01:39:47.122000

    Sequence period:	 2018-04-26T23:07:22.880000 -- 2018-04-27T01:36:40.490000
    Sequence Shape:	 [  98. 1095. 1018.] pix
    Axis Types:		 ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

Or we can plot the data of our IRISMapCubeSequence object in an image or animation
by doing that: ::

    >>> my_sequence.plot()

We can customize the visualization by using standard matplotlib kwargs relevant to the type of
visualization produces by the plot method. For example, for a 2D image/animation, we can use
``vmin`` and ``vmax`` to set the floor and ceiling of the color map like so: ::

    >>> my_sequence.plot(vmin=0, vmax=300)

Now that we have created our IRISMapCubeSequence object, we can use one of the followed methods
to manipulate the data.

Exposure Time Correction
------------------------

We can apply the exposure time correction to data and to the uncertainty and
this method adjusts the unit for each IRISMapCube objects inside our IRISMapCubeSequence. ::

    >>> my_sequence.apply_exposure_time_correction()

We also can undo the exposure time correction by doing: ::

    >>> my_sequence.apply_exposure_time_correction(undo=True)

The correction is only applied (undone) if the object's unit doesn't (does) already
include inverse time. This can be overridden so that correction is applied (undone)
regardless of unit by setting ``force=True``. Use one of the two lines above to apply
(undone) by using the force kwarg. ::

    >>> my_sequence.apply_exposure_time_correction(force=True)
    >>> my_sequence.apply_exposure_time_correction(undo=True, force=True)

.. _NDCubeSequence: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html
.. _IRIS: http://iris.lmsal.com/search/
