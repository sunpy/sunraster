.. _IRISMapCubeSequence:

===================
IRISMapCubeSequence
===================

This class represents Slit Jaw Images from IRIS as a list of IRISMapCube, where these
sub-cubes are 2D or 3D described by their own World Coordinate System (WCS).

This class inherits from NDCubeSequence_, so many methods of IRISMapCubeSequence are linked
to their description in NDCubeSequence_ documentation.

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

Class Structure
---------------

Sequence Attributes
^^^^^^^^^^^^^^^^^^^

So now, ``my_sequence`` is an IRISMapCubeSequence object and we can access to a lot of
information like:

- ``my_sequence.data`` : In this attribute, we can find the data array of our sequence.
- ``my_sequence.meta`` : The meta data contains a lot of information which are not inside
  the data array.
- ``my_sequence.common_axis``: We can found some explanation in NDCubeSequence.Common_Axis_
  paragraph.

Cube Attributes
^^^^^^^^^^^^^^^

As an IRISMapCubeSequence object is a list of IRISMapCube objects, we can do the same things
that we have seen previously in the IRISMapCube documentation. To access to the information
stored in a IRISMapCube object, we can select the first IRISMapCube object by writting
``my_sequence[0]`` or the Nth object with ``my_sequence[N-1]``. As they are IRISMapCube
objects, we can access to their information by doing:

- ``my_sequence[0].data`` : The data attribute of the first IRISMapCube object.
- ``my_sequence[1].wcs`` : The WCS attribute of the second IRISMapCube object.
- ``my_sequence[2].mask`` : The mask attribute of the third IRISMapCube object.
- ``my_sequence[3].uncertainty`` : The uncertainty attribute of the fourth IRISMapCube object.

Dimensions
----------

As the IRISMapCubeSequence object inherits from NDCubeSequence_, we can use the two
properties of NDCubeSequence_ which allow us to the data shape and the axis types of
our IRISMapCubeSequence object. These properties are described in the
NDCubeSequence.Dimensions_ section.

Cropping and Indexing
---------------------

One of the most powerful capability of IRISMapCubeSequence, coming from NDCubeSequence,
is the slicing process. To slice the sequence, we can slice the IRISMapCubeSequence with
an Array-like Indexing or we can crop it by the Real World Coordinates. As the
IRISMapCubeSequence object inherits from NDCubeSequence_, we can use the described
processes in the NDCubeSequence.Slicing_ section.

Data Manipulation
-----------------

Now we have our IRISMapCubeSequence object and we know how access to all the information
it contains, we can manipulate the data with the below presented methods.

Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

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

Dust particle mask
^^^^^^^^^^^^^^^^^^

If we take time to look inside the data, we can see that some spots are not relevant with
the data we want to study. These spots are dust particles that can come from space or a
not perfect construction of the satellite. In IRISMapCubeSequence, we can use a method
that will modify each IRISMapCube mask inside our object with the dust particles positions.
We can use this method by doing: ::

    >>> my_sequence.apply_dust_mask()

Now, all our ``my_sequence[N].mask`` contains the dust particles positions and we can use
them to select only the data we want. If we want to remove the dust particles positions from
all masks, we can call again this method with the ``undo`` kwarg. ::

    >>> my_sequence.apply_dust_mask(undo=True)

If we don't remember or we don't know if the dust particles positions are applied or not
in our Nth cube mask, we can check an attribute of our object. ::

    >>> my_sequence[N].dust_masked

If the result is ``True``, the dust particles positions are applied in our this cube.
If the result is ``False``, the dust particles positions are not applied.

By default, all cubes are modified in the same time so if one is ``True``, all other
value must be ``True``. We can, by hand, modify only one mask if we want. To do modify
the Nth mask, we just need to do one of the following lines: ::

    >>> my_sequence[N].apply_dust_mask()
    >>> my_sequence[N].apply_dust_mask(undo=True)

Visualization
-------------

There is two different ways to visualize our IRISMapCubeSequence object. The first one
is to use the representation of our object and the second one is to plot the data of
our object.

Representation
^^^^^^^^^^^^^^

In this part, we can have a look of our IRISMapCubeSequence object by using its
representation property. ::

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

Here, we can found some information about our object like when the observation was made,
the observation ID, its dimensions, etc ...

Plotting
^^^^^^^^

As the IRISMapCubeSequence object inherits from NDCubeSequence_, we can use the plotting
method of NDCubeSequence_ which allow us to see the data in plots or animations. This
method is described in the NDCubeSequence.Plotting_ section.

.. _NDCubeSequence: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html
.. _IRIS: http://iris.lmsal.com/search/
.. _NDCubeSequence.Common_Axis: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#common-axis
.. _NDCubeSequence.Dimensions: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#dimensions
.. _NDCubeSequence.Slicing: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#slicing
.. _NDCubeSequence.Plotting: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting
