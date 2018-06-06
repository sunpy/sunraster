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

Now, that we have create our IRISMapCubeSequence object, we can check if the creation has been
successfully done. We can use the representation property of our object, where are stored
many information, to do that. This property is called just by writing the name of our
object in the console. ::

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

Class Structure
---------------

Sequence Attributes
^^^^^^^^^^^^^^^^^^^

Our IRISMapCubeSequence object, called ``my_sequence`` is now created and checked, we can
be interested in the structure of our object.

- ``my_sequence.data`` : In this attribute, we can find the data array of our sequence.
- ``my_sequence.meta`` : The meta data contains a lot of information which are not inside
  the data array.
- ``my_sequence.common_axis``: We can found some explanation in NDCubeSequence.Common_Axis_
  paragraph.

Data
""""

The ``data`` attribute of a IRISMapCubeSequence object contains the list of all cubes
inside the object. We can see this list by doing: ::

    >>> my_sequence.data
    [
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
         ,
         IRISMapCube
         ---------
         Observatory:		 IRIS
         Instrument:			 SJI
         Bandpass:			 1400.0
         Obs. Start:			 2018-04-26T23:07:22.780000
         Obs. End:			 2018-04-27T01:39:47.122000
         Instance Start:		 2018-04-26T23:08:25.360000
         Instance End:		 2018-04-27T01:37:42.990000
         Total Frames in Obs.:	 49
         IRIS Obs. id:		 3690015104
         IRIS Obs. Description:	 Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x
         Cube dimensions:		 [  49. 1095. 1018.] pix
         Axis Types:			 ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')
         ]


The return is a list with the representation of each cube inside. For an upcoming work, we
will need to select only on cube. We can do that by indexing our sequence from 0 to N-1,
where N is the number of cubes. To select the first cube do: ::

    >>> my_sequence[0]
        IRISMapCube
        ---------
        Observatory:		 IRIS
        Instrument:			 SJI
        Bandpass:			 1330.0
        Obs. Start:			 2018-04-26T23:07:22.780000
        Obs. End:			 2018-04-27T01:39:47.122000
        Instance Start:		 2018-04-26T23:07:22.880000
        Instance End:		 2018-04-26T23:07:22.880000
        Total Frames in Obs.:	 49
        IRIS Obs. id:		 3690015104
        IRIS Obs. Description:	 Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x
        Cube dimensions:		 [1095. 1018.] pix
        Axis Types:			 ('custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

The return is the representation of the first IRISMapCube object so now we can work on
this cube with everything we saw previously for an IRISMapCube object (see :ref:`Cube_Attribute`).

Meta
""""

The ``meta`` attribute is storing the meta dictionary of the first IRISMapCube object.
We can access to it by doing: ::

    >>> my_sequence.meta
    {'ENDOBS': datetime.datetime(2018, 4, 27, 1, 39, 47, 122000),
     'INSTRUME': 'SJI',
     'NBFRAMES': 49,
     'OBSID': '3690015104',
     'OBS_DESC': 'Very large sit-and-stare 0.3x175 1s  C II   Si IV   Mg II h/k Deep x',
     'STARTOBS': datetime.datetime(2018, 4, 26, 23, 7, 22, 780000),
     'TELESCOP': 'IRIS',
     'TWAVE1': 1400.0}

And we can also select only one key (eg. ``OBSID``) with the line: ::

    >>> my_sequence.meta['OBSID']
    '3690015104'

.. _Cube_Attribute:

Cube Attributes
^^^^^^^^^^^^^^^

As an IRISMapCubeSequence object is a list of IRISMapCube objects, we can do the same things
that we have seen previously in the IRISMapCube documentation. To access to the information
stored in a IRISMapCube object, we can select the first IRISMapCube object by writting
``my_sequence[0]`` or the Nth object with ``my_sequence[N-1]``. As they are IRISMapCube
objects, we can access to their information by just replacing ``my_cube`` by the name
of the cube we want to inspect.

Dimensions
----------

The IRISMapCubeSequence object inherits from NDCubeSequence_, so we can use the two
properties of NDCubeSequence_ which allow us to get the data shape and the axis types of
our IRISMapCubeSequence object. However, to stay consistent with the ``IRISMapCube.dimensions``
and ``IRISMapCube.world_axis_physical_types`` methods, these methods have been rewritten for
the IRISMapCubeSequence objects to have the same format as IRISMapCube. To see that,
we can do: ::

    >>> my_sequence.dimensions
    <Quantity [  98., 1095., 1018.] pix>

    >>> my_sequence.world_axis_physical_types
    ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

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

If we take time to look inside the data, we see that some pixels are obscured by dust,
and so do not reflect the emission from the Sun at that location. The ``apply_dust_mask``
method of IRISMapCubeSequence can be used to add the locations of the dust pixels to the mask
of each IRISMapCube inside and so that we can easily exclude them from our analysis. We can
use this method by doing: ::

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

As the IRISMapCubeSequence object inherits from NDCubeSequence_, we can use the plotting
method of NDCubeSequence_ which allow us to see the data in plots or animations. This
method is described in the NDCubeSequence.Plotting_ section.

.. _NDCubeSequence: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html
.. _IRIS: http://iris.lmsal.com/search/
.. _NDCubeSequence.Common_Axis: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#common-axis
.. _NDCubeSequence.Dimensions: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#dimensions
.. _NDCubeSequence.Slicing: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#slicing
.. _NDCubeSequence.Plotting: http://docs.sunpy.org/projects/ndcube/en/stable/ndcubesequence.html#plotting
