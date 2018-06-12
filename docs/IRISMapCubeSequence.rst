.. _IRISMapCubeSequence:

===================
IRISMapCubeSequence
===================

This class provides access to IRIS Slit-jaw image (SJI) files as a list of
IRISMapCube objects, where these sub-cubes are 2D or 3D described by their
 own World Coordinate System (WCS).

This class inherits from NDCubeSequence_, so many methods of IRISMapCubeSequence are linked
to their description in NDCubeSequence_ documentation.

Initialisation
--------------

The initialisation of an IRISMapCubeSequence is very similar to that of :ref:`IRISMapCube`,
via the function ``read_iris_sji_level2_fits``. This function takes as argument a file name (string) or a list of filenames.

Let us assume that we want to load an IRISMapCubeSequence object with two fits files called ``my_fits_file_0`` and ``my_fits_file_1``. We can create an IRISMapCubeSequence
with as many fits files as necessary, but ideally with from the same observation (e.g. observations interrupted by eclipses). To load an IRISMapCubeSequence object as ``my_sequence``, we do: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files)

As for IRISMapCube, if you don't have a lot of RAM memory or if you are loading a huge file,
we recommend to use the ``memmap=True`` keyword. By using it, you will only load the data when needed. However, some methods that require all data in memory will not be accessible. You can use memmap by doing: ::

    >>> from irispy import read_iris_sji_level2_fits
    >>> my_fits_files = [my_fits_file_0, my_fits_file_1]
    >>> my_sequence = read_iris_sji_level2_fits(my_fits_files, memmap=True)

Once we have our ``IRISMapCube`` object, we can quickly inspect it. For that, we can use the ``__repr__`` property of our object, which shows some metadata. This property is called just by writing the name of our object in the command line: ::

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

Here we can find some information about our object such as when the observation was made,
the OBS ID, its dimensions, etc.

Class Structure
---------------

Sequence Attributes
^^^^^^^^^^^^^^^^^^^

Our IRISMapCubeSequence object, here called ``my_sequence`` is now created and inspected, we can look into its structure.

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


The return value is a list with the representation of each cube inside. For an upcoming work, we
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

The return value is the first IRISMapCube object so now we can work on
this cube with everything we saw previously for an IRISMapCube object (see :ref:`Cube_Attribute`).

Meta
""""

Metadata for the first IRISMapCube object is saved under ``.meta``: ::

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

As an IRISMapCubeSequence object is a list of IRISMapCube objects, we can do the same operations
listed in the :ref:`IRISMapCube` documentation. To access to the information
stored in a IRISMapCube object, we can select the first IRISMapCube object by writting
``my_sequence[0]`` or the Nth object with ``my_sequence[N-1]``. As they are IRISMapCube
objects, we can access to their information by just replacing ``my_cube`` by the name
of the cube we want to inspect.

Dimensions
----------

The IRISMapCubeSequence object inherits from NDCubeSequence_, so we can use the two
properties of NDCubeSequence_ which allow us to get the data shape and the axis types of
our IRISMapCubeSequence object. However, to be consistent with the methods ``IRISMapCube.dimensions`` and ``IRISMapCube.world_axis_physical_types`` methods,
these methods have been rewritten for the IRISMapCubeSequence objects to have the same format as IRISMapCube. To see that, we can do: ::

    >>> my_sequence.dimensions
    <Quantity [  98., 1095., 1018.] pix>

    >>> my_sequence.world_axis_physical_types
    ('time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon')

Cropping and Indexing
---------------------

One of the most powerful capabilities of IRISMapCubeSequence, coming from NDCubeSequence,
is the slicing. To slice the sequence, we can slice IRISMapCubeSequence using
array-like indices or by coordinates. As IRISMapCubeSequence object inherits from NDCubeSequence_, we can use the processes described in the NDCubeSequence.Slicing_ section.

Manipulating the Data
---------------------

We can manipulate an IRISMapCubeSequence object with the methods listed below.


Exposure Time Correction
^^^^^^^^^^^^^^^^^^^^^^^^

This method scales the data from data number (DN) units to DN per second, thereby correcting for any changes in exposure time during an observation and allowing a better comparison between different observations. It works in the same way as for IRISMapCube, see :ref:`Exposure_Time_Correction`.


Dust particle mask
^^^^^^^^^^^^^^^^^^

This method adds the dust pixels to the invalid pixel mask. It works similarly as for IRISMapCube, see :ref:`Dust_Particle_Mask`. The main difference is that the ``mask`` and ``dust_masked`` attributes are defined for each object in the sequence (``my_sequence[N].mask``).

Visualisation
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
