********************
Welcome to sunraster
********************

Welcome to the ``sunraster`` User Guide.

``sunraster`` is a free, open-source, community-developed, SunPy-affiliated package that provides tools to manipulate and visualize slit spectrograph data.
We are always glad to welcome new contributors and users.

The ``sunraster`` classes link observations with various forms of supporting data including: measurement uncertainties; units; a data mask to mark pixels with unreliable or un-physical data values; WCS (World Coordinate System) transformations that describe the position, wavelengths and times represented by the pixels; and general metadata.
These classes also provide methods for applying and removing exposure time corrections to/from the observations.
Moreover, because the data unit is linked to the object, it is always obvious what unit(s) the data is in.
This saves scientists the hassle of performing important, but laborious and repetitive data conversions and avoid confusion by always tracking the unit(s) of the data through those conversions.

The ``sunraster`` classes inherit more fundamental functionalities from the `ndcube`_ package.
These include a powerful, generic slicing API (application programmable interface) allowing users to manipulate the same data object as though it were 3D (time, position along slit, wavelength) or 4D (raster scan number, slit step, position along slit, wavelength), which is very useful when dealing with scanning slit-spectrograph data.
The API simultaneously slices not only the data, but the uncertainties, data mask, and WCS transformations leading to faster and less error-prone data analysis.
The ``sunraster`` classes also inherit the ability to crop by real world coordinates --- useful when locating a region of interest using information from other observatories --- and a visualization suite which allows users to easily and intuitively visually inspect their data.

This guide explains the capabilities offered by ``sunraster`` and how to utilize them.
It will describe the different data classes, as well as how to install ``sunraster``, contact the development team, and contribute to the package.

.. _ndcube: https://docs.sunpy.org/projects/ndcube/en/stable/

Contents
========

* :ref:`genindex`
* :ref:`modindex`

.. toctree::
   :maxdepth: 2

   installation
   data_types/index
   api
   changelog
