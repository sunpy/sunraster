An Introduction to sunraster
============================

.. warning::
    
    sunraster is still under heavy development and has not yet seen its first
    release. If you want to help develop and test sunraster, follow these
    `installation instructions
    <https://github.com/sunpy/irispy/wiki/sunraster-Installation-Instructions>`_.

What is sunraster?
------------------

sunraster is a free, open-source, SunPy-affiliated package that provides
tools to manipulate and visualize slit spectrograph data using the Python
programming language.  The sunraster classes link observations
with various forms of supporting data including: measurement
uncertainties; units; a data mask to mark pixels with
unreliable or unphysical data values; WCS (World Coordinate System)
transformations that describe the position, wavelengths and times
represented by the pixels; and general metadata.  These classes also
provide methods for applying and removing exposure time correction to/from 
the observations. Moreover, because the data unit is linked
to the object, it is always obvious what unit the data is in.  This
saves scientists the hassle of performing important, but laborious and
repetitive data conversions and avoid confusion by always tracking the
unit of the data through those conversions.  This leads to more
efficient and accurate science.

sunraster classes inherit more fundamental functionality from the
`ndcube`_ package.  A powerful, generic slicing API (application
programmable interface) allows users to manipulate the same data
object as if it were 3D (time, position along slit, wavelength) or 
4D (raster scan number, slit step, position along slit, wavelength), 
which is very useful when dealing with scanning slit-spectrograph data. 
The API simultaneously slices not only the data, but the uncertainties,
data mask, and WCS transformations leading to faster and less
error-prone data analysis.  The sunraster classes also inherit the
ability to crop by real world coordinates --- useful when locating a
region of interest using information from other observatories --- and a
visualization suite which allows users to easily and intuitively
visually inspect their data.

This guide will explains the capabilities offered by sunraster
and how to utilize them.  It will describe the different data classes,
as well as how to install sunraster, contact the development team, and
contribute to the package.  sunraster is open-source and
community-developed and we are always glad to welcome new contributors
and users.

.. _ndcube: http://docs.sunpy.org/projects/ndcube/en/stable/
