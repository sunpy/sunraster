An Introduction to IRIS and IRISpy
==================================

.. warning::
    
    IRISPy is still under heavy development and has not yet seen its first
    release. If you want to help develop and test IRISPy, follow these
    `installation instructions
    <https://github.com/sunpy/irispy/wiki/IRISpy-Installation-Instructions>`_.

The IRIS Instrument
-------------------

The `Interface Region Imaging Spectrograph`_ (IRIS) is satellite-borne
solar scanning slit spectrograph, funded as part of NASA's
Small Explorer program and launched in June 2013.  It provides
simultaneous UV images and spectra of the regions between the visible
surface of the Sun and its outer atmosphere (photosphere,
chromosphere, transition region, and corona) with 0.33 – 0.4 arcsec
spatial resolution, two-second temporal resolution and 1 km/s velocity
resolution over a field-of-view of up to 175 arcsec × 175 arcsec.
IRIS images are provided by its Slit-Jaw Imager (SJI) in four
passbands: C ii 1330 A, Si iv 1400 A, Mg ii k 2796 A and Mg ii wing
2830 A.  Meanwhile IRIS's spectrograph provides spectra in several
spectral windows in the ranges 1332–1358 A, 1389–1407 A and
2783–2834 A.  For more detail, read the `instrument paper`_.

As a scanning slit-spectrograph, IRIS disperses the
sunlight by passing it through a narrow slit and onto a CCD.  The 
spectrograph can operate in two basic modes: sit-and-stare, where the
slit is aligned with a single position on the Sun; and raster, where
the slit moves sequentially across the Sun perpendicular to the long
axis of the slit in a pre-determined number of steps and step size.
When the last position is reached, the slit returns to the origin
position and starts again.  IRIS is has many different observing
programs involving different numbers of raster steps, step 
sizes, exposure times, etc., making it a powerful and flexible tool for
solar physics.

The complexity of IRIS leads to a variety of different
data products with different dimensionalities and ways in which
they are used by scientists.  Therefore a powerful, yet flexible suite of data
analysis tools is required to enable users to efficiently, reliably and
effectively pursue their science goals.  This is the aim of IRISpy.

What is IRISpy?
---------------

IRISPy is a free, open-source, SunPy-affiliated package that provides
tools to read, manipulate and visualize IRIS data using the Python
programming language.  IRISpy provides a set of classes for handling
both SJI and spectrograph observations.  These link the observations
with various forms of supporting data including: measurement
uncertainties; units; a data mask to mark pixels with
unreliable or unphysical data values; WCS (World Coordinate System)
transformations that describe the position, wavelengths and times
represented by the pixels; and general metadata.  These classes also
provide methods for applying a number of calibration routines
including exposure time correction and conversion between data number,
photons, and energy units.  Moreover, because the data unit is linked
to the object, it is always obvious what unit the data is in.  This
saves scientists the hassle of performing important, but laborious and
repetitive data conversions and avoid confusion by always tracking the
unit of the data through those conversions.  This leads to more
efficient and accurate science.

IRISpy classes inherit more fundamental functionality from the
`ndcube`_ package.  A powerful, generic slicing API (application
programmable interface) allows users to manipulate the same data
object as if it were 3D (time, latitude, wavelength) or 4D (raster
number, longitude, latitude, wavelength), which is very useful
when dealing with scanning slit-spectrograph data.  The 
API simultaneously slices not only the data, but the uncertainties,
data mask, and WCS transformations leading to faster and less
error-prone data analysis.  The IRISpy classes also inherit the
ability to crop by real world coordinates, useful when locating a
region of interest using information from other observatories, and a
visualization suite which allows users to easily and intuitively
visually inspect their data.

This guide will explain in detail the capabilities offered by IRISpy
and how to utilize them.  It will describe the different data classes,
as well as how to install IRISpy, contact the development team, and
contribute to the package.  IRISpy is open-source and
community-developed and we are always glad to welcome new contributors
and users.

.. _Interface Region Imaging Spectrograph: http://iris.lmsal.com/
.. _instrument paper: https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0196&file_type=pdf
.. _ndcube: http://docs.sunpy.org/projects/ndcube/en/stable/
