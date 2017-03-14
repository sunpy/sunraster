# IRISpy Requirements
------------------------

This document is a draft of required/desired functionalities for IRISpy.

## IRISRaster Object

An object to read in, store and manipulate IRIS spectrograph files
(raster or sit-and-stare) from the same OBS campaign.

#### Data Description

IRIS spectrograph data have the following characteristics:
* They contain a number of non-contiguous 3D data cubes, one for each
spectral window in OBS.
* Data from sit-and-stare OBS are stored in a single FITS file with a
single set of WCS values.
* Data from rastering OBS are broken into multiple files, one for each
  raster scan.  Each file has its own set of WCS values.
* Data cubes have the following dimensions:
 * spectral dimension (Dimension along which wavelength varies.)
 * slit dimension (Dimension along length of slit.)
 * raster dimension (Dimension along which time and/or raster slit
 position varies.)
* Data cubes have a pseudo-dimension labelling the spectral window.
  This is referred to as a pseudo-dimension because data may want to
  be selected or manipulated based on this attribute.  But the values
  along this pseudo-dimension are not numerical.
* The spectral dimension is associated with the following coordinate
  systems:
 * wavelength (or physical equivalent, e.g. frequency.)
* The slit dimension is associated with the following coordinate
  systems:
 * helio-projected latitude.  This is dependent on the helio-projected
    longitude also associated with the raster dimension.
 * slit pixel number
* The raster dimension is associated with the following coordinate
  systems:
 * helioprojected longitude.  This is dependent on the
   helio-projected latitude also associated with the slit
   dimension.
 * time
 * raster scan position number, i.e. the slit position number within the
   scan sequence.  Therefore spectra taken at the same position
   within each raster scan will have the same value.  This makes it
   possible to identify spectra taken at nearly the same longitude
   (if solar rotation tracking off) or over nearly the same feature
   on the solar (if solar rotation tracking on).

#### Required Functionalities

* Read in and store IRIS raster/sit-and-stare files of the same OBS
campaign.
* Store non-contiguous data cubes, one for each each spectral window.
* Truncate all spectral windows using criteria applied to common
 dimension(s) with a single command.  Common dimensions are the slit
 and raster dimensions.  The spectral dimension is not common.
 Spectral window is a common pseudo-dimension.
* Dimensions should be indexable by the coordinate systems associated
with them.
* Both data and coordinate systems should be unit aware (desired).
* Convert value and unit of data from DN to photon counts and vice
versa using gain and yield info.
* Convert value and unit of data from DN (or photon counts) to DN/s
(or counts/s) by applying an exposure time correction.
* Calculate fractional uncertainty of intensity data.
* Convert intensity data from DN to radiance.
* Calculate and apply wavelngth correction due to orbital variation.


## IRISSJIMap

An object to read in, store and manipulate IRIS slit-jaw imager (SJI)
data files from the same OBS campaign.

#### Data Description

SJI provides spatially resolved images of the chromosphere and
transition region in a number of passbands.

#### Required Functionalities


#### Implementation Suggestions

* Make IRISSJIMap a subclass of SunPy Map or MapCube.

## IRIS Tools

### IRIS Quicklook tool

### Get IRIS Response
