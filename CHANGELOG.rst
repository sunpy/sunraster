0.4.1 (2022-05-24)
==================

Breaking Changes
----------------

- Increased the minimum version of ``sunpy`` to 4.0.0

0.4.0 (2022-03-08)
==================

Breaking Changes
----------------

- Removed IRIS reader, you will want to install and use ``irispy-lmsal`` instead.
- Removed support for Python 3.7. (`#198 <https://github.com/sunpy/sunraster/pull/198>`__)


0.3.0 (2021-11-19)
==================

Breaking Changes
----------------

- In IRIS spectrograph read, move all metadata to the meta objects of the raster cubes. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- Remove extra_coords keyword from `~sunraster.spectrogram.SpectrogramCube.__init__` in accordance with new ndcube 2.0 API.
  Extra coords can by added through the ndcube ExtraCoords.add API which is new in ndcube 2.0. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- In IRIS spectrograph reader, all extra coords except time have been moved to the meta object. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- Removed ``lon`` and ``lat`` properties from all objects in sunraster. (`#184 <https://github.com/sunpy/sunraster/pull/184>`__)


New Features
------------

- Create new property `~sunraster.spectrogram.SpectrogramSequence.celestial`, on `~sunraster.spectrogram.SpectrogramSequence` to return a `~astropy.coordinates.SkyCoord` holding the celestial world coords of the pixels. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- Create new property `~sunraster.spectrogram.SpectrogramCube.celestial`, on `~sunraster.spectrogram.SpectrogramCube` to return a `~astropy.coordinates.SkyCoord` holding the celestial world coords of the pixels. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- Create a new `~sunraster.instr.iris.IRISSGMeta` metadata object. (`#182 <https://github.com/sunpy/sunraster/pull/182>`__)
- Added a sliceable ``Meta`` class for axis-associated metadata. (`#184 <https://github.com/sunpy/sunraster/pull/184>`__)


0.2.0 (2021-01-28)
==================

Features
--------

- Include a base time to output of `sunraster.SpectrogramCube.time` when time is derived from WCS and a recognized base time can be found in meta. (`#168 <https://github.com/sunpy/sunraster/pull/168>`__)
- Add optional instrument_axes attribute to SpectrogramCube to enable users to keep track of axes (including through slicing) when axes may have a significance not fully described by the world axis physical types. (`#169 <https://github.com/sunpy/sunraster/pull/169>`__)
- Create new Metadata classes for defining mapping of metadata from instrument-specific files to a general metedata API. Includes a specific mapping for SolO/SPICE. (`#171 <https://github.com/sunpy/sunraster/pull/171>`__)
- Replace RasterSequence world_axis_physical_type properties with versions using NDCubeSequence.array_axis_physical_types. (`#173 <https://github.com/sunpy/sunraster/pull/173>`__)
- Provide functions to read SPICE file. Also refactor Meta class to be dict-like. (`#173 <https://github.com/sunpy/sunraster/pull/173>`__)
- Enable SPICE FITS reader to handle multiple files. (`#178 <https://github.com/sunpy/sunraster/pull/178>`__)

Bug Fixes
---------

- Bump min ndcube version to fix bug caused when OS is bot 64-bit. (`#162 <https://github.com/sunpy/sunraster/pull/162>`__)
- Stop `~sunraster.spectrogram_sequence.SpectrogramSequence` crashing when time coord not 1-D. (`#178 <https://github.com/sunpy/sunraster/pull/178>`__)
- Allow SPICE FITS reader to read handle dumbbell windows. (`#178 <https://github.com/sunpy/sunraster/pull/178>`__)
- Ensure args are passed correctly to NDCube constructor by SpectrogramCube by entering them as kwargs instead of ordered args. (`#179 <https://github.com/sunpy/sunraster/pull/179>`__)

Trivial/Internal Changes
------------------------

- Altered names of some SPICEMeta properties. (`#178 <https://github.com/sunpy/sunraster/pull/178>`__)
