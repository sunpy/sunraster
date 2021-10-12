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
