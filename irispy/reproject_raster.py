"""
Reproject a series of IRIS raster scans in to one NDCube.
"""
import numpy as np

import astropy.units as u
import sunpy.coordinates
from astropy.time import Time
from astropy.wcs import WCS
from reproject.interpolation import reproject_interp

from irispy.spectrograph import read_iris_spectrograph_level2_fits
from ndcube import NDCube


def stack_spectrogram_sequence(cube_sequence):
    """
    Given a sequence of IRIS rasters stack them into a single `ndcube.NDCube`.

    .. warning::

        This is intended to be used for plotting only, it will not preserve the
        flux in the image.

    Parameters
    ----------
    cube_sequence : `irispy.spectrogram.IRISSpectrogramCubeSequence`

    Returns
    -------
    `ndcube.NDCube`: A 4D cube with a new time dimension

    """

    if len(cube_sequence.data) == 1:
        raise ValueError("No point doing this to one raster")
    target_wcs = cube_sequence[0].wcs
    target_shape = cube_sequence[0].data.shape

    new_arrays = [cube_sequence[0].data]
    times = [cube_sequence[0].extra_coords['time']['value'][0]]
    for cube in cube_sequence.data[1:]:
        array, footprint = reproject_interp((cube.data, cube.wcs),
                                            target_wcs, shape_out=target_shape,
                                            independent_celestial_slices=True)
        new_arrays.append(array)
        times.append(cube.extra_coords['time']['value'][0])

    times = Time(times)
    dts = times[1:] - times[:-1]

    if u.allclose(dts[0].to(u.s), dts.to(u.s), atol=0.5*u.s):
        dt = dts[0]
    else:
        raise ValueError("Can't handle tabular wcs")

    out_wcs = WCS(naxis=4)
    out_wcs.wcs.crpix = list(target_wcs.wcs.crpix) + [0]
    out_wcs.wcs.crval = list(target_wcs.wcs.crval) + [0]
    out_wcs.wcs.cdelt = list(target_wcs.wcs.cdelt) + [dt.to(u.s).value]
    out_wcs.wcs.ctype = list(target_wcs.wcs.ctype) + ['TIME']
    out_wcs.wcs.cunit = list(target_wcs.wcs.cunit) + ['s']

    pc = np.identity(4)
    pc[1:, 1:] = target_wcs.wcs.pc

    out_wcs.wcs.pc = pc

    out_array = np.stack(new_arrays, axis=0)
    return NDCube(out_array, out_wcs)
