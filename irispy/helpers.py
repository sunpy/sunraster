import copy
import warnings


def _truncate(object, parameter, axis=None):
    """
    Helper function to truncate the IRISRaster instance.

    Parameters
    ----------

    object : `~irispy.spectrograph.IRISRaster`.

    parameter : `slice` or `tuple`.

    axis: String.
        Can be time/raster_position/slit_axis.
        Return a truncated IRISRaster object of the given axis.
    """
    new_raster = copy.deepcopy(object)
    spectral_windows = new_raster.spectral_windows['name']
    if axis == None:
        raise ValueError("Axis must be defined")
    for window in spectral_windows:
        if axis == 'time' or axis == 'raster_axis' or axis == 'raster_position':
            new_raster.data[window] = new_raster.data[window].sel(raster_axis=parameter)
        if axis == 'slit_axis' or axis == 'slit_position':
            new_raster.data[window] = new_raster.data[window].sel(slit_axis=parameter)

    return new_raster