from collections import namedtuple

import ndcube.utils.sequence
import numpy as np

SequenceSlice = namedtuple("SequenceSlice", "sequence_index slit_step_item")
"""
Define SequenceSlice named tuple of length 2. Its attributes are:
sequence_index: an int giving the index of a cube within an NDCubeSequence.
slit_step_item: slice of int index of to be to be applied to the common
axis of the cube.
"""

SequenceItem = namedtuple("SequenceItem", "sequence_index cube_item")
"""
Define SequenceItem named tuple of length 2. Its attributes are:
sequence_index: an int giving the index of a cube within an NDCubeSequence.
cube_item: item (int, slice, tuple) to be applied to cube identified
by sequence_index attribute.
"""


def _slice_sequence_as_SnS(sequence, item):
    """
    Enables RasterSequence instance to be indexed as if it were a single sit-
    and-stare cube.

    Parameters
    ----------
    sequence: `sunraster.RasterSequence`
        The sequence to be sliced.

    item: `int`, `slice` or `tuple` of `int` and/or `slice`.
        The cube-like item.  If tuple, length must be <= number
        of dimensions in single raster within sequence.

    Example
    -------
    >>> # Say we have three Cubes each cube has slit_step_axis=1 is time and shape=(3,3,3)
    >>> data_list = [cubeA, cubeB, cubeC] # doctest: +SKIP
    >>> cs = NDCubeSequence(data_list, meta=None, slit_step_axis=1) # doctest: +SKIP
    >>> # return zeroth time slice of cubeB in via normal CubeSequence indexing.
    >>> cs[1,:,0,:] # doctest: +SKIP
    >>> # Return same slice using this function
    >>> _index_sequence_as_cube(cs, (slice(0, cubeB.shape[0]), 0,
    ...                             (slice(0, cubeB.shape[2])) # doctest: +SKIP
    """
    # Convert item to a list of single Raster items for each relevant raster.
    n_slit_steps = np.array([c.data.shape[sequence._slit_step_axis] for c in sequence.data])
    sequence_items = convert_SnS_item_to_sequence_items(item, sequence._slit_step_axis,
                                                        n_slit_steps)
    # Use sequence items to slice NDCubeSequence.
    return ndcube.utils.sequence.slice_sequence_by_sequence_items(sequence, sequence_items)


def convert_SnS_item_to_sequence_items(SnS_item, slit_step_axis, n_slit_steps):
    """
    Converts an item input to RasterSequence.slice_as_SnS to a list of
    SequenceSlice objects.

    Parameters
    ----------
    SnS_item: `int`, `slice`, of `tuple` of `int and/or `slice`.
        Item compatible with input to RasterSequence.slice_as_SnS.

    slit_step_axis: `int`
        Data axis of Rasters corresponding to the slit step axis.

    n_slit_steps: `numpy.ndarray`
        Number of slit steps in each individual Raster.

    Returns
    -------
    sequence_items: `list` of SequenceItem `namedtuple`.
        The slice/index items for each relevant Raster within the RasterSequence
        which together represent the original input slice/index item.
    """
    invalid_item_error_message = "Invalid index/slice input."
    # Case 1: Item is int and slit step axis is 0.
    if isinstance(SnS_item, int):
        if slit_step_axis != 0:
            raise ValueError("Input can only be indexed with an int if "
                             "RasterSequence's slit step axis is 0. common "
                             "axis = {}".format(slit_step_axis))
        else:
            # Derive list of SequenceSlice objects that describes the
            # SnS_item in regular slicing notation.
            sequence_slices = [_convert_SnS_index_to_sequence_slice(
                SnS_item, n_slit_steps)]
            all_axes_item = None
    # Case 2: Item is slice and slit step axis is 0.
    elif isinstance(SnS_item, slice):
        if slit_step_axis != 0:
            raise ValueError("Input can only be sliced with a single slice if "
                             "CubeSequence's slit step axis is 0. common "
                             "axis = {}".format(slit_step_axis))
        else:
            # Derive list of SequenceSlice objects that describes the
            # SnS_item in regular slicing notation.
            # First ensure None types within slice are replaced with appropriate ints.
            sequence_slices = _convert_SnS_slice_to_sequence_slices(
                SnS_item, n_slit_steps)
            all_axes_item = None
    # Case 3: Item is tuple.
    elif isinstance(SnS_item, tuple):
        # Check item is long enough to include slit step axis.
        if len(SnS_item) < slit_step_axis + 1:
            raise ValueError("Input item not long enough to include slit step axis."
                             "Must have length > {}".format(slit_step_axis))
        # Based on type of slice/index in the slit step axis position of
        # the SnS_item, derive list of SequenceSlice objects that
        # describes the SnS_item in regular slicing notation.
        if isinstance(SnS_item[slit_step_axis], int):
            sequence_slices = [_convert_SnS_index_to_sequence_slice(
                SnS_item[slit_step_axis], n_slit_steps)]
        elif isinstance(SnS_item[slit_step_axis], slice):
            sequence_slices = _convert_SnS_slice_to_sequence_slices(
                SnS_item[slit_step_axis], n_slit_steps)
        else:
            raise ValueError(invalid_item_error_message)
        all_axes_item = SnS_item
    else:
        raise TypeError("Unrecognized item type.")
    # Convert the sequence slices, that only describe the slicing along
    # the sequence axis and slit step axis to sequence items which
    # additionally describe how the non-common cube axes should be sliced.
    sequence_items = [_convert_sequence_slice_to_sequence_item(sequence_slice, slit_step_axis,
                                                               SnS_item=all_axes_item)
                      for sequence_slice in sequence_slices]
    return sequence_items


def _convert_SnS_index_to_sequence_slice(SnS_index, n_slit_steps):
    """
    Converts a cube-like index of an NDCubeSequence to indices along the
    sequence and common axes.
    Parameters
    ----------
    SnS_index: `int`
        Cube-like index of NDCubeSequence
    n_slit_steps: iterable of `int`
        Length of each cube along slit step axis.
    Returns
    -------
    sequence_slice: SequenceSlice `namedtuple`.
        First element gives index of cube along sequence axis.
        Second element each index along slit step axis of relevant cube.
    """
    # Derive cumulative lengths of cubes along slit step axis.
    cumul_slit_steps = np.cumsum(n_slit_steps)
    # If SnS_index is within 0th cube in sequence, it is
    # simple to determine the sequence and slit step axis indices.
    try:
        index_in_0th_raster = SnS_index < cumul_slit_steps[0]
    except TypeError as err:
        none_not_int_error_messages = [
            "'>' not supported between instances of 'int' and 'NoneType'",
            "unorderable types: int() > NoneType()"]
        if err.args[0] in none_not_int_error_messages:
            index_in_0th_raster = True
        else:
            raise err
    if index_in_0th_raster:
        sequence_index = 0
        slit_step_index = SnS_index
    # Else use more in-depth method.
    else:
        # Determine the index of the relevant cube within the sequence
        # from the cumulative slit step axis cube lengths.
        sequence_index = np.where(cumul_slit_steps <= SnS_index)[0][-1]
        if SnS_index > cumul_slit_steps[-1] - 1:
            # If the cube is out of range then return the last slit step axis index.
            slit_step_index = n_slit_steps[-1]
        else:
            # Else use simple equation to derive the relevant slit step axis index.
            slit_step_index = SnS_index - cumul_slit_steps[sequence_index]
        # sequence_index should be plus one as the sequence_index earlier is
        # previous index if it is not already the last cube index.
        if sequence_index < cumul_slit_steps.size - 1:
            sequence_index += 1
    # Return sequence and cube indices.  Ensure they are int, rather
    # than np.int64 to avoid confusion in checking type elsewhere.
    if slit_step_index is not None:
        slit_step_index = int(slit_step_index)
    return SequenceSlice(int(sequence_index), slit_step_index)


def _convert_SnS_slice_to_sequence_slices(SnS_slice, n_slit_steps):
    """
    Converts slit step axis slice input to NDCubeSequence.index_as_cube to a list
    of sequence indices.
    Parameters
    ----------
    SnS_slice: `slice`
        Slice along slit step axis in NDCubeSequence.index_as_cube item.
    n_slit_steps: iterable of `int`
        Length of each cube along slit step axis.
    Returns
    -------
    sequence_slices: `list` of SequenceSlice `namedtuple`.
        List sequence slices (sequence axis, slit step axis) for each element
        along slit step axis represented by input SnS_slice.
    """
    # Ensure any None attributes in input slice are filled with appropriate ints.
    cumul_slit_steps = np.cumsum(n_slit_steps)
    # Determine sequence indices of cubes included in cube-like slice.
    SnS_indices = np.arange(cumul_slit_steps[-1])[SnS_slice]
    n_SnS_indices = len(SnS_indices)
    one_step_sequence_slices = np.empty(n_SnS_indices, dtype=object)
    # Define array of ints for all indices along slit step axis.
    # This is restricted to range of interest below.
    sequence_int_indices = np.zeros(n_SnS_indices, dtype=int)
    for i in range(n_SnS_indices):
        one_step_sequence_slices[i] = _convert_SnS_index_to_sequence_slice(SnS_indices[i],
                                                                           n_slit_steps)
        sequence_int_indices[i] = one_step_sequence_slices[i].sequence_index
    unique_index = np.sort(np.unique(sequence_int_indices, return_index=True)[1])
    unique_sequence_indices = sequence_int_indices[unique_index]
    # Convert start and stop cube-like indices to sequence indices.
    first_sequence_index = _convert_SnS_index_to_sequence_slice(SnS_slice.start, n_slit_steps)
    last_sequence_index = _convert_SnS_index_to_sequence_slice(SnS_slice.stop, n_slit_steps)
    # Since the last index of any slice represents
    # 'up to but not including this element', if the last sequence index
    # is the first element of a new cube, elements from the last cube
    # will not appear in the sliced sequence.  Therefore for ease of
    # slicing, we can redefine the final sequence index as the penultimate
    # cube and its slit step axis index as beyond the range of the
    # penultimate cube's length along the slit step axis.
    if (last_sequence_index.sequence_index > first_sequence_index.sequence_index and
            last_sequence_index.slit_step_item == 0):
        last_sequence_index = SequenceSlice(last_sequence_index.sequence_index - 1,
                                            n_slit_steps[last_sequence_index.sequence_index - 1])
    # Iterate through relevant cubes and determine slices for each.
    # Do last cube outside loop as its end index may not correspond to
    # the end of the cube's slit step axis.
    if SnS_slice.step is None:
        step = 1
    else:
        step = SnS_slice.step
    sequence_slices = []
    slit_step_start_index = first_sequence_index.slit_step_item
    j = 0
    while j < len(unique_sequence_indices) - 1:
        # Let i be the index along the sequence axis of the next relevant cube.
        i = unique_sequence_indices[j]
        # Determine last slit step axis index for this cube.
        slit_step_last_index = n_slit_steps[i] - ((n_slit_steps[i] - slit_step_start_index) % step)
        # Generate SequenceSlice for this cube and append to list.
        sequence_slices.append(SequenceSlice(
            i, slice(slit_step_start_index,
                     min(slit_step_last_index + 1, n_slit_steps[i]), step)))
        # Determine first slit step axis index for next cube.
        if n_slit_steps[i] == slit_step_last_index:
            slit_step_start_index = step - 1
        else:
            slit_step_start_index = \
                step - (((n_slit_steps[i] - slit_step_last_index) % step) +
                        cumul_slit_steps[unique_sequence_indices[j + 1] - 1] -
                        cumul_slit_steps[i])
        # Iterate counter.
        j += 1
    # Create slice for last cube manually.
    sequence_slices.append(SequenceSlice(
        unique_sequence_indices[j],
        slice(slit_step_start_index, last_sequence_index.slit_step_item, step)))
    return sequence_slices


def _convert_sequence_slice_to_sequence_item(sequence_slice, slit_step_axis, SnS_item=None):
    """
    Converts a SequenceSlice named tuple to a SequenceItem named tuple.

    Parameters
    ----------
    sequence_slice: SequenceSlice `namedtuple`
        0th element gives index of along sequence axis.
        1st element each index along slit step axis of relevant cube.
        Must be same format as output from _convert_SnS_index_to_sequence_slice.

    slit_step_axis: `int`
        Slit step axis as defined in RasterSequence.

    SnS_item: `None` or `tuple` of `slice` and/or `int` objects (Optional)
        The original item input to `RasterSequence.slice_as_SnS` including the
        slices/indices of non-slit-step axes of rasters within sequence.  If None, a
        tuple of slice(None) objects is generated  long enough so that the last
        element in the tuple corresponds to the common axis and is set to the
        1st (0-based counting) the sequence_index input, above.  This tuple is
        then set to the cube_item attribute of the output `SequenceSlice` object.

    Returns
    -------
    sequence_item: SequenceItem `namedtuple`.
        Describes sequence index of a Raster within a RasterSequence and the
        slice/index item to be applied to the whole Raster.
    """
    if SnS_item is None and slit_step_axis == 0:
        sequence_item = SequenceItem(sequence_slice.sequence_index,
                                     sequence_slice.slit_step_item)
    else:
        # Create mutable version of SnS_item.
        try:
            SnS_item_list = list(SnS_item)
        except TypeError as err:
            if err.args[0] == "'NoneType' object is not iterable":
                SnS_item_list = []
            else:
                raise err
        # Make sure SnS_item is long enough to include common axis
        while len(SnS_item_list) <= slit_step_axis:
            SnS_item_list.append(slice(None))
        # Create new sequence slice
        SnS_item_list[slit_step_axis] = sequence_slice.slit_step_item
        sequence_item = SequenceItem(sequence_slice.sequence_index, tuple(SnS_item_list))
    return sequence_item
