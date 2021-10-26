import copy

import numpy as np
import pytest

from sunraster.meta import Meta


def assert_metas_equal(test_input, expected_output):
    if type(test_input) is not type(expected_output):
        raise AssertionError(
            "input and expected are of different type. " f"input: {type(test_input)}; expected: {type(expected_output)}"
        )
    multi_element_msg = "more than one element is ambiguous"
    if isinstance(test_input, Meta) and isinstance(expected_output, Meta):
        # Check keys are the same.
        assert test_input.keys() == expected_output.keys()

        # Check shapes are the same.
        if test_input.shape is None or expected_output.shape is None:
            assert test_input.shape == expected_output.shape
        else:
            assert np.allclose(test_input.shape, expected_output.shape)

        # Check values and axes are the same.
        for test_value, expected_value in zip(test_input.values(), expected_output.values()):
            try:
                assert test_value == expected_value
            except ValueError as err:
                if multi_element_msg in err.args[0]:
                    assert np.allclose(test_value, expected_value)
        # Check axes are the same.
        for test_axis, expected_axis in zip(test_input.axes.values(), expected_output.axes.values()):
            assert (test_axis is None and expected_axis is None) or all(test_axis == expected_axis)
    else:
        if not (test_input is None and expected_output is None):
            assert test_input.keys() == expected_output.keys()
            for key in list(test_input.keys()):
                assert test_input[key] == expected_output[key]


@pytest.fixture
def basic_meta_values():
    return {
        "a": "hello",
        "b": list(range(10, 25, 10)),
        "c": np.array([[1, 2, 3, 4], [10, 20, 30, 40], [100, 200, 300, 400]]),
        "d": list(range(3, 13, 3)),
    }


@pytest.fixture
def basic_comments():
    return {
        "a": "Comment A",
        "b": "Comment B",
        "c": "Comment C",
    }


@pytest.fixture
def basic_axes():
    return {
        "b": 0,
        "c": (1, 2),
        "d": (2,),
    }


@pytest.fixture
def basic_data_shape():
    return (2, 3, 4, 5)


@pytest.fixture
def basic_meta(basic_meta_values, basic_comments, basic_axes, basic_data_shape):
    return Meta(basic_meta_values, basic_comments, basic_axes, basic_data_shape)


@pytest.fixture
def no_shape_meta():
    return Meta({"a": "hello"})


def test_shape(basic_meta, basic_data_shape):
    meta = basic_meta
    shape = np.asarray(basic_data_shape)
    assert all(meta.shape == shape)


def test_slice_axis_with_no_meta(basic_meta):
    meta = basic_meta
    output = meta[:, :, :, 0]
    expected = copy.deepcopy(meta)
    expected._data_shape = meta._data_shape[:-1]
    assert_metas_equal(output, expected)


def test_slice_away_independent_axis(basic_meta):
    meta = basic_meta
    # Get output
    item = 0
    output = meta[item]
    # Build expected result.
    values = dict([(key, value) for key, value in meta.items()])
    values["b"] = values["b"][0]
    comments = meta.comments
    axes = dict([(key, axis) for key, axis in meta.axes.items()])
    del axes["b"]
    axes["c"] -= 1
    axes["d"] -= 1
    shape = meta.shape[1:]
    print(values, comments, axes, shape)
    expected = Meta(values, comments, axes, shape)
    # Compare output and expected.
    assert_metas_equal(output, expected)


def test_slice_dependent_axes(basic_meta):
    meta = basic_meta
    print(meta["a"])
    # Get output
    output = meta[:, 1:3, 1]
    print(meta["a"])
    # Build expected result.
    values = dict([(key, value) for key, value in meta.items()])
    values["c"] = values["c"][1:3, 1]
    values["d"] = values["d"][1]
    comments = meta.comments
    axes = dict([(key, axis) for key, axis in meta.axes.items()])
    axes["c"] = 1
    del axes["d"]
    shape = np.array([2, 2, 5])
    expected = Meta(values, comments, axes, shape)
    # Compare output and expected.
    assert_metas_equal(output, expected)


def test_slice_by_str(basic_meta):
    meta = basic_meta
    assert meta["a"] == "hello"
    assert meta["b"] == list(range(10, 25, 10))


def test_add1(basic_meta):
    meta = basic_meta
    name = "z"
    value = 100
    comment = "Comment E"
    meta.add(name, value, comment, None)
    assert name in meta.keys()
    assert meta[name] == value
    assert meta.comments[name] == comment
    assert meta.axes.get(name, None) is None


def test_add2(basic_meta):
    meta = basic_meta
    name = "z"
    value = list(range(2))
    axis = 0
    meta.add(name, value, None, axis)
    assert name in meta.keys()
    assert meta[name] == value
    assert meta.comments.get(name, None) is None
    assert meta.axes[name] == np.array([axis])


def test_add_overwrite(basic_meta):
    meta = basic_meta
    name = "a"
    value = "goodbye"
    meta.add(name, value, None, None, overwrite=True)
    assert meta[name] == value


def test_add_overwrite_error(basic_meta):
    meta = basic_meta
    with pytest.raises(KeyError):
        meta.add("a", "world", None, None)


def test_add_axis_without_shape(no_shape_meta):
    meta = no_shape_meta
    with pytest.raises(TypeError):
        meta.add("z", [100], axis=0)


def test_remove(basic_meta):
    meta = basic_meta
    name = "b"
    meta.remove(name)
    assert name not in meta.keys()
    assert name not in meta.comments.keys()
    assert name not in meta.axes.keys()
