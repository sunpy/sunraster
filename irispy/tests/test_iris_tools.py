# -*- coding: utf-8 -*-
from __future__ import print_function

import pytest

import numpy as np
import numpy.testing as np_test
import astropy.units as u

import irispy.iris_tools as iris_tools

DETECTOR_TYPE_KEY = "detector type"

# Arrays of DN
SOURCE_DATA_DN = np.array([[ 0.563,  1.132, -1.343], [-0.719,  1.441, 1.566]])
SOURCE_DATA_DN_1 = np.array([[1, 2, 3], [4, 5, 6]])

# Arrays relating SOURCE_DATA_DN to photons in NUV and FUV
SOURCE_DATA_PHOTONS_NUV = np.array([[ 10.134, 20.376, -24.174], [-12.942, 25.938, 28.188]])
SOURCE_DATA_PHOTONS_FUV = np.array([[ 2.252, 4.528, -5.372], [-2.876, 5.764, 6.264]])

# Arrays relating SOURCE_DATA_DN_1 and photons in SJI
SOURCE_DATA_PHOTONS_SJI_1 = np.array([[18, 36, 54], [72, 90, 108]])

single_exposure_time = 2.
EXPOSURE_TIME = np.zeros(3)+single_exposure_time


@pytest.mark.parametrize("test_input, expected_output", [
    ({DETECTOR_TYPE_KEY: "FUV1"}, "FUV"),
    ({DETECTOR_TYPE_KEY: "NUV"}, "NUV"),
    ({DETECTOR_TYPE_KEY: "SJI"}], "SJI")
])
def test_get_detector_type(test_input, expected_output):
    assert iris_tools.get_detector_type(test_input) == expected_output

@pytest.mark.parametrize("data_arrays, old_unit, new_unit, expected_data_arrays, expected_unit", [
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["FUV"], u.photon,
     [SOURCE_DATA_PHOTONS_FUV, SOURCE_DATA_PHOTONS_FUV], u.photon),
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"], u.photon,
     [SOURCE_DATA_PHOTONS_NUV, SOURCE_DATA_PHOTONS_NUV], u.photon),
    ([SOURCE_DATA_DN_1, SOURCE_DATA_DN_1], iris_tools.DN_UNIT["SJI"], u.photon,
     [SOURCE_DATA_PHOTONS_SJI_1, SOURCE_DATA_PHOTONS_SJI_1], u.photon),
    ([SOURCE_DATA_PHOTONS_FUV, SOURCE_DATA_PHOTONS_FUV], u.photon, iris_tools.DN_UNIT["FUV"],
    [SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["FUV"]),
    ([SOURCE_DATA_PHOTONS_NUV, SOURCE_DATA_PHOTONS_NUV], u.photon, iris_tools.DN_UNIT["NUV"],
    [SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"]),
    ([SOURCE_DATA_PHOTONS_SJI_1, SOURCE_DATA_PHOTONS_SJI_1], u.photon, iris_tools.DN_UNIT["SJI"],
    [SOURCE_DATA_DN_1, SOURCE_DATA_DN_1], iris_tools.DN_UNIT["SJI"])
])
def test_convert_between_DN_and_photons(input_arrays, old_unit, new_unit,
                                        expected_arrays, expected_unit):
    output_arrays, output_unit = iris_tools.convert_between_DN_and_photons(input_arrays,
                                                                           old_unit, new_unit)
    for i, output_array in enumerate(output_arrays):
        np_test.assert_allclose(output_array, expected_arrays[i])
    assert output_unit == expected_unit

@pytest.mark.parametrize(
    "input_arrays, old_unit, exposure_time, force, expected_arrays, expected_unit",[
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon, EXPOSURE_TIME, False
         [SOURCE_DATA_DN/single_exposure_time, SOURCE_DATA_DN/single_exposure_time],
         u.photon/u.s),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"], EXPOSURE_TIME, False
         [SOURCE_DATA_DN/single_exposure_time, SOURCE_DATA_DN/single_exposure_time],
         iris_tools.DN_UNIT["NUV"]/u.s),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon/u.s, EXPOSURE_TIME, True,
         [SOURCE_DATA_DN/single_exposure_time, SOURCE_DATA_DN/single_exposure_time],
         u.photon/u.s/u.s),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"]/u.s, EXPOSURE_TIME, True,
         [SOURCE_DATA_DN/single_exposure_time, SOURCE_DATA_DN/single_exposure_time],
         iris_tools.DN_UNIT["NUV"]/u.s/u.s)
    ])
def test_calculate_exposure_time_correction(input_arrays, old_unit, exposure_time,
                                            expected_arrays, expected_unit):
    output_arrays, output_unit = iris_tools.calculate_exposure_time_correction(
        input_arrays, old_unit, exposure_time)
    for i, output_array in enumerate(output_arrays):
        np_test.assert_allclose(output_array, expected_arrays[i])
    assert output_unit == expected_unit

@pytest.mark.parametrize("input_arrays, old_unit, exposure_time, force", [
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon/u.s, EXPOSURE_TIME, False),
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"]/u.s, EXPOSURE_TIME, False)
])
def test_calculate_exposure_time_correction_error():
    assert pytest.raises(ValueError, iris_tools.calculate_exposure_time_correction,
                         input_arrays, old_unit, exposure_time, force)

@pytest.mark.parametrize(
    "input_arrays, old_unit, exposure_time, force, expected_arrays, expected_unit", [
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon/u.s, EXPOSURE_TIME, False,
         [SOURCE_DATA_DN * single_exposure_time, SOURCE_DATA_DN * single_exposure_time], u.photon),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"]/u.s, EXPOSURE_TIME, False,
         [SOURCE_DATA_DN * single_exposure_time, SOURCE_DATA_DN * single_exposure_time],
         iris_tools.DN_UNIT["NUV"]),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon, EXPOSURE_TIME, True,
         [SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon*u.s),
        ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["FUV"], EXPOSURE_TIME, True,
         [SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["FUV"]*u.s)
])
def test_uncalculate_exposure_time_correction(input_arrays, old_unit, exposure_time,
                                              expected_arrays, expected_unit):
        output_arrays, output_unit = iris_tools.uncalculate_exposure_time_correction(
            input_arrays, old_unit, exposure_time)
        for i, output_array in enumerate(output_arrays):
            np_test.assert_allclose(output_array, expected_arrays[i])
        assert output_unit == expected_unit

@pytest.mark.parametrize("input_arrays, old_unit, exposure_time, force", [
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], u.photon, EXPOSURE_TIME, False),
    ([SOURCE_DATA_DN, SOURCE_DATA_DN], iris_tools.DN_UNIT["NUV"], EXPOSURE_TIME, False)
])
def test_uncalculate_exposure_time_correction_error():
    assert pytest.raises(ValueError, iris_tools.uncalculate_exposure_time_correction,
                         input_arrays, old_unit, exposure_time, force)

def test_get_iris_response_response_version():
    assert pytest.raises(ValueError, iris_tools.get_iris_response, response_version=4)


def test_get_iris_response_not_equal_to_one():
    assert pytest.raises(ValueError, iris_tools.get_iris_response, pre_launch=True, response_version=3)


def test_get_iris_response_response_file():
    assert pytest.raises(KeyError, iris_tools.get_iris_response, response_file="hello.py")


# def test_get_iris_response():
# 	"""
# 	"""

# def test_gaussian1d_on_linear_bg():
# 	"""
# 	"""

# def test_calculate_orbital_wavelength_variation():
# 	"""
# 	"""