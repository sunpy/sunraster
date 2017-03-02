# -*- coding: utf-8 -*-
from __future__ import print_function
import pytest
import irispy.iris_tools as iris_tools
import numpy as np
import numpy.testing as np_test

source_data = np.array([[ 0.563,  1.132, -1.343],
						[-0.719,  1.441, 1.566]])

source_data1 = np.array([[1, 2, 3],
						[4, 5, 6]])

def test_convert_DN_to_photons_NUV():
	"""
	"""
	expected_output = np.array([[ 10.134,  20.376, -24.174],
       					[-12.942,  25.938,  28.188]])

	photons_count = iris_tools.convert_DN_to_photons(source_data, 'NUV')

	np_test.assert_allclose(photons_count, expected_output)


def test_convert_DN_to_photons_FUV():
	"""
	"""
	expected_output = np.array([[ 2.252,  4.528, -5.372],
       					[-2.876,  5.764,  6.264]])

	photons_count = iris_tools.convert_DN_to_photons(source_data, 'FUV')

	np_test.assert_allclose(photons_count, expected_output)


# def test_convert_DN_to_photons_SJI():
# 	"""
# 	"""


def test_convert_photons_to_DN_NUV():
	"""
	"""
	expected_output = np.array([[ 0.05555556,  0.11111111,  0.16666667],
    						[ 0.22222222,  0.27777778,  0.33333333]])

	DN = iris_tools.convert_photons_to_DN(source_data1, 'NUV')

	np_test.assert_allclose(DN, expected_output)


def test_convert_photons_to_DN_FUV():
	"""
	"""
	expected_output = np.array([[ 0.25,  0.5 ,  0.75],
       							[ 1.  ,  1.25,  1.5 ]])

	DN = iris_tools.convert_photons_to_DN(source_data1, 'FUV')

	np_test.assert_allclose(DN, expected_output)


# def test_convert_photons_to_DN_SJI():
# 	"""
# 	"""

# def test_calculate_intensity_fractional_uncertainty():
# 	"""
# 	"""

# def test_get_iris_response():
# 	"""
# 	"""

# def test_gaussian1d_on_linear_bg():
# 	"""
# 	"""

# def test_calculate_orbital_wavelength_variation():
# 	"""
# 	"""