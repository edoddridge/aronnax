'''Unit tests for Aronnax'''

from contextlib import contextmanager
import os.path as p
import re

import numpy as np
from scipy.io import FortranFile

import aronnax as aro

import pytest

def test_time_series_generator_func():
    '''Test the time series variable generator with a function.'''
    
    # choose dt and n_time_steps to give 360 days of time series
    n_time_steps = 360
    dt = 86400

    locally_generated_ts = ts_func(n_time_steps, dt)
    generator_function_ts = aro.time_series_variable(n_time_steps, dt, [ts_func])

    np.testing.assert_array_equal(locally_generated_ts, generator_function_ts)

    # expect the time series generator to raise and exception error if
    # the list `func` has more than one element.
    with pytest.raises(AssertionError):
        generator_function_ts = aro.time_series_variable(n_time_steps, dt, 
                                [ts_func, ts_func])

def ts_func(n_time_steps, dt):
    time_vector = np.arange(n_time_steps, dtype = np.float64)*dt

    # for the final fifteen days of each 30 day month set time series to zero
    # otherwise set it to one
    time_series = np.ones(n_time_steps, dtype = np.float64)
    time_series[(np.mod(time_vector, 30.)>15)] = 0    

    return time_series


def test_time_series_generator_const():
    '''Test the time series variable generator with a constant.'''
    
    # choose dt and n_time_steps to give 360 days of time series
    n_time_steps = 360
    dt = 86400
    const_value = 5.443

    locally_generated_ts = np.ones(n_time_steps, dtype=np.float64)*const_value
    generator_function_ts = aro.time_series_variable(n_time_steps, dt, [const_value])

    np.testing.assert_array_equal(locally_generated_ts, generator_function_ts)

    # expect the time series generator to raise and exception error if
    # the list `func` has more than one element.
    with pytest.raises(AssertionError):
        generator_function_ts = aro.time_series_variable(n_time_steps, dt, 
                                [const_value, const_value])




