'''Tests for the error handling of the output reading functions.'''

import aronnax as aro

import pytest


def test_interpret_raw_file_unrecognised_name():
    '''An unknown file name gives a clear error, not a TypeError from
    unpacking None.'''
    with pytest.raises(ValueError, match="Unrecognised output file name"):
        aro.interpret_raw_file('output/diagnostic.h.csv', 10, 10, 1)


def test_open_mfdataarray_unrecognised_name():
    grid = aro.Grid(10, 10, 1, 1e3, 1e3)
    with pytest.raises(ValueError, match="Unrecognised output file name"):
        aro.open_mfdataarray(['output/bogus.0000000001'], grid)


def test_open_mfdataarray_no_files():
    '''An empty file list (for example from a glob that matched nothing)
    gives a clear error, not an IndexError.'''
    grid = aro.Grid(10, 10, 1, 1e3, 1e3)
    with pytest.raises(ValueError, match="at least one file"):
        aro.open_mfdataarray([], grid)
