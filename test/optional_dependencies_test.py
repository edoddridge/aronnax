'''Tests that the optional analysis dependencies are only needed on use.'''

import subprocess as sub
import sys

SCRIPT = """
import sys
# A None entry in sys.modules makes the import of that module fail.
sys.modules['dask'] = None
sys.modules['dask.array'] = None
sys.modules['xarray'] = None

import aronnax as aro

grid = aro.Grid(10, 10, 1, 1e3, 1e3)
try:
    aro.open_mfdataarray(['snap.h.0000000001'], grid)
except ImportError as err:
    assert 'dask and xarray' in str(err), err
else:
    raise SystemExit('expected an ImportError')
"""


def test_import_without_dask_and_xarray():
    '''Importing aronnax must not require dask or xarray, and loading
    output without them must give an informative error.'''
    sub.check_call([sys.executable, '-c', SCRIPT])
