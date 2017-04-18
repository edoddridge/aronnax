"""!
The main file with the class definitions

Core
==============

This file contains all of the classes for the module.
"""

from contextlib import contextmanager
import os.path as p
import re

import numpy as np
from scipy.io import FortranFile

class Grid(object):
    """Make a grid object containing all of the axes.

        :param int nx: Number of grid points in the x direction
        :param int ny: Number of grid points in the y direction
        :param float dx: Grid size in x direction in metres
        :param float dy: Grid size in y direction in metres
        :param float x0: x value at lower left corner of domain
        :param float y0: y value at lower left corner of domain"""

    def __init__(self,nx,ny,dx,dy,x0=0,y0=0):
        """Instantiate a grid object for Aronnax."""

        # axes for vorticity points
        self.xp1 = np.linspace(x0,nx*dx+x0,nx+1)
        self.yp1 = np.linspace(y0,ny*dy+y0,ny+1)

        # Axes for tracer points.
        self.x = (self.xp1[1:] + self.xp1[:-1])/2.
        self.y = (self.yp1[1:] + self.yp1[:-1])/2.

        # Size
        self.nx = nx
        self.ny = ny

@contextmanager
def fortran_file(*args, **kwargs):
    f = FortranFile(*args, **kwargs)
    try:
        yield f
    finally:
        f.close()

def interpret_raw_file(name, nx, ny, layers):
    """Read an output file dumped by the Aronnax core.

    Each such file contains one array, whose size depends on what,
    exactly, is in it, and on the resolution of the simulation.
    Hence, the parameters nx, ny, and layers, as well as the file
    naming convetion, suffice to interpret the content (assuming it
    was generated on the same system)."""
    # Note: This depends on inspection of the output writing code in
    # the Aronnax core, to align array sizes and dimensions.  In
    # particular, Fortran arrays are indexed in decreasing order of
    # rate of change as one traverses the elements sequentially,
    # whereas Python (and all other programming languages I am aware
    # of) indexes in increasing order.
    file_part = p.basename(name)
    dx = 0; dy = 0; layered = True
    if file_part.startswith("snap.h"):
        pass
    if file_part.startswith("snap.u"):
        dx = 1
    if file_part.startswith("snap.v"):
        dy = 1
    if file_part.startswith("snap.eta"):
        layered = False
    if file_part.startswith("wind_x"):
        dx = 1
        layered = False
    if file_part.startswith("wind_y"):
        dy = 1
        layered = False
    if file_part.startswith("av.h"):
        pass
    if file_part.startswith("av.u"):
        dx = 1
    if file_part.startswith("av.v"):
        dy = 1
    if file_part.startswith("av.eta"):
        layered = False
    with fortran_file(name, 'r') as f:
        if layered:
            return f.read_reals(dtype=np.float64) \
                    .reshape(layers, ny+dy, nx+dx).transpose()
        else:
            return f.read_reals(dtype=np.float64) \
                    .reshape(ny+dy, nx+dx).transpose()

### General input construction helpers

def depths(grid, *h_funcs):
    X,Y = np.meshgrid(grid.x, grid.y)
    initH = np.ones((len(h_funcs), grid.ny, grid.nx))
    for i, f in enumerate(h_funcs):
        if isinstance(f, (int, long, float)):
            initH[i,:,:] = f
        else:
            initH[i,:,:] = f(X, Y)
    return initH

def wind_x(grid, func):
    X,Y = np.meshgrid(grid.xp1, grid.y)
    if isinstance(func, (int, long, float)):
        wind_x = np.ones(grid.ny, grid.nx+1) * func
    else:
        wind_x = func(X, Y)
    return wind_x

def wind_y(grid, func):
    X,Y = np.meshgrid(grid.y, grid.xp1)
    if isinstance(func, (int, long, float)):
        wind_y = np.ones(grid.ny+1, grid.nx) * func
    else:
        wind_y = func(X, Y)
    return wind_y

### Specific construction helpers

def f_plane_u(grid, coeff):
    """Define an f-plane approximation to the Coriolis force (u component)."""
    return np.ones((grid.nx+1, grid.ny), dtype=np.float64) * coeff

def f_plane_v(grid, coeff):
    """Define an f-plane approximation to the Coriolis force (v component)."""
    return np.ones((grid.nx, grid.ny+1), dtype=np.float64) * coeff

def beta_plane_u(grid, f0, beta):
    """Define a beta-plane approximation to the Coriolis force (u component)."""
    _, Y = np.meshgrid(grid.xp1, grid.y)
    fu = f0 + Y*beta
    return fu

def beta_plane_v(grid, f0, beta):
    """Define a beta-plane approximation to the Coriolis force (v component)."""
    _, Y = np.meshgrid(grid.x, grid.yp1)
    fv = f0 + Y*beta
    return fv

def rectangular_pool(grid):
    """The wet mask file for a maximal rectangular pool."""
    nx = grid.nx; ny = grid.ny
    wetmask = np.ones((nx, ny), dtype=np.float64)
    wetmask[ 0, :] = 0
    wetmask[-1, :] = 0
    wetmask[ :, 0] = 0
    wetmask[ :,-1] = 0
    return wetmask

specifier_rx = re.compile(r':(.*):(.*)')

ok_generators = {
    'depths': depths,
    'beta_plane_u': beta_plane_u,
    'beta_plane_v': beta_plane_v,
    'f_plane_u': f_plane_u,
    'f_plane_v': f_plane_v,
    'rectangular_pool': rectangular_pool,
    'wind_x': wind_x,
    'wind_y': wind_y,
}

def interpret_data_specifier(string):
    m = re.match(specifier_rx, string)
    if m:
        name = m.group(1)
        arg_str = m.group(2)
        if len(arg_str) > 0:
            args = [float(a) for a in arg_str.split(',')]
        else:
            args = []
        return (ok_generators[name], args)
    else:
        return None

def interpret_requested_data(requested_data, shape, config):
    """Interpret a flexible input data specification.

    The requested_data can be one of

    - TODO A string giving the path to a NetCDF file, whose content
      will be interpolated to match the desired grid specification;

    - A string giving the path to a raw Fortran array file, whose
      content will be used as-is;

    - TODO A numpy array in memory, whose content will be used as-is,
      or TODO interpolated; or

    - A specifier for auto-generating the required data, in this format:
      :<generator_func_name>:arg1,arg2,...argn
    """
    grid = Grid(config.getint("grid", "nx"), config.getint("grid", "ny"),
                config.getfloat("grid", "dx"), config.getfloat("grid", "dy"))
    if isinstance(requested_data, basestring):
        candidate = interpret_data_specifier(requested_data)
        if candidate is not None:
            (func, args) = candidate
            return func(grid, *args)
        else:
            # Assume Fortran file name
            with fortran_file(requested_data, 'r') as f:
                return f.read_reals(dtype=np.float64)
    else:
        if shape == "3d":
            return depths(grid, *requested_data)
        if shape == "2dx":
            return wind_x(grid, requested_data)
        else:
            raise Exception("TODO implement custom generation for other input shapes")
