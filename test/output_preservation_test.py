from contextlib import contextmanager
import os
import os.path as p
import subprocess as sub

self_path = p.dirname(p.abspath(__file__))

import numpy as np
from scipy.io import FortranFile

import MIMutils as mim

@contextmanager
def fortran_file(*args, **kwargs):
    f = FortranFile(*args, **kwargs)
    try:
        yield f
    finally:
        f.close()

@contextmanager
def working_directory(path):
    old_path = os.getcwd()
    sub.check_call(["mkdir", "-p", path])
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_path)

def compile_mim(nx, ny, layers):
    mim_path = p.join(p.dirname(self_path), "MIM.f90")
    sub.check_call(
        "cat %s " % (mim_path,) +
        "| sed 's/^    integer, parameter :: nx =.*$/    integer, parameter :: nx = %d/'" % (nx,) +
        "| sed 's/^    integer, parameter :: ny =.*$/    integer, parameter :: ny = %d/'" % (ny,) +
        "| sed 's/^    integer, parameter :: layers =.*$/    integer, parameter :: layers = %d/'" % (layers,) +
        "> MIM.f90", shell=True)
    sub.check_call(["gfortran", "-Ofast", "MIM.f90", "-o", "MIM"])

def write_f_plane(nx, ny, coeff):
    """Write files defining an f-plane approximation to the Coriolis force."""
    with fortran_file('fu.bin', 'w') as f:
        f.write_record(np.ones((nx, ny+1), dtype=np.float64) * coeff)
    with fortran_file('fv.bin', 'w') as f:
        f.write_record(np.ones((nx+1, ny), dtype=np.float64) * coeff)

def write_rectangular_pool(nx, ny):
    """Write the wet mask file for a maximal rectangular pool."""
    with fortran_file('wetmask.bin', 'w') as f:
        wetmask = np.ones((nx, ny), dtype=np.float64)
        wetmask[ 0, :] = 0
        wetmask[-1, :] = 0
        wetmask[ :, 0] = 0
        wetmask[ :,-1] = 0
        f.write_record(wetmask)

def write_input_f_plane_red(nx, ny, layers):
    assert layers == 1
    write_f_plane(nx, ny, 10e-4)
    write_rectangular_pool(nx, ny)
    with fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny), dtype=np.float64)*400)

def run_experiment(write_input, nx, ny, layers):
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])
    with working_directory("input"):
        write_input(nx, ny, layers)
    compile_mim(nx, ny, layers)
    sub.check_call(["MIM"])

def test_f_plane_red():
    with working_directory(p.join(self_path, "f_plane_red")):
        run_experiment(write_input_f_plane_red, 10, 10, 1)
        sub.check_call(["diff", "-ru", "good-output/", "output/"])

def write_input_f_plane(nx, ny, layers):
    assert layers == 2
    write_f_plane(nx, ny, 10e-4)
    write_rectangular_pool(nx, ny)
    with fortran_file('initH.bin', 'w') as f:
        initH = np.ones((2,nx,ny), dtype=np.float64)
        initH[0,:,:] = 400
        initH[1,:,:] = 2000 - initH[0,:,:]
        f.write_record(initH)

def test_f_plane():
    with working_directory(p.join(self_path, "f_plane")):
        run_experiment(write_input_f_plane, 10, 10, 2)
        sub.check_call(["diff", "-ru", "good-output/", "output/"])

def write_input_beta_plane_bump_red(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 1e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    f0 = 1e-5
    beta = 2e-11

    with fortran_file('fu.bin', 'w') as f:
        X,Y = np.meshgrid(grid.xp1,grid.y)
        fu = f0 + Y*beta
        f.write_record(fu.astype(np.float64))

    with fortran_file('fv.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.yp1)
        fv = f0 + Y*beta
        f.write_record(fv.astype(np.float64))

    write_rectangular_pool(nx, ny)

    with fortran_file('initH.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        initH = 500 + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))
        f.write_record(initH.astype(np.float64))

def test_gaussian_bump_red():
    with working_directory(p.join(self_path, "beta_plane_bump_red")):
        run_experiment(write_input_beta_plane_bump_red, 10, 10, 1)
        sub.check_call(["diff", "-ru", "good-output/", "output/"])
