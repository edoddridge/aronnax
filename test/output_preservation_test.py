from contextlib import contextmanager
import os
import os.path as p
import subprocess as sub
import time

import glob

import numpy as np
from scipy.io import FortranFile

import aronnax as mim

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

### General helpers

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

def tweak_parameters(nx, ny, layers):
    sub.check_call(
        "cat parameters.in " +
        "| sed 's/^ nx =.*,$/ nx = %d,/'" % (nx,) +
        "| sed 's/^ ny =.*,$/ ny = %d,/'" % (ny,) +
        "| sed 's/^ layers =.*,$/ layers = %d,/'" % (layers,) +
        "> parameters.new", shell=True)
    sub.check_call(["mv", "parameters.new", "parameters.in"])

def run_experiment(write_input, nx, ny, layers, mim_exec=None, valgrind=False):
    if mim_exec is None:
        mim_exec = "MIM_test"
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])
    with working_directory(root_path):
        sub.check_call(["make", mim_exec])
    with working_directory("input"):
        write_input(nx, ny, layers)
    tweak_parameters(nx, ny, layers)
    then = time.time()
    if valgrind or 'MIM_TEST_VALGRIND_ALL' in os.environ:
        sub.check_call(["valgrind", "--error-exitcode=5", p.join(root_path, mim_exec)])
    else:
        sub.check_call([p.join(root_path, mim_exec)])
    run_time = time.time() - then
    print "MIM execution took", run_time
    return run_time

def interpret_mim_raw_file(name, nx, ny, layers):
    """Read an output file dumped by the MIM core.

    Each such file contains one array, whose size depends on what,
    exactly, is in it, and on the resolution of the simulation.
    Hence, the parameters nx, ny, and layers, as well as the file
    naming convetion, suffice to interpret the content (assuming it
    was generated on the same system)."""
    # Note: This depends on inspection of the output writing code in
    # the MIM core, to align array sizes and dimensions.  In
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

def array_relative_error(a1, a2):
    """Return the elementwise absolute difference between the inputs,
scaled by the maximum value that occurs in the input."""
    denom = max(np.amax(np.absolute(a1)), np.amax(np.absolute(a2)))
    if denom == 0:
        # Both input arrays are all zeros, so there is no relative error.
        return 0
    else:
        return np.absolute(a1 - a2) / denom

def assert_outputs_close(nx, ny, layers, rtol):
    outfiles = sorted(os.listdir("output/"))
    assert outfiles == sorted(os.listdir("good-output/"))
    for outfile in outfiles:
        ans = interpret_mim_raw_file(p.join("output/", outfile), nx, ny, layers)
        good_ans = interpret_mim_raw_file(p.join("good-output/", outfile), nx, ny, layers)
        assert np.amax(array_relative_error(ans, good_ans)) < rtol

def assert_volume_conservation(nx,ny,layers,rtol):
    hfiles = sorted(glob.glob("output/snap.h.*"))

    h_0 = interpret_mim_raw_file(hfiles[0], nx, ny, layers)
    h_final = interpret_mim_raw_file(hfiles[-1], nx, ny, layers)

    volume_0 = np.zeros((layers))
    volume_final = np.zeros((layers))

    for k in xrange(layers):
        volume_0[k] = np.sum(h_0[:,:,k])
        volume_final[k] = np.sum(h_final[:,:,k])

        assert np.abs((volume_0[k] - volume_final[k])/volume_0[k]) < rtol


### Input construction helpers

def write_f_plane(nx, ny, coeff):
    """Write files defining an f-plane approximation to the Coriolis force."""
    with fortran_file('fu.bin', 'w') as f:
        f.write_record(np.ones((nx+1, ny), dtype=np.float64) * coeff)
    with fortran_file('fv.bin', 'w') as f:
        f.write_record(np.ones((nx, ny+1), dtype=np.float64) * coeff)

def write_beta_plane(grid, f0, beta):
    """Write files defining a beta-plane approximation to the Coriolis force."""
    with fortran_file('fu.bin', 'w') as f:
        _, Y = np.meshgrid(grid.xp1, grid.y)
        fu = f0 + Y*beta
        f.write_record(fu.astype(np.float64))
    with fortran_file('fv.bin', 'w') as f:
        _, Y = np.meshgrid(grid.x, grid.yp1)
        fv = f0 + Y*beta
        f.write_record(fv.astype(np.float64))

def write_rectangular_pool(nx, ny):
    """Write the wet mask file for a maximal rectangular pool."""
    with fortran_file('wetmask.bin', 'w') as f:
        wetmask = np.ones((nx, ny), dtype=np.float64)
        wetmask[ 0, :] = 0
        wetmask[-1, :] = 0
        wetmask[ :, 0] = 0
        wetmask[ :,-1] = 0
        f.write_record(wetmask)

### The test cases themselves

def write_input_f_plane_red_grav(nx, ny, layers):
    assert layers == 1
    write_f_plane(nx, ny, 10e-4)
    write_rectangular_pool(nx, ny)
    with fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny), dtype=np.float64) * 400)

def test_f_plane_red_grav():
    with working_directory(p.join(self_path, "f_plane_red_grav")):
        run_experiment(write_input_f_plane_red_grav, 10, 10, 1)
        assert_outputs_close(10, 10, 1, 1e-15)
        assert_volume_conservation(10, 10, 1, 1e-5)

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
        assert_outputs_close(10, 10, 2, 1e-15)
        assert_volume_conservation(10, 10, 2, 1e-5)

def write_input_beta_plane_bump_red_grav(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 1e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    write_beta_plane(grid, 1e-5, 2e-11)
    write_rectangular_pool(nx, ny)

    with fortran_file('initH.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        initH = 500 + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))
        f.write_record(initH.astype(np.float64))

def test_gaussian_bump_red_grav():
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        run_experiment(write_input_beta_plane_bump_red_grav, 10, 10, 1)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)

def write_input_beta_plane_bump(nx, ny, layers):
    assert layers == 2
    xlen = 1e6
    ylen = 1e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    write_beta_plane(grid, 1e-5, 2e-11)
    write_rectangular_pool(nx, ny)

    with fortran_file('initH.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        initH = np.ones((2,ny,nx))
        initH[0,:,:] = 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))
        initH[1,:,:] = 2000. - initH[0,:,:]
        f.write_record(initH.astype(np.float64))

def test_gaussian_bump():
    with working_directory(p.join(self_path, "beta_plane_bump")):
        run_experiment(write_input_beta_plane_bump, 10, 10, 2)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)

def write_input_beta_plane_gyre_red_grav(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 2e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    write_beta_plane(grid, 1e-5, 2e-11)
    write_rectangular_pool(nx, ny)

    with fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny), dtype=np.float64) * 400)

    with fortran_file('wind_x.bin', 'w') as f:
        _, Y = np.meshgrid(grid.xp1, grid.y)
        wind_x = 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
        f.write_record(wind_x.astype(np.float64))

def test_beta_plane_gyre_red_grav():
    with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
        run_experiment(write_input_beta_plane_gyre_red_grav, 10, 10, 1, valgrind=True)
        assert_outputs_close(10, 10, 1, 2e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)

def write_input_beta_plane_gyre(nx, ny, layers):
    assert layers == 2
    xlen = 1e6
    ylen = 2e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    write_beta_plane(grid, 1e-5, 2e-11)
    write_rectangular_pool(nx, ny)

    with fortran_file('initH.bin', 'w') as f:
        _, Y = np.meshgrid(grid.x, grid.y)
        initH = np.ones((2, ny, nx))
        initH[0,:,:] = 600.
        initH[1,:,:] = 2000. - initH[0,:,:]
        f.write_record(initH.astype(np.float64))

    with fortran_file('wind_x.bin', 'w') as f:
        _, Y = np.meshgrid(grid.xp1, grid.y)
        wind_x = 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
        f.write_record(wind_x.astype(np.float64))

def test_beta_plane_gyre():
    with working_directory(p.join(self_path, "beta_plane_gyre")):
        run_experiment(write_input_beta_plane_gyre, 10, 10, 2, valgrind=True)
        assert_outputs_close(10, 10, 2, 3e-12)
        assert_volume_conservation(10, 10, 2, 1e-5)
