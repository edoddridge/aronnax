from contextlib import contextmanager
import os
import os.path as p
import subprocess as sub
import time

import glob

import numpy as np
from scipy.io import FortranFile

import aronnax as aro

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

def run_experiment(write_input, nx, ny, layers, aro_exec=None, valgrind=False):
    if aro_exec is None:
        aro_exec = "aronnax_test"
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])
    with working_directory(root_path):
        sub.check_call(["make", aro_exec])
    with working_directory("input"):
        write_input(nx, ny, layers)
    tweak_parameters(nx, ny, layers)
    then = time.time()
    if valgrind or 'ARONNAX_TEST_VALGRIND_ALL' in os.environ:
        sub.check_call(["valgrind", "--error-exitcode=5", p.join(root_path, aro_exec)])
    else:
        sub.check_call([p.join(root_path, aro_exec)])
    run_time = time.time() - then
    print "Aronnax execution took", run_time
    return run_time

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
        ans = aro.interpret_raw_file(p.join("output/", outfile), nx, ny, layers)
        good_ans = aro.interpret_raw_file(p.join("good-output/", outfile), nx, ny, layers)
        assert np.amax(array_relative_error(ans, good_ans)) < rtol

def assert_volume_conservation(nx,ny,layers,rtol):
    hfiles = sorted(glob.glob("output/snap.h.*"))

    h_0 = aro.interpret_raw_file(hfiles[0], nx, ny, layers)
    h_final = aro.interpret_raw_file(hfiles[-1], nx, ny, layers)

    volume_0 = np.zeros((layers))
    volume_final = np.zeros((layers))

    for k in xrange(layers):
        volume_0[k] = np.sum(h_0[:,:,k])
        volume_final[k] = np.sum(h_final[:,:,k])

        assert np.abs((volume_0[k] - volume_final[k])/volume_0[k]) < rtol


### The test cases themselves

def write_input_f_plane_red_grav(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 1e6
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)
    aro.write_f_plane(nx, ny, 10e-4)
    aro.write_rectangular_pool(nx, ny)
    aro.write_initial_heights(grid, [400.0])

def test_f_plane_red_grav():
    with working_directory(p.join(self_path, "f_plane_red_grav")):
        run_experiment(write_input_f_plane_red_grav, 10, 10, 1)
        assert_outputs_close(10, 10, 1, 1e-15)
        assert_volume_conservation(10, 10, 1, 1e-5)

def write_input_f_plane(nx, ny, layers):
    assert layers == 2
    xlen = 1e6
    ylen = 1e6
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)
    aro.write_f_plane(nx, ny, 10e-4)
    aro.write_rectangular_pool(nx, ny)
    aro.write_initial_heights(grid, [400.0, 1600.0])

def test_f_plane():
    with working_directory(p.join(self_path, "f_plane")):
        run_experiment(write_input_f_plane, 10, 10, 2)
        assert_outputs_close(10, 10, 2, 1e-15)
        assert_volume_conservation(10, 10, 2, 1e-5)

def write_input_beta_plane_bump_red_grav(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 1e6
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)

    aro.write_beta_plane(grid, 1e-5, 2e-11)
    aro.write_rectangular_pool(nx, ny)
    def bump(X, Y):
        return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))
    aro.write_initial_heights(grid, [bump])

def test_gaussian_bump_red_grav():
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        run_experiment(write_input_beta_plane_bump_red_grav, 10, 10, 1)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)

def write_input_beta_plane_bump(nx, ny, layers):
    assert layers == 2
    xlen = 1e6
    ylen = 1e6
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)

    aro.write_beta_plane(grid, 1e-5, 2e-11)
    aro.write_rectangular_pool(nx, ny)
    def bump(X, Y):
        return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))
    aro.write_initial_heights(grid, [bump, lambda X, Y: 2000. - bump(X, Y)])

def test_gaussian_bump():
    with working_directory(p.join(self_path, "beta_plane_bump")):
        run_experiment(write_input_beta_plane_bump, 10, 10, 2)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)

def write_input_beta_plane_gyre_red_grav(nx, ny, layers):
    assert layers == 1
    xlen = 1e6
    ylen = 2e6
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)

    aro.write_beta_plane(grid, 1e-5, 2e-11)
    aro.write_rectangular_pool(nx, ny)
    aro.write_initial_heights(grid, [400.0])

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
    grid = aro.Grid(nx, ny, xlen / nx, ylen / ny)

    aro.write_beta_plane(grid, 1e-5, 2e-11)
    aro.write_rectangular_pool(nx, ny)
    aro.write_initial_heights(grid, [600.0, 1400.0])

    with fortran_file('wind_x.bin', 'w') as f:
        _, Y = np.meshgrid(grid.xp1, grid.y)
        wind_x = 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
        f.write_record(wind_x.astype(np.float64))

def test_beta_plane_gyre():
    with working_directory(p.join(self_path, "beta_plane_gyre")):
        run_experiment(write_input_beta_plane_gyre, 10, 10, 2, valgrind=True)
        assert_outputs_close(10, 10, 2, 3e-12)
        assert_volume_conservation(10, 10, 2, 1e-5)
