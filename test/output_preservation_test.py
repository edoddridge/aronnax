import os
import os.path as p
import subprocess as sub
import time

import glob

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))

### General helpers

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
        relerr = np.amax(array_relative_error(ans, good_ans))
        if (relerr >= rtol or np.isnan(relerr)):
            print outfile
            print ans
            print good_ans

            plt.figure()
            plt.pcolormesh(ans[:,:,0])
            plt.colorbar()
            plt.title(outfile)
            plt.savefig('current_output.png')
            plt.close()

            plt.figure()
            plt.pcolormesh(good_ans[:,:,0])
            plt.colorbar()
            plt.title(outfile)
            plt.savefig('blessed_output.png')
            plt.close()

        assert relerr < rtol

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

test_executable = "aronnax_test"

def test_f_plane_red_grav():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "f_plane_red_grav")):
        drv.simulate(exe=test_executable,
            nx=10, ny=10, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 1, 1e-15)
        assert_volume_conservation(10, 10, 1, 1e-5)

def test_f_plane():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "f_plane")):
        drv.simulate(exe=test_executable,
            nx=10, ny=10, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 1e-15)
        assert_volume_conservation(10, 10, 2, 1e-5)

def bump(X, Y):
    return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

def test_gaussian_bump_red_grav():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        drv.simulate(initHfile=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)

def test_gaussian_bump():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump")):
        drv.simulate(initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     nx=10, ny=10, exe=test_executable, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)

def test_gaussian_bump_debug_test():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_debug_test")):
        drv.simulate(initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     nx=10, ny=10, exe=test_executable, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)

def test_beta_plane_gyre_red_grav():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 20
    layers = 1
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
        drv.simulate(zonalWindFile=wind, valgrind=False,
                     nx=nx, ny=ny, exe=test_executable, dx=xlen/nx, dy=ylen/ny)
        assert_outputs_close(nx, ny, layers, 4e-13)
        assert_volume_conservation(nx, ny, layers, 1e-5)

def test_beta_plane_gyre():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 10
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre")):
        drv.simulate(zonalWindFile=wind, valgrind=False,
                     nx=nx, ny=ny, exe="aronnax_test", dx=xlen/nx, dy=ylen/ny)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)

def test_beta_plane_gyre_free_surf():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 20
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre_free_surf")):
        drv.simulate(zonalWindFile=wind, valgrind=False,
                     nx=nx, ny=ny, exe=test_executable, dx=xlen/nx, dy=ylen/ny)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)

def test_periodic_BC():
    nx = 50
    ny = 20
    layers = 2

    dx = 5e4
    dy = 5e4

    grid = aro.Grid(nx, ny, layers, dx, dy)

    rho0 = 1035

    def wetmask(X, Y):
        mask = np.ones(X.shape, dtype=np.float64)
        return mask

    def eta(X, Y):
        return 0. + 1.*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

    with working_directory(p.join(self_path, "periodic_BC")):
        drv.simulate(initHfile=[500.,500.],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     exe="aronnax_external_solver_test",
                     wetMaskFile=wetmask, 
                     fUfile=-1e-4,
                     fVfile=-1e-4,
                     initEtaFile=eta,
                     depthFile=1000., nTimeSteps=801,
                     dumpFreq=10000)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 3e-5)

def test_periodic_BC_red_grav():
    nx = 50
    ny = 20
    layers = 2

    dx = 5e4
    dy = 5e4

    grid = aro.Grid(nx, ny, layers, dx, dy)

    rho0 = 1035

    def wetmask(X, Y):
        mask = np.ones(X.shape, dtype=np.float64)
        return mask

    def layer_1(X, Y):
        return 500. + 300*np.exp(-((6e5-X)**2 + (6e5-Y)**2)/(2*1e5**2))

    def layer_2(X, Y):
        return 1000. - layer_1(X, Y)

    with working_directory(p.join(self_path, "periodic_BC_red_grav")):
        drv.simulate(initHfile=[layer_1, layer_2],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     exe=test_executable, wetMaskFile=wetmask, 
                     fUfile=-1e-4,
                     fVfile=-1e-4,
                     nTimeSteps=801,
                     dumpFreq=10000)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)