import os
import os.path as p
import time

import glob

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from builtins import range

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))

### General helpers

def array_relative_error(a1, a2):
    """Return the elementwise absolute difference between the inputs,
scaled by the maximum value that occurs in the input."""
    denom = max(np.amax(np.absolute(a1)), np.amax(np.absolute(a2)))
    difference = np.amax(np.absolute(a1 - a2))
    if denom < 1e-15:
        # Both input arrays are all zeros, so there is no relative error.
        return 0
    else:
        if difference < 1e-15:
            # absolute difference between arrays is ~ machine precision
            return 0
        else:
            return difference / denom

def assert_outputs_close(nx, ny, layers, rtol):

    # collct all files with a number in them (i.e. all outputs
    # other than the diagnostic files)
    outfiles = sorted(glob.glob("output/*.0*"))
    good_outfiles = sorted(glob.glob("good-output/*.0*"))

    # assert p.basename(outfiles) == p.basename(good_outfiles)
    test_passes = True
    for i, outfile in enumerate(outfiles):
        ans = aro.interpret_raw_file(outfile, nx, ny, layers)
        good_ans = aro.interpret_raw_file(good_outfiles[i], nx, ny, layers)
        relerr = np.amax(array_relative_error(ans, good_ans))
        if (relerr >= rtol or np.isnan(relerr)):
            print('test failed at {0}'.format(outfile))
            print('rtol = {0}: relerr = {1}'.format(rtol, relerr))
            # print ans
            # print good_ans
            test_passes = False

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.pcolormesh(ans[0,:,:])
            plt.colorbar()
            plt.title(outfile)
            ax.set_aspect('equal', 'datalim')
            plt.savefig('current_output_{0}.png'.format(outfile[12:]))
            plt.close()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.pcolormesh(good_ans[0,:,:])
            plt.colorbar()
            plt.title(outfile)
            ax.set_aspect('equal', 'datalim')
            plt.savefig('blessed_output_{0}.png'.format(outfile[12:]))
            plt.close()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.pcolormesh(ans[0,:,:] - good_ans[0,:,:], cmap='RdBu_r')
            plt.colorbar()
            plt.title('current - blessed at {0}'.format(outfile[12:]))
            ax.set_aspect('equal', 'datalim')
            plt.savefig('difference_{0}.png'.format(outfile[12:]))
            plt.close()

    assert test_passes

def assert_volume_conservation(nx,ny,layers,rtol):
    hfiles = sorted(glob.glob("output/snap.h.*"))

    h_0 = aro.interpret_raw_file(hfiles[0], nx, ny, layers)
    h_final = aro.interpret_raw_file(hfiles[-1], nx, ny, layers)

    volume_0 = np.zeros((layers))
    volume_final = np.zeros((layers))

    for k in range(layers):
        volume_0[k] = np.sum(h_0[k,:,:])
        volume_final[k] = np.sum(h_final[k,:,:])

        assert np.abs((volume_0[k] - volume_final[k])/volume_0[k]) < rtol

def assert_diagnostics_similar(variables, rtol):

    diags = {variable:np.loadtxt('output/diagnostic.{}.csv'.format(variable),
                delimiter=',', skiprows=1, usecols=(2,3,4)) for variable in variables}

    diags_blessed = {variable:np.loadtxt('good-output/diagnostic.{}.csv'.format(
        variable), delimiter=',', skiprows=1, usecols=(2,3,4)) for variable in variables}

    for variable in variables:
        # in case the diagnostic value is zero, add a very small number to it in
        # the denominator.
        assert np.max((diags[variable] - diags_blessed[variable])/
            (diags[variable]+1e-10)) < rtol
    

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
        assert_diagnostics_similar(['h'], 1e-10)

def test_f_plane():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "f_plane")):
        drv.simulate(exe=test_executable,
            nx=10, ny=10, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 1e-15)
        assert_volume_conservation(10, 10, 2, 1e-5)
        assert_diagnostics_similar(['h'], 1e-10)

def bump(X, Y):
    return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

def test_gaussian_bump_red_grav():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        drv.simulate(init_h_file=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)
        # assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_gaussian_bump_red_grav_continuation_of_single_core_MPI_2X():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        drv.simulate(init_h_file=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10,
                     niter0=101, n_time_steps=400, n_proc_x=2)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)


def test_gaussian_bump_red_grav_MPI_2X():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        drv.simulate(init_h_file=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10,
                     n_proc_x=2)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)
        # assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_gaussian_bump_red_grav_continuation_MPI_2X():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        drv.simulate(init_h_file=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10,
                     niter0=101, n_time_steps=400, n_proc_x=2)
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)


def test_gaussian_bump():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump")):
        drv.simulate(init_h_file=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     nx=10, ny=10, exe=test_executable, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)

def test_gaussian_bump_continuation():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump")):
        drv.simulate(init_h_file=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     nx=10, ny=10, exe=test_executable, 
                     dx=xlen/10, dy=ylen/10,
                     niter0=201, n_time_steps=200)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v', 'eta'], 1e-10)

def test_gaussian_bump_debug_test():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_debug_test")):
        drv.simulate(init_h_file=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     nx=10, ny=10, exe=test_executable, dx=xlen/10, dy=ylen/10)
        assert_outputs_close(10, 10, 2, 2e-13)
        assert_volume_conservation(10, 10, 2, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v', 'eta'], 1e-10)

def test_gaussian_bump_red_grav_debug_test():
    xlen = 1e6
    ylen = 1e6
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav_debug_test")):
        drv.simulate(init_h_file=[bump], exe=test_executable,
                     nx=10, ny=10, dx=xlen/10, dy=ylen/10,
                     n_proc_y=2
                     )
        assert_outputs_close(10, 10, 1, 1.5e-13)
        assert_volume_conservation(10, 10, 1, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_beta_plane_gyre_red_grav():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 20
    layers = 1
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
        drv.simulate(zonal_wind_file=[wind], valgrind=False,
                     nx=nx, ny=ny, exe=test_executable, dx=xlen/nx, dy=ylen/ny,
                     n_proc_x=2)
        assert_outputs_close(nx, ny, layers, 4e-13)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_beta_plane_gyre():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 10
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre")):
        drv.simulate(zonal_wind_file=[wind], valgrind=False,
                     nx=nx, ny=ny, exe="aronnax_test", dx=xlen/nx, dy=ylen/ny)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v', 'eta'], 1e-10)

def test_beta_plane_gyre_free_surf():
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 20
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)
    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))
    with working_directory(p.join(self_path, "beta_plane_gyre_free_surf")):
        drv.simulate(zonal_wind_file=[wind], valgrind=False,
                     nx=nx, ny=ny, exe=test_executable, dx=xlen/nx, dy=ylen/ny)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v', 'eta'], 1e-8)

def test_beta_plane_gyre_free_surf_Hypre_MPI(n_proc_x=1, n_proc_y=1):
    xlen = 1e6
    ylen = 2e6
    nx = 10; ny = 20
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

    with working_directory(p.join(self_path, "beta_plane_gyre_free_surf_hypre")):
        if (n_proc_x == 1 and n_proc_y == 1):
            drv.simulate(zonal_wind_file=[wind], valgrind=False,
                         nx=nx, ny=ny, exe='aronnax_external_solver_test',
                         dx=xlen/nx, dy=ylen/ny)
        elif (n_proc_x != 1 and n_proc_y == 1):
            drv.simulate(zonal_wind_file=[wind], valgrind=False,
                         nx=nx, ny=ny, exe='aronnax_external_solver_test',
                         dx=xlen/nx, dy=ylen/ny,
                         n_proc_x=n_proc_x)
        elif (n_proc_x == 1 and n_proc_y != 1):
            drv.simulate(zonal_wind_file=[wind], valgrind=False,
                         nx=nx, ny=ny, exe='aronnax_external_solver_test',
                         dx=xlen/nx, dy=ylen/ny,
                         n_proc_y=n_proc_y)
        else:
            drv.simulate(zonal_wind_file=[wind], valgrind=False,
                         nx=nx, ny=ny, exe='aronnax_external_solver_test',
                         dx=xlen/nx, dy=ylen/ny,
                         n_proc_x=n_proc_x, n_proc_y=n_proc_y)

        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v', 'eta'], 1e-8)

def test_beta_plane_gyre_free_surf_Hypre_MPI_2X():
    test_beta_plane_gyre_free_surf_Hypre_MPI(n_proc_x=2)

def test_beta_plane_gyre_free_surf_Hypre_MPI_2Y():
    test_beta_plane_gyre_free_surf_Hypre_MPI(n_proc_y=2)

def test_beta_plane_gyre_free_surf_Hypre_MPI_2X_2Y():
    test_beta_plane_gyre_free_surf_Hypre_MPI(n_proc_x=2, n_proc_y=2)


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
        drv.simulate(init_h_file=[layer_1, layer_2],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     exe=test_executable, wet_mask_file=[wetmask], 
                     f_u_file=[-1e-4],
                     f_v_file=[-1e-4],
                     n_time_steps=801,
                     dump_freq=10000, n_proc_y=2)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)


def test_relative_wind():
    nx = 320
    ny = 320
    layers = 1
    xlen = 1280e3
    ylen = 1280e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    def wetmask(X, Y):
        # water everywhere, doubly periodic
        wetmask = np.ones(X.shape, dtype=np.float64)
        return wetmask


    def wind_x(X, Y):
        rMx=300e3     # radius of circle where Max wind stress
        r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
        theta = np.arctan2(Y-640e3, X-640e3)
        tau = np.sin(np.pi*r/rMx/2.)
        tau_x = tau*np.sin(theta)
        tau_x[r>2.*rMx] = 0

        return tau_x

    def wind_y(X, Y):
        rMx=300e3     # radius of circle where Max wind stress
        r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
        theta = np.arctan2(Y-640e3, X-640e3)
        tau = np.sin(np.pi*r/rMx/2.)
        tau_y = -tau*np.cos(theta)
        tau_y[r>2.*rMx] = 0
        return tau_y

    with working_directory(p.join(self_path, 
        "relative_wind")):
        drv.simulate(
                init_h_file=[400],
                zonal_wind_file=[wind_x], meridional_wind_file=[wind_y],
                wind_mag_time_series_file=[0.08],
                exe=test_executable, wet_mask_file=[wetmask],
                nx=nx, ny=ny, dx=dx, dy=dy)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_relative_wind_upwind_advection():
    nx = 320
    ny = 320
    layers = 1
    xlen = 1280e3
    ylen = 1280e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    def wetmask(X, Y):
        # water everywhere, doubly periodic
        wetmask = np.ones(X.shape, dtype=np.float64)
        return wetmask


    def wind_x(X, Y):
        rMx=300e3     # radius of circle where Max wind stress
        r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
        theta = np.arctan2(Y-640e3, X-640e3)
        tau = np.sin(np.pi*r/rMx/2.)
        tau_x = tau*np.sin(theta)
        tau_x[r>2.*rMx] = 0

        return tau_x

    def wind_y(X, Y):
        rMx=300e3     # radius of circle where Max wind stress
        r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
        theta = np.arctan2(Y-640e3, X-640e3)
        tau = np.sin(np.pi*r/rMx/2.)
        tau_y = -tau*np.cos(theta)
        tau_y[r>2.*rMx] = 0
        return tau_y

    with working_directory(p.join(self_path, 
        "relative_wind")):
        drv.simulate(
                init_h_file=[400],
                zonal_wind_file=[wind_x], meridional_wind_file=[wind_y],
                wind_mag_time_series_file=[0.08],
                exe=test_executable, wet_mask_file=[wetmask],
                nx=nx, ny=ny, dx=dx, dy=dy, h_advec_scheme=2)
        assert_outputs_close(nx, ny, layers, 3e-11)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_periodic_BC_Hypre():
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
        drv.simulate(init_h_file=[500.,500.],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     exe="aronnax_external_solver_test",
                     wet_mask_file=[wetmask],
                     f_u_file=[-1e-4],
                     f_v_file=[-1e-4],
                     init_eta_file=[eta],
                     depth_file=[1000.], n_time_steps=801,
                     dump_freq=10000)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 3e-5)

def test_vertical_thickness_diffusion():
    nx = 2
    ny = 2

    dx = 1e4
    dy = 1e4

    grid = aro.Grid(nx, ny, 1, dx, dy)

    kv = 1e-5
    initH = 10.

    def wetmask(X, Y):
        mask = np.ones(X.shape, dtype=np.float64)
        return mask

    # 2 years
    dt = 200.
    n_time_steps = int(0.25*365*86400/dt)
    diag_freq = n_time_steps*dt/50.


    with working_directory(p.join(self_path, "vertical_diffusion")):
        drv.simulate(init_h_file=[initH],
                     nx=nx, ny=ny, dx=dx, dy=dy,
                     exe="aronnax_test",
                     wet_mask_file=[wetmask],
                     f_u_file=[-1e-4],
                     f_v_file=[-1e-4],
                     kv=kv,
                     dt=dt,
                     n_time_steps=n_time_steps,
                     diag_freq=diag_freq,
                     dump_freq=1e9)

        simuated_h_evo = np.loadtxt('output/diagnostic.h.csv',
                delimiter=',', skiprows=1, usecols=(0,2))
        expected_h_evo = np.sqrt(simuated_h_evo[:,0]*dt*kv*2. + initH**2)
        assert np.max(abs(simuated_h_evo[:,1] - expected_h_evo)) < 0.0005


def test_vertical_thickness_diffusion_Hypre_3_layers():
    nx = 2
    ny = 2
    layers = 3

    dx = 1e4
    dy = 1e4

    grid = aro.Grid(nx, ny, layers, dx, dy)

    kv = 1e-5
    initH = [10., 100., 10.]

    def wetmask(X, Y):
        mask = np.ones(X.shape, dtype=np.float64)
        return mask

    # 2 years
    dt = 200.
    n_time_steps = int(0.25*365*86400/dt)
    diag_freq = n_time_steps*dt/50.


    with working_directory(p.join(self_path, "vertical_diffusion")):
        drv.simulate(init_h_file=[10., 100., 10.],
                     nx=nx, ny=ny, dx=dx, dy=dy,
                     layers=3,
                     exe="aronnax_external_solver_test",
                     wet_mask_file=[wetmask],
                     depth_file=[120.],
                     f_u_file=[-1e-4],
                     f_v_file=[-1e-4],
                     kv=kv,
                     dt=dt,
                     n_time_steps=n_time_steps,
                     diag_freq=diag_freq,
                     dump_freq=1e9,
                     red_grav=0)

        simuated_h_evo = np.loadtxt('output/diagnostic.h.csv',
                delimiter=',', skiprows=1, usecols=(0,2,6,10))
        expected_h_evo = np.sqrt(simuated_h_evo[:,0]*dt*kv*2. + initH[0]**2)

        # plt.plot(simuated_h_evo[:,0]*dt, expected_h_evo, label='expected')
        # plt.plot(simuated_h_evo[:,0]*dt, simuated_h_evo[:,1], label='simulated top')
        # plt.plot(simuated_h_evo[:,0]*dt, simuated_h_evo[:,2], label='simulated mid')
        # plt.plot(simuated_h_evo[:,0]*dt, simuated_h_evo[:,3], label='simulated bottom')
        # plt.legend()
        # plt.savefig('h_evo.pdf')
        assert np.max(abs(simuated_h_evo[:,1] - simuated_h_evo[:,3])) < 0.0005


def test_outcropping_Hypre():

    nx = 60
    ny = 30
    layers = 2
    xlen = 120e3
    ylen = 120e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    def wetmask(X, Y):
        # start with water everywhere and add some land
        wetmask = np.ones(X.shape, dtype=np.float64)

        # clean up the edges
        wetmask[0,:] = 0
        wetmask[-1,:] = 0

        return wetmask

    def bathymetry(X,Y):
        depth = 50 + 50*Y/Y.max()
        return depth

    def init_h1(X, Y):
        depth = bathymetry(X,Y)
        proto_h1 = 150 - 200*Y/Y.max()
        # leave 0.1 m for bottom layer
        h1 = np.minimum(proto_h1, depth) - 0.1
        # set outrcropped region to 0.1 m thick
        h1 = np.maximum(h1, 0.1 + 0*X)
        # crude smoothing to reduce the sharp gradients
        h1[1:,:] = (h1[1:,:] + h1[:-1,:])/2
        h1[:,1:] = (h1[:,1:] + h1[:,:-1])/2.
        return h1

    def init_h2(X,Y):
        depth = bathymetry(X,Y)
        h1 = init_h1(X,Y)
        h2 = depth - h1
        return h2


    with working_directory(p.join(self_path, 'outcropping')):
        drv.simulate(
                init_h_file=[init_h1, init_h2],
                zonal_wind_file=[0.1], meridional_wind_file=[0],
                wind_depth=40,
                wet_mask_file=[wetmask],
                depth_file=[bathymetry],
                nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                exe='aronnax_external_solver_test',
                dt=10, n_time_steps=1000)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 3e-5)

def test_Hypre_spatial_wind():
    nx = 16
    ny = 32
    layers = 2
    xlen = 640e3
    ylen = 1280e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    wind_n_records = 2

    def wetmask(X, Y):
        # water everywhere, doubly periodic
        wetmask = np.ones(X.shape, dtype=np.float64)
        wetmask[0,:]  = 0
        wetmask[-1,:] = 0
        # wetmask[:,0] = 0
        # wetmask[:,-1] = 0
        return wetmask


    def wind_x(X, Y, wind_n_records):
        tau_x = np.zeros((wind_n_records, X.shape[0], X.shape[1]))

        rMx=300e3     # radius of circle where Max wind stress

        for i in range(wind_n_records):
            X0 = 640e3*i/(wind_n_records-1)
            r = np.sqrt((Y-640e3)**2 + (X-X0)**2)
            theta = np.arctan2(Y-640e3, X-X0)
            tau = np.sin(np.pi*r/rMx/2.)
            tau_x[i,:,:] = tau*np.sin(theta)
            r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
            tau_x[i, r>2.*rMx] = 0

        return tau_x

    def wind_y(X, Y, wind_n_records):
        tau_y = np.zeros((wind_n_records, X.shape[0], X.shape[1]))
        rMx=300e3     # radius of circle where Max wind stress

        for i in range(wind_n_records):
            X0 = 640e3*i/(wind_n_records-1)
            r = np.sqrt((Y-640e3)**2 + (X-X0)**2)
            theta = np.arctan2(Y-640e3, X-X0)
            tau = np.sin(np.pi*r/rMx/2.)
            tau_y[i,:,:] = -tau*np.cos(theta)
            r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
            tau_y[i, r>2.*rMx] = 0
        return tau_y

    with working_directory(p.join(self_path,
        "shifting_wind")):
        drv.simulate(
                init_h_file=[400, 400],
                depth_file=[800],
                wind_n_records=3,
                zonal_wind_file=[wind_x], meridional_wind_file=[wind_y],
                wind_mag_time_series_file=[0.1],
                exe=test_executable, wet_mask_file=[wetmask],
                nx=nx, ny=ny, dx=dx, dy=dy)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)

def test_Hypre_spatial_wind_interp():
    nx = 16
    ny = 32
    layers = 2
    xlen = 640e3
    ylen = 1280e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    wind_n_records = 2

    def wetmask(X, Y):
        # water everywhere, doubly periodic
        wetmask = np.ones(X.shape, dtype=np.float64)
        wetmask[0,:]  = 0
        wetmask[-1,:] = 0
        # wetmask[:,0] = 0
        # wetmask[:,-1] = 0
        return wetmask


    def wind_x(X, Y, wind_n_records):
        tau_x = np.zeros((wind_n_records, X.shape[0], X.shape[1]))

        rMx=300e3     # radius of circle where Max wind stress

        for i in range(wind_n_records):
            X0 = 640e3*i/(wind_n_records-1)
            r = np.sqrt((Y-640e3)**2 + (X-X0)**2)
            theta = np.arctan2(Y-640e3, X-X0)
            tau = np.sin(np.pi*r/rMx/2.)
            tau_x[i,:,:] = tau*np.sin(theta)
            r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
            tau_x[i, r>2.*rMx] = 0

        return tau_x

    def wind_y(X, Y, wind_n_records):
        tau_y = np.zeros((wind_n_records, X.shape[0], X.shape[1]))
        rMx=300e3     # radius of circle where Max wind stress

        for i in range(wind_n_records):
            X0 = 640e3*i/(wind_n_records-1)
            r = np.sqrt((Y-640e3)**2 + (X-X0)**2)
            theta = np.arctan2(Y-640e3, X-X0)
            tau = np.sin(np.pi*r/rMx/2.)
            tau_y[i,:,:] = -tau*np.cos(theta)
            r = np.sqrt((Y-640e3)**2 + (X-640e3)**2)
            tau_y[i, r>2.*rMx] = 0
        return tau_y

    with working_directory(p.join(self_path,
        "shifting_wind_interp")):
        drv.simulate(
                init_h_file=[400, 400],
                depth_file=[800],
                wind_n_records=3,
                zonal_wind_file=[wind_x], meridional_wind_file=[wind_y],
                wind_mag_time_series_file=[0.1],
                exe=test_executable, wet_mask_file=[wetmask],
                nx=nx, ny=ny, dx=dx, dy=dy)
        assert_outputs_close(nx, ny, layers, 3e-12)
        assert_volume_conservation(nx, ny, layers, 1e-5)
        assert_diagnostics_similar(['h', 'u', 'v'], 1e-10)
