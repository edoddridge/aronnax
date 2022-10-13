import os.path as p

import glob

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt

def test_f_plane_Hypre_bot_drag(bot_drag=1e-5, layers=1):

    test_executable = "aronnax_external_solver_test"

    dt = 100

    nx = 10
    ny = 10
    dx = 1e3
    dy = 1e3

    rho0 = 1035.

    grid = aro.Grid(nx, ny, layers, dx, dy)


    def init_U(X, Y, *arg):
        init_u = np.zeros(Y.shape,dtype=np.float64)
        init_u[:,:] = 3e-5

        return init_u

    def dbl_periodic_wetmask(X, Y):
        return np.ones(X.shape,dtype=np.float64)

    with working_directory(p.join(self_path, "bot_drag")):

        if layers ==1:
            drv.simulate(init_h_file=[400.], 
            init_u_file=[init_U],
            depth_file=[layers*400],
            exe=test_executable,
            wet_mask_file=[dbl_periodic_wetmask],
                     nx=nx, ny=ny, dx=dx, dy=dy,
                     dt = dt,
                     dump_freq = 200, diag_freq = dt,
                     bot_drag=bot_drag,
                     n_time_steps=400)
        else:
            drv.simulate(init_h_file=[400. for i in range(layers)], 
            init_u_file=[init_U for i in range(layers)],
            depth_file=[layers*400],
            exe=test_executable,
            wet_mask_file=[dbl_periodic_wetmask],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     dt = dt,
                     dump_freq = 200, diag_freq = dt,
                     bot_drag=bot_drag,
                     n_time_steps=400)

        hfiles = sorted(glob.glob("output/snap.h.*"))
        ufiles = sorted(glob.glob("output/snap.u.*"))
        vfiles = sorted(glob.glob("output/snap.v.*"))


        model_iteration = np.zeros(len(hfiles))

        momentum = np.zeros(len(hfiles))
        momentum_expected = np.zeros(len(hfiles))

        for counter,ufile in enumerate(ufiles):

            h = aro.interpret_raw_file(hfiles[counter], nx, ny, layers)
            u = aro.interpret_raw_file(ufile, nx, ny, layers)
            v = aro.interpret_raw_file(vfiles[counter], nx, ny, layers)

            model_iteration[counter] = float(ufile[-10:])

            momentum[counter] = (nx * ny * dx * dy * rho0 * 
                                    (np.mean(h[-1,:,:])*(np.mean(u[-1,:,:]) +
                                                        np.mean(v[-1,:,:]))))


        opt.assert_volume_conservation(nx, ny, layers, 1e-9)

        init_h = 400
        X, Y = np.meshgrid(grid.x, grid.yp1)
        init_u = init_U(X, Y)


        momentum_expected[:] = (nx * ny * dx * dy * rho0 * 
                                (np.mean(init_h)*np.mean(init_u[:,:]))*
                                    np.exp(-model_iteration*dt*bot_drag))

        test_passes = True

        try:
            np.testing.assert_allclose(momentum, momentum_expected, rtol=2e-3, atol=0)
            return
        except AssertionError as error:
            test_passes = False

            # plot output for visual inspection
            plt.figure()
            plt.plot(model_iteration*dt/(86400), momentum, '-', alpha=1,
                    label='Simulated momentum')
            plt.plot(model_iteration*dt/(86400), momentum_expected, '-', alpha=1,
                    label='Expected momentum')
            plt.legend()
            plt.xlabel('Time (days)')
            plt.ylabel('Momentum')
            plt.savefig('f_plane_momentum_test.png', dpi=150)


            plt.figure()
            plt.plot(model_iteration,momentum/momentum_expected)
            plt.xlabel('timestep')
            plt.ylabel('simulated/expected')
            plt.savefig('momentum_ratio.png')
            plt.close()

            plt.figure()
            plt.plot(model_iteration*dt/(86400),
                100.*(momentum - momentum_expected)/momentum_expected)
            plt.xlabel('Time (days)')
            plt.ylabel('percent error')
            plt.ylim(-20,80)
            plt.savefig('momentum_percent_error.png')
            plt.close()

        assert test_passes

def test_f_plane_Hypre_bot_drag_2Layers():
    test_f_plane_Hypre_bot_drag(bot_drag=1e-5, layers=2)

