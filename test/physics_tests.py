import os.path as p

import glob

import numpy as np
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt



def f_plane_init_u_test(physics, aro_exec, dt):
    nx = 100
    ny = 100
    layers = 1

    dx = 10e3
    dy = 10e3

    rho0 = 1035.

    grid = aro.Grid(nx, ny, layers, dx, dy)



    def init_U(X, Y, *arg):
        init_u = np.zeros(Y.shape,dtype=np.float64)
        init_u[int(grid.nx/2),int(grid.ny/2)] = 3e-5

        if not arg:
            plt.figure()
            plt.pcolormesh(init_u)
            plt.colorbar()
            plt.savefig('init_u.png')
        return init_u

    def init_V(X, Y, *arg):
        init_v = np.zeros(X.shape,dtype=np.float64)
        init_v[int(nx/2),int(ny/2)] = 3e-5

        if not arg:
            plt.figure()
            plt.pcolormesh(init_v)
            plt.colorbar()
            plt.savefig('init_v.png')
        return init_v


    with working_directory(p.join(self_path, "physics_tests/f_plane_{0}_init_u".format(physics))):
        drv.simulate(initHfile=[400.], 
            initUfile=[init_U], initVfile=[init_V], valgrind=False,
                     nx=nx, ny=ny, exe="aronnax_test", dx=dx, dy=dy)

        hfiles = sorted(glob.glob("output/snap.h.*"))
        ufiles = sorted(glob.glob("output/snap.u.*"))
        vfiles = sorted(glob.glob("output/snap.v.*"))


        model_iteration = np.zeros(len(hfiles))

        energy = np.zeros(len(hfiles))
        energy_expected = np.zeros(len(hfiles))

        momentum = np.zeros(len(hfiles))
        momentum_expected = np.zeros(len(hfiles))

        volume = np.zeros(len(hfiles))

        for counter,ufile in enumerate(ufiles):

            h = aro.interpret_raw_file(hfiles[counter], nx, ny, layers)
            u = aro.interpret_raw_file(ufile, nx, ny, layers)
            v = aro.interpret_raw_file(vfiles[counter], nx, ny, layers)

            model_iteration[counter] = float(ufile[-10:])

            # plt.figure()
            # plt.pcolormesh(grid.xp1,grid.y,u[:,:,0].transpose())
            # plt.colorbar()
            # plt.savefig('u.{0}.png'.format(ufile[-10:]),dpi=150)
            # plt.close()

            # plt.figure()
            # plt.pcolormesh(grid.x,grid.y,h[:,:,0].transpose())
            # plt.colorbar()
            # plt.savefig('h.{0}.png'.format(ufile[-10:]),dpi=150)
            # plt.close()

            energy[counter] = (dx * dy * rho0 * (np.sum(np.absolute(h * (u[1:,...]**2 + u[:-1,...]**2)/4.)/2.) + np.sum(np.absolute(h * (v[:,1:,:]**2 + v[:,:-1,:]**2)/4.)/2.)) +  
              dx * dy * rho0 * 0.01 * np.sum(np.absolute(h - 400.)))

            momentum[counter] = dx * dy * rho0 * (np.sum(np.absolute(h * (u[1:,...] + u[:-1,...])/2.)) + np.sum(np.absolute(h * (v[:,1:,:] + v[:,:-1,:])/2.)))

            volume[counter] = np.sum(h)

            # plt.figure()
            # plt.pcolormesh(grid.xp1, grid.y, np.transpose(u[:,:,0]))
            # plt.colorbar()
            # plt.savefig('output/u.{0}.png'.format(model_iteration[counter]),dpi=100)
            # plt.close()

        opt.assert_volume_conservation(nx, ny, layers, 1e-9)

        X, Y = np.meshgrid(grid.xp1, grid.y)
        init_u = init_U(X, Y, True)

        X, Y = np.meshgrid(grid.x, grid.yp1)
        init_v = init_V(X, Y, True)

        energy_expected[:] = 2. * (dx * dy * rho0 * np.sum(np.absolute(400. * ((init_u[1:,...] + init_u[:-1,...])**2)/4.)/2.))

        momentum_expected[:] = 2. * dx * dy * rho0 * (np.sum(np.absolute(400. * (init_u[1:,...] + init_u[:-1,...])/2.)))

        # print momentum[0]/momentum_expected[0]

        #assert np.amax(array_relative_error(ans, good_ans)) < rtol
        plt.figure()
        #plt.plot(model_iteration, energy_expected, '-o', alpha=0.5,
        #        label='Expected energy')
        plt.plot(model_iteration, energy, '-', alpha=1,
                label='simulated energy')
        plt.legend()
        plt.xlabel('time step')
        plt.ylabel('energy')
        plt.savefig('f_plane_energy_test.png', dpi=150)

        plt.figure()
        plt.plot(model_iteration,energy/energy_expected)
        plt.xlabel('timestep')
        plt.ylabel('simulated/expected')
        plt.savefig('energy_ratio.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,volume)
        plt.ylabel('Volume')
        plt.xlabel('timestep')
        plt.savefig('volume.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration, momentum, '-', alpha=1,
                label='simulated momentum')
        plt.legend()
        plt.xlabel('time step')
        plt.ylabel('momentum')
        plt.savefig('f_plane_momentum_test.png', dpi=150)

        plt.figure()
        plt.plot(model_iteration,momentum/momentum_expected)
        plt.xlabel('timestep')
        plt.ylabel('simulated/expected')
        plt.savefig('momentum_ratio.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,
            100.*(momentum - momentum_expected)/momentum_expected)
        plt.xlabel('timestep')
        plt.ylabel('percent error')
        plt.ylim(-20,80)
        plt.savefig('momentum_percent_error.png')
        plt.close()


def f_plane_wind_test(physics, aro_exec, nx, ny, dt):

    layers = 1

    dx = 10e3
    dy = 10e3

    grid = aro.Grid(nx, ny, layers, dx, dy)

    rho0 = 1035.

    def wind_x(X, Y, *arg):
        wind_x = np.zeros(Y.shape,dtype=np.float64)
        wind_x[int(grid.nx/2),int(grid.ny/2)] = 1e-5

        if not arg:
            plt.figure()
            plt.pcolormesh(X/1e3, Y/1e3, wind_x)
            plt.colorbar()
            plt.savefig('wind_x.png')
        return wind_x

    def wind_y(X, Y, *arg):
        wind_y = np.zeros(X.shape,dtype=np.float64)
        wind_y[int(grid.nx/2),int(grid.ny/2)] = 1e-5

        if not arg:
            plt.figure()
            plt.pcolormesh(X/1e3, Y/1e3, wind_y)
            plt.colorbar()
            plt.savefig('wind_y.png')
        return wind_y


    with opt.working_directory(p.join(self_path, "physics_tests/f_plane_{0}_wind".format(physics))):
        drv.simulate(initHfile=[400.],
            zonalWindFile=wind_x, meridionalWindFile=wind_y, valgrind=False,
                     nx=nx, ny=ny, exe="aronnax_test", dx=dx, dy=dy, dt=dt)


        hfiles = sorted(glob.glob("output/snap.h.*"))
        ufiles = sorted(glob.glob("output/snap.u.*"))
        vfiles = sorted(glob.glob("output/snap.v.*"))

        # expect the momentum to grow according to u*h*rho0 = delta_t * wind
        # F = m * a
        # m * v = h * rho0 * xlen * ylen * v 
        #       = m * a * dt 
        #       = F * dt 
        #       = wind * dx * dy * dt 


        momentum = np.zeros(len(hfiles),dtype=np.float64)
        model_iteration = np.zeros(len(hfiles),dtype=np.float64)

        momentum_expected = np.zeros(len(hfiles),dtype=np.float64)

        volume = np.zeros(len(hfiles))

        for counter,ufile in enumerate(ufiles):

            h = aro.interpret_raw_file(hfiles[counter], nx, ny, layers)
            u = aro.interpret_raw_file(ufile, nx, ny, layers)
            v = aro.interpret_raw_file(vfiles[counter], nx, ny, layers)

            model_iteration[counter] = float(ufile[-10:])

            # plt.figure()
            # plt.pcolormesh(grid.xp1,grid.y,u[:,:,0].transpose())
            # plt.colorbar()
            # plt.savefig('u.{0}.png'.format(ufile[-10:]),dpi=150)
            # plt.close()

            # plt.figure()
            # plt.pcolormesh(grid.x,grid.y,h[:,:,0].transpose())
            # plt.colorbar()
            # plt.savefig('h.{0}.png'.format(ufile[-10:]),dpi=150)
            # plt.close()

            momentum[counter] = dx * dy * rho0 * (np.sum(np.absolute(h * (u[1:,...] + u[:-1,...])/2.)) + np.sum(np.absolute(h * (v[:,1:,:] + v[:,:-1,:])/2.)))

            momentum_expected[counter] =  2.* dx * dy * 1e-5 * (model_iteration[counter] + 2) * dt

            volume[counter] = np.sum(dx * dy * h)

            # plt.figure()
            # plt.pcolormesh(grid.xp1, grid.y, np.transpose(u[:,:,0]))
            # plt.colorbar()
            # plt.savefig('output/u.{0}.png'.format(model_iteration[counter]),dpi=100)
            # plt.close()


        opt.assert_volume_conservation(nx, ny, layers, 1e-9)

        plt.figure()
        plt.plot(model_iteration, momentum_expected, '-', alpha=1,
                label='Expected momentum')
        plt.plot(model_iteration, momentum, '-', alpha=1,
                label='simulated momentum')
        plt.legend()
        plt.xlabel('time step')
        plt.ylabel('momentum')
        plt.savefig('f_plane_momentum_test.png', dpi=150)

        plt.figure()
        plt.plot(model_iteration,momentum/momentum_expected)
        plt.xlabel('timestep')
        plt.ylabel('simulated/expected')
        plt.title('final ratio = {0}'.format(str(momentum[-1]/momentum_expected[-1])))
        plt.savefig('ratio.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,
            100.*(momentum - momentum_expected)/momentum_expected)
        plt.xlabel('timestep')
        plt.ylabel('percent error')
        plt.ylim(-2,2)
        plt.savefig('percent_error.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,momentum - momentum_expected)
        plt.xlabel('timestep')
        plt.ylabel('simulated - expected')
        plt.savefig('difference.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,volume)
        plt.ylabel('Volume')
        plt.xlabel('timestep')
        plt.savefig('volume.png')
        plt.close()

# create inputs
def write_f_plane_wind_input(nx,ny,layers):
    aro.write_f_plane(nx, ny, 10e-4)
    aro.write_rectangular_pool(nx, ny)
    with aro.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)

    with aro.fortran_file('wind_x.bin','w') as f:
        wind_x = np.zeros((ny,nx+1),dtype=np.float64)
        wind_x[int(nx/2),int(ny/2)] = 1e-5
        f.write_record(wind_x)

        plt.figure()
        plt.pcolormesh(wind_x)
        plt.colorbar()
        plt.savefig('wind_x.png')

    with aro.fortran_file('wind_y.bin','w') as f:
        wind_y = np.zeros((ny+1,nx),dtype=np.float64)
        wind_y[int(nx/2),int(ny/2)] = 1e-5
        f.write_record(wind_y)

        plt.figure()
        plt.pcolormesh(wind_y)
        plt.colorbar()
        plt.savefig('wind_y.png')

def write_f_plane_init_u_input(nx,ny,layers):
    aro.write_f_plane(nx, ny, 10e-4)
    aro.write_rectangular_pool(nx, ny)
    with aro.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)
    with aro.fortran_file('init_u.bin','w') as f:
        init_u = np.zeros((ny,nx+1),dtype=np.float64)
        init_u[int(nx/2),int(ny/2)] = 3e-5
        f.write_record(init_u)

        plt.figure()
        plt.pcolormesh(init_u)
        plt.colorbar()
        plt.savefig('init_u.png')
    with aro.fortran_file('init_v.bin','w') as f:
        init_v = np.zeros((ny,nx+1),dtype=np.float64)
        init_v[int(nx/2),int(ny/2)] = 3e-5
        f.write_record(init_v)

        plt.figure()
        plt.pcolormesh(init_v)
        plt.colorbar()
        plt.savefig('init_v.png')


if __name__ == '__main__':
    #f_plane_wind_test('red_grav', aro_exec = "aronnax_core",
    #    nx = 200, ny = 200, dt = 600.)
    #f_plane_wind_test('n_layer', aro_exec = "aronnax_external_solver",
    #    nx = 50, ny = 50, dt = 100.)

    f_plane_init_u_test('red_grav', aro_exec = "aronnax_core", dt = 600.)
    #f_plane_init_u_test('n_layer', aro_exec = "aronnax_external_solver", dt = 100.)
