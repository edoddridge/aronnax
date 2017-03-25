import os.path as p

import glob

import numpy as np
import matplotlib.pyplot as plt

import MIMutils as mim

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
<<<<<<< HEAD
import output_preservation_test as opt



def f_plane_red_grav_init_u_test():
    nx = 100
    ny = 100
    layers = 1

    dx = 10e3
    dy = 10e3

    grid = mim.Grid(nx,ny,dx,dy)

    rho0 = 1035.

    dt = 60.

    with opt.working_directory(p.join(self_path, "physics_tests/f_plane_red_grav_init_u")):
        mim_exec = "MIM_test"
        opt.run_experiment(
              write_f_plane_red_grav_init_u_input, nx, ny, layers, mim_exec)

        hfiles = sorted(glob.glob("output/snap.h.*"))
        ufiles = sorted(glob.glob("output/snap.u.*"))
        vfiles = sorted(glob.glob("output/snap.v.*"))


        model_iteration = np.zeros(len(hfiles))

        energy = np.zeros(len(hfiles))
        energy_expected = np.zeros(len(hfiles))

        momentum = np.zeros(len(hfiles))

        volume = np.zeros(len(hfiles))

        for counter,ufile in enumerate(ufiles):

            h = opt.interpret_mim_raw_file(hfiles[counter], nx, ny, layers)
            u = opt.interpret_mim_raw_file(ufile, nx, ny, layers)
            v = opt.interpret_mim_raw_file(vfiles[counter], nx, ny, layers)

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

            energy[counter] = (dx * dy * rho0 * (np.sum(np.absolute(h * ((u[1:,...]**2 + u[:-1,...]**2))/4.)/2.) + np.sum(np.absolute(h * ((v[:,1:,:]**2 + v[:,:-1,:]**2))/4.)/2.)) +  
              dx * dy * rho0 * 0.01 * np.sum(np.absolute(h - 400.)))

            momentum[counter] = dx * dy * rho0 * (np.sum(np.absolute(h * (u[1:,...] + u[:-1,...])/2.)) + np.sum(np.absolute(h * (v[:,1:,:] + v[:,:-1,:])/2.)))

            volume[counter] = np.sum(h)
        
        opt.assert_volume_conservation(nx, ny, layers, 1e-9)

        init_u = np.zeros((nx+1,ny),dtype=np.float64)
        init_u[50,50] = 0.2

        energy_expected[:] = (dx * dy * rho0 * np.sum(np.absolute(400. * ((init_u[1:,...] + init_u[:-1,...])**2)/2.)/2.))


        #assert np.amax(array_relative_error(ans, good_ans)) < rtol
        plt.figure()
        #plt.plot(model_iteration, energy_expected, '-o', alpha=0.5,
        #        label='Expected energy')
        plt.plot(model_iteration, energy, '-*', alpha=0.5,
                label='simulated energy')
        plt.legend()
        plt.xlabel('time step')
        plt.ylabel('energy')
        plt.savefig('f_plane_energy_test.png', dpi=150)

        plt.figure()
        plt.plot(model_iteration,energy/energy_expected)
        plt.xlabel('timestep')
        plt.ylabel('simulated/expected')
        plt.savefig('ratio.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration,volume)
        plt.ylabel('Volume')
        plt.xlabel('timestep')
        plt.savefig('volume.png')
        plt.close()

        plt.figure()
        plt.plot(model_iteration, momentum, '-*', alpha=0.5,
                label='simulated momentum')
        plt.legend()
        plt.xlabel('time step')
        plt.ylabel('momentum')
        plt.savefig('f_plane_momentum_test.png', dpi=150)


def write_f_plane_red_grav_wind_input(nx,ny,layers):
    opt.write_f_plane(nx, ny, 10e-4)
    opt.write_rectangular_pool(nx, ny)
    with opt.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)
    with opt.fortran_file('wind_x.bin','w') as f:
        wind_x = np.zeros((ny,nx+1),dtype=np.float64)
        wind_x[50,50] = 0.2
        f.write_record(wind_x)

        plt.figure()
        plt.pcolormesh(wind_x)
        plt.colorbar()
        plt.savefig('wind_x.png')

def write_f_plane_red_grav_init_u_input(nx,ny,layers):
    opt.write_f_plane(nx, ny, 10e-4)
    opt.write_rectangular_pool(nx, ny)
    with opt.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)
    with opt.fortran_file('init_u.bin','w') as f:
        init_u = np.zeros((ny,nx+1),dtype=np.float64)
        init_u[50,50] = 0.2
        f.write_record(init_u)

        plt.figure()
        plt.pcolormesh(init_u)
        plt.colorbar()
        plt.savefig('init_u.png')
=======
import output_preservation_test as opt
>>>>>>> parent of 03580ac... functions to write inputs
