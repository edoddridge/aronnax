from __future__ import print_function

import os.path as p
import subprocess as sub
from builtins import range

import glob

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))


# The examples
#
# N.B. These configurations are equivalent to several of those in the test
# suite, but run at a higher resolution.

executable = "aronnax_external_solver"


def bump(X, Y):

    iniH = 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

    plt.figure()
    CS = plt.contour(X/1e3,Y/1e3,iniH, np.arange(490., 520, 2.), colors='k')
    plt.clabel(CS, inline=1, inline_spacing=17,
     fontsize=10, fmt=r'%3.0f')

    im = plt.pcolormesh(X/1e3,Y/1e3, iniH, vmin = 500, vmax = 520)
    im.set_edgecolor('face')
    plt.colorbar()
    plt.title('Initial depth of upper layer (m)')
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')

    plt.savefig('initial_h_bump.png', bbox_inches='tight', dpi=200)

    return iniH



def beta_plane_bump_red_grav():
    print("Running reduced gravity beta-plane bump example")
    xlen = 1e6
    ylen = 1e6
    nx = 100
    ny = 100
    layers = 1 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_bump")):
        drv.simulate(init_h_file=[bump],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        with working_directory("figures"):
            plt_output(grid, 'red-grav-bump')

def beta_plane_bump_n_layer():
    print("Running n-layer beta-plane bump example")
    xlen = 1e6
    ylen = 1e6
    nx = 100
    ny = 100
    layers = 2 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    with working_directory(p.join(self_path, "n_layer/beta_plane_bump")):
        drv.simulate(init_h_file=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        with working_directory("figures"):
            plt_output(grid, 'n-layer-bump')

def beta_plane_gyre_red_grav():
    print("Running reduced gravity beta-plane gyre example")
    xlen = 1e6
    ylen = 2e6
    nx = 100 
    ny = 200
    layers = 1 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    def wind(_, Y):
        
        wind_forcing = 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

        plt.figure()
        plt.plot(Y[:,1]/1e3, wind_forcing[:,1])
        plt.xlabel('y coordinate (km)')
        plt.ylabel('Wind forcing (N/m^2)')
        plt.savefig('twin_gyre_wind_forcing.png', 
            bbox_inches='tight', dpi=200)
        plt.close()

        return wind_forcing

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_gyre")):
        drv.simulate(zonal_wind_file=[wind],
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny,
                     exe=executable)
        with working_directory("figures"):
            plt_output(grid, 'red-grav-twin-gyre', colour_lim=20)

def beta_plane_gyre_n_layer():
    print("Running n-layer beta-plane gyre example")
    xlen = 1e6
    ylen = 2e6
    nx = 100
    ny = 200
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    def wind(_, Y):
        
        wind_forcing = 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

        plt.figure()
        plt.plot(Y[:,1]/1e3, wind_forcing[:,1])
        plt.xlabel('Latitude (km)')
        plt.ylabel('Wind forcing (N/m^2)')
        plt.savefig('twin_gyre_wind_forcing.png', 
            bbox_inches='tight', dpi=200)
        plt.close()

        return wind_forcing

    with working_directory(p.join(self_path, "n_layer/beta_plane_gyre")):
        drv.simulate(zonal_wind_file=[wind],
                     exe=executable,
                     nx=nx, ny=ny, layers=layers,
                     dx=xlen/nx, dy=ylen/ny,
                     n_time_steps = 20001, dump_freq = 6e5, av_freq = 48e5
                     # NB: this takes quite a long time to run.
                     )
        with working_directory("figures"):
            plt_output(grid, 'n-layer-twin-gyre', colour_lim=20)

def plt_output(grid, sim_name, colour_lim=2):
    h_files = sorted(glob.glob("../output/snap.h.*"))
    v_files = sorted(glob.glob("../output/snap.v.*"))


    # plot each state of the run

    for i in range(len(v_files)):
        h = aro.interpret_raw_file(h_files[i], grid.nx, grid.ny, grid.layers)
        v = aro.interpret_raw_file(v_files[i], grid.nx, grid.ny, grid.layers)

        X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)

        plt.figure()
        CS = plt.contour(X,Y,h[0,:,:],colors='k')
        plt.clabel(CS, inline=1, inline_spacing=17,
         fontsize=10, fmt=r'%3.0f')
        X,Y = np.meshgrid(grid.x/1e3, grid.yp1/1e3)

        im = plt.pcolormesh(X,Y,v[0,:,:]*100., cmap='RdBu_r'
            ,vmin = -colour_lim, vmax = colour_lim)
        im.set_edgecolor('face')
        CB = plt.colorbar()
        CB.set_label('y component of velocity (cm / s)')
        plt.axes().set_aspect('equal')
        plt.xlabel('x coordinate (km)')
        plt.ylabel('y coordinate (km)')
        plt.title('timestep {0}'.format(v_files[i][-10:]))

        plt.savefig('state_{0}.png'.format(v_files[i][-10:]), dpi=150,
            bbox_inches='tight')
        if i==len(v_files)-1:
            plt.savefig('{0}.pdf'.format(sim_name), dpi=150,
                bbox_inches='tight')
        plt.close()

    try:
        sub.check_call(["convert", "-delay", "30", "-loop", "0", "state_*.png", "{0}.gif".format(sim_name)])
    except:
        print("failed to make animation")


if __name__ == '__main__':
    beta_plane_bump_red_grav()
    beta_plane_bump_n_layer()
    beta_plane_gyre_red_grav()
    beta_plane_gyre_n_layer()
