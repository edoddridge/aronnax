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


# The examples
#
# N.B. These configurations are equivalent to several of those in the test
# suite, but run at a higher resolution.

executable = "aronnax_external_solver"


def bump(X, Y):
    return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

def beta_plane_bump_red_grav():
    print "Running reduced gravity beta-plane bump example"
    xlen = 1e6
    ylen = 1e6
    nx = 100
    ny = 100
    layers = 1 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_bump")):
        drv.simulate(initHfile=[bump],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        with working_directory("figures"):
            plt_output(grid)

def beta_plane_bump_n_layer():
    print "Running n-layer beta-plane bump example"
    xlen = 1e6
    ylen = 1e6
    nx = 100
    ny = 100
    layers = 2 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    with working_directory(p.join(self_path, "n_layer/beta_plane_bump")):
        drv.simulate(initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        with working_directory("figures"):
            plt_output(grid)

def beta_plane_gyre_red_grav():
    print "Running reduced gravity beta-plane gyre example"
    xlen = 1e6
    ylen = 2e6
    nx = 100 
    ny = 200
    layers = 1 # also set in aronnax.conf file
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(Y)))

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_gyre")):
        drv.simulate(zonalWindFile=[wind],
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny,
                     exe=executable)
        with working_directory("figures"):
            plt_output(grid, colour_lim=20)

def beta_plane_gyre_n_layer():
    print "Running n-layer beta-plane gyre example"
    xlen = 1e6
    ylen = 2e6
    nx = 100
    ny = 200
    layers = 2
    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

    with working_directory(p.join(self_path, "n_layer/beta_plane_gyre")):
        drv.simulate(zonalWindFile=[wind],
                     exe=executable,
                     nx=nx, ny=ny, layers=layers,
                     dx=xlen/nx, dy=ylen/ny,
                     # nTimeSteps = 20001, dumpFreq = 6e5, avFreq = 48e5
                     # uncomment previous line to reproduce simulation shown in manual
                     )
        with working_directory("figures"):
            plt_output(grid, colour_lim=20)

def plt_output(grid, colour_lim=2):
    h_files = sorted(glob.glob("../output/snap.h.*"))
    v_files = sorted(glob.glob("../output/snap.v.*"))


    # plot each state of the run

    for i in xrange(len(v_files)):
        h = aro.interpret_raw_file(h_files[i], grid.nx, grid.ny, grid.layers)
        v = aro.interpret_raw_file(v_files[i], grid.nx, grid.ny, grid.layers)

        X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)

        plt.figure()
        CS = plt.contour(X,Y,np.transpose(h[:,:,0]),colors='k')
        plt.clabel(CS, inline=1, inline_spacing=17,
         fontsize=10, fmt=r'%3.0f')
        X,Y = np.meshgrid(grid.x/1e3, grid.yp1/1e3)

        plt.pcolormesh(X,Y,np.transpose(v[:,:,0])*100., cmap='RdBu_r'
            ,vmin = -colour_lim, vmax = colour_lim)
        CB = plt.colorbar()
        CB.set_label('y component of velocity (cm / s)')
        plt.axes().set_aspect('equal')
        plt.xlabel('x coordinate (km)')
        plt.ylabel('y coordinate (km)')
        plt.title('timestep {0}'.format(v_files[i][-10:]))

        plt.savefig('state_{0}.png'.format(v_files[i][-10:]), dpi=150,
            bbox_inches='tight')
        plt.close()

    try:
        sub.check_call(["convert", "-delay", "30", "-loop", "0", "*.png", "animation.gif"])
    except:
        print "failed to make animation"


if __name__ == '__main__':
    beta_plane_bump_red_grav()
    beta_plane_bump_n_layer()
    beta_plane_gyre_red_grav()
    beta_plane_gyre_n_layer()
