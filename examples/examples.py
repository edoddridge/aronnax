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

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_bump")):
        drv.simulate(initHfile=[bump],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        # todo: plot output


def beta_plane_bump_n_layer():
    print "Running n-layer beta-plane bump example"
    xlen = 1e6
    ylen = 1e6
    nx = 100
    ny = 100

    with working_directory(p.join(self_path, "n_layer/beta_plane_bump")):
        drv.simulate(initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)],
                     exe=executable,
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny)
        # todo: plot output


def beta_plane_gyre_red_grav():
    print "Running reduced gravity beta-plane gyre example"
    xlen = 1e6
    ylen = 2e6
    nx = 100 
    ny = 200
    grid = aro.Grid(nx, ny, 1, xlen / nx, ylen / ny)

    def wind(_, Y):
        return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(Y)))

    with working_directory(p.join(self_path,
                            "reduced_gravity/beta_plane_gyre")):
        drv.simulate(zonalWindFile=[wind],
                     nx=nx, ny=ny,
                     dx=xlen/nx, dy=ylen/ny,
                     exe=executable)
        # todo: plot output

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
                     dx=xlen/nx, dy=ylen/ny)
        # todo: plot output


if __name__ == '__main__':
    beta_plane_bump_red_grav()
    beta_plane_bump_n_layer()
    beta_plane_gyre_red_grav()
    beta_plane_gyre_n_layer()
