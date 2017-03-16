import os.path as p

import numpy as np
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Davis_et_al_2014'))

import output_preservation_test as opt

import subprocess as sub

import MIMutils as mim

import time

mim_exec = p.join(root_path, "MIM")


def run_davis_et_al_2014(nx,ny,layers,nTimeSteps,dt):
    with opt.working_directory(p.join(self_path, "Davis_et_al_2014")):
        run_experiment(write_input_davis_et_al_2014, nx, ny, layers, nTimeSteps, dt)


# Generic function to run an experiment
def run_experiment(write_input, nx, ny, layers, nTimeSteps,dt):
    """ Run MIM"""
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])
    sub.check_call(["rm", "-rf", "run_finished.txt"])

    with opt.working_directory(root_path):
        sub.check_call(["make", "MIM"])
    with opt.working_directory("input"):
        write_input(nx, ny, layers,nTimeSteps,dt)
    tweak_parameters(nx, ny, layers,nTimeSteps, dt)

    then = time.time()
    sub.check_call([mim_exec])
    print "MIM execution took", time.time() - then


def tweak_parameters(nx, ny, layers, nTimeSteps,dt):
    """Alter the parameters.in file"""
    sub.check_call(
        "cat parameters.in " +
        "| sed 's/^ nx =.*,$/ nx = %d,/'" % (nx,) +
        "| sed 's/^ ny =.*,$/ ny = %d,/'" % (ny,) +
        "| sed 's/^ layers =.*,$/ layers = %d,/'" % (layers,) +
        "| sed 's/^ nTimeSteps =.*,$/ nTimeSteps = %d,/'" % (nTimeSteps,) +
        "| sed 's/^ dt =.*,$/ dt = %d,/'" % (dt,) +
        "> parameters.new", shell=True)
    sub.check_call(["mv", "parameters.new", "parameters.in"])


# Creat the inputs
def write_input_davis_et_al_2014(nx, ny, layers, nTimeSteps,dt):
    assert layers == 1
    xlen = 1.5e6
    ylen = 2.7e6
    grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

    opt.write_f_plane(nx,ny, 14.5842318e-5) # Coriolis at North Pole
    
    write_davis_wetmask(grid)
    write_davis_wind(grid)
    write_davis_sponge(grid)
    write_wind_time_series(nTimeSteps,dt)


    with opt.fortran_file('initH.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        initH = 400.*np.ones(X.shape)
        f.write_record(initH.astype(np.float64))


def write_davis_wetmask(grid):
    """Write the wet mask file for a recreation of Davis et al. (2014)."""

    with opt.fortran_file('wetmask.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        # start with land everywhere and carve out space for water
        wetmask = np.zeros(X.shape, dtype=np.float64)
        # circular gyre region
        wetmask[((Y-1950e3)**2 + (X-750e3)**2) < 750e3**2] = 1
        # 150 km wide channel
        wetmask[(X-750e3)**2 < 75e3**2] = 1
        # sponge region
        wetmask[Y<780e3] = 1
        # clean up the edges
        wetmask[ 0, :] = 0
        wetmask[-1, :] = 0
        wetmask[ :, 0] = 0
        wetmask[ :,-1] = 0
        f.write_record(wetmask)

def write_davis_wind(grid):
    """produce the wind forcing files for a recreation of Davis et al. (2014)."""
    
    L = 750e3
        
    with opt.fortran_file('tau_x.bin', 'w') as f:
        X,Y = np.meshgrid(grid.xp1,grid.y)
        r = np.sqrt((Y-1950e3)**2 + (X-750e3)**2)
        theta = np.arctan2(Y-1950e3,X-750e3)

        tau_x = np.sin(theta)*(r/4. + 
            np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
        tau_x = tau_x/np.max(np.absolute(tau_x[-1,:]))
        tau_x[Y<1250e3] = tau_x[82,50]*Y[Y<1250e3]/(r[Y<1250e3]-750e3)**2

        f.write_record(tau_x.astype(np.float64))

    with opt.fortran_file('tau_y.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.yp1)
        r = np.sqrt((Y-1950e3)**2 + (X-750e3)**2)
        theta = np.arctan2(Y-1950e3,X-750e3)

        tau_y = -np.cos(theta)*(np.pi*r/(8.*L) + 
            np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
        tau_y = tau_y/np.max(np.absolute(tau_y[:,-1]))
        tau_y[Y<1250e3] = -tau_x[82,50]*X[Y<1250e3]/(r[Y<1250e3] - 750e3)**2

        f.write_record(tau_y.astype(np.float64))
