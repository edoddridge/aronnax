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


def run_davis_et_al_2014(nx,ny,layers,nTimeSteps,dt,simulation=None):
    with opt.working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
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
    xlen = 1530e3
    ylen = 2730e3
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

        plt.pcolormesh(X,Y,initH)
        plt.colorbar()
        plt.savefig('initH.png',dpi=150)
        plt.close()


def write_davis_wetmask(grid):
    """Write the wet mask file for a recreation of Davis et al. (2014)."""

    with opt.fortran_file('wetmask.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        # start with land everywhere and carve out space for water
        wetmask = np.zeros(X.shape, dtype=np.float64)
        # circular gyre region
        wetmask[((Y-1965e3)**2 + (X-765e3)**2) < 750e3**2] = 1
        # 150 km wide channel
        wetmask[(X-765e3)**2 < 75e3**2] = 1
        # sponge region
        wetmask[Y<780e3] = 1
        # clean up the edges
        wetmask[ 0, :] = 0
        wetmask[-1, :] = 0
        wetmask[ :, 0] = 0
        wetmask[ :,-1] = 0
        f.write_record(wetmask.astype(np.float64))

        plt.pcolormesh(X,Y,wetmask)
        plt.colorbar()
        plt.savefig('wetmask.png',dpi=150)
        plt.close()

def write_davis_wind(grid):
    """produce the wind forcing files for a recreation of Davis et al. (2014)."""
    
    L = 750e3
        
    with opt.fortran_file('tau_x.bin', 'w') as f:
        X,Y = np.meshgrid(grid.xp1,grid.y)
        r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
        theta = np.arctan2(Y-1965e3,X-765e3)

        norm = (L*(2+np.pi)/2.)/(4*np.pi);
        
        tau_x = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
        tau_x = tau_x*np.sin(theta)/norm;

        tau_x[Y<1000e3] = 0
        #-2*L*((np.pi-2)/(np.pi+2))*(-Y[Y<1200e3]/r[Y<1200e3]**2)
        #tau_x = np.sin(theta) * (r**2 * np.pi / (8. * L) + r * np.sin(r * np.pi / L)/4. + L * np.cos(r * np.pi / L)/(4. * r**2 * np.pi))
        #np.sin(theta)*(r/(4.*L) + np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
        #tau_x[Y<1200e3] = (tau_x[82,50]*(1950e3-Y[Y<1200e3]))/(r[Y<1200e3])**2
        #tau_x = tau_x/np.max(np.absolute(tau_x[:,:]))
        f.write_record(tau_x.astype(np.float64))

        plt.pcolormesh(X,Y,tau_x,cmap='RdBu_r')
        plt.colorbar()
        plt.savefig('tau_x.png',bbox_inches='tight')
        plt.close()


    with opt.fortran_file('tau_y.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.yp1)
        r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
        theta = np.arctan2(Y-1965e3,X-765e3)

        norm = (L*(2+np.pi)/2.)/(4*np.pi);
        
        tau_y = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
        tau_y = -tau_y*np.cos(theta)/norm;

        tau_y[Y<1000e3] = 0

        # tau_y = -np.cos(theta)*(np.pi*r/(8.*L) + 
        #     np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
        # tau_y[Y<1250e3] = -tau_x[82,50]*(750e3 - X[Y<1250e3])/(r[Y<1250e3])**2
        # tau_y = tau_y/np.max(np.absolute(tau_y[:,-1]))
        f.write_record(tau_y.astype(np.float64))

        plt.pcolormesh(X,Y,tau_y,cmap='RdBu_r')
        plt.colorbar()
        plt.savefig('tau_y.png',bbox_inches='tight')
        plt.close()

def write_davis_sponge(grid):
    """Produce the sponge file used by Davis et al. (2014)."""
    with opt.fortran_file('sponge_h_timescale.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        # start with land everywhere and carve out space for water
        sponge_h_timescale = np.zeros(X.shape, dtype=np.float64)
        sponge_h_timescale[Y<300e3] = 1/(3.*30.*86400.) # six month relaxation time
        f.write_record(sponge_h_timescale)

        plt.pcolormesh(X,Y,sponge_h_timescale*86400.*30.)
        plt.colorbar()
        plt.savefig('sponge_h_timescale.png',dpi=150)
        plt.close()


    with opt.fortran_file('sponge_h.bin', 'w') as f:
        X,Y = np.meshgrid(grid.x,grid.y)
        # start with land everywhere and carve out space for water
        sponge_h = 400.*np.ones(X.shape, dtype=np.float64)
        f.write_record(sponge_h)

        plt.pcolormesh(X,Y,sponge_h)
        plt.colorbar()
        plt.savefig('sponge_h.png',dpi=150)
        plt.close()

def write_wind_time_series(nTimeSteps,dt):
    with opt.fortran_file('wind_time_series.bin', 'w') as f:
        wind_time_series = 0.02375*np.ones(nTimeSteps,dtype=np.float64)
        time = np.arange(nTimeSteps)*dt
        wind_time_series[(np.mod(time,12.*30.*86400.)>8.*30.*86400.)] = 0.0125
        f.write_record(wind_time_series)
        plt.plot(time/86400./30/12, wind_time_series)
        plt.savefig('wind_time_series.png')
        plt.close()



if __name__ == '__main__':
    run_davis_et_al_2014(102, 182, 1, 1244160, 1000, 'control')
    run_davis_et_al_2014(102, 182, 1, 155520, 1000, 'control_final_five') # need to copy final outputs from the previous run
