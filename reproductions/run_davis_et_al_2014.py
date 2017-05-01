import os.path as p

import numpy as np
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Davis_et_al_2014'))

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

import output_preservation_test as opt


def davis_wetmask(X, Y):
    """The wet mask for a recreation of Davis et al. (2014)."""

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

    plt.pcolormesh(X,Y,wetmask)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('wetmask.png',dpi=150)
    plt.close()

    return wetmask

def davis_wind_x(X, Y):
    L = 750e3

    r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
    theta = np.arctan2(Y-1965e3,X-765e3)

    norm = (L*(2+np.pi)/2.)/(4*np.pi);
    
    tau_x = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
    tau_x = tau_x*np.sin(theta)/norm;

    tau_x[Y<1200e3] = 0
    #-2*L*((np.pi-2)/(np.pi+2))*(-Y[Y<1200e3]/r[Y<1200e3]**2)
    #tau_x = np.sin(theta) * (r**2 * np.pi / (8. * L) + r * np.sin(r * np.pi / L)/4. + L * np.cos(r * np.pi / L)/(4. * r**2 * np.pi))
    #np.sin(theta)*(r/(4.*L) + np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
    #tau_x[Y<1200e3] = (tau_x[82,50]*(1950e3-Y[Y<1200e3]))/(r[Y<1200e3])**2
    #tau_x = tau_x/np.max(np.absolute(tau_x[:,:]))

    plt.pcolormesh(X,Y,tau_x,cmap='RdBu_r')
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig('tau_x.png',bbox_inches='tight')
    plt.close()

    return tau_x

def davis_wind_y(X, Y):
    L = 750e3

    r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
    theta = np.arctan2(Y-1965e3,X-765e3)

    norm = (L*(2+np.pi)/2.)/(4*np.pi);
    
    tau_y = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
    tau_y = -tau_y*np.cos(theta)/norm;

    tau_y[Y<1200e3] = 0

    # tau_y = -np.cos(theta)*(np.pi*r/(8.*L) + 
    #     np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
    # tau_y[Y<1250e3] = -tau_x[82,50]*(750e3 - X[Y<1250e3])/(r[Y<1250e3])**2
    # tau_y = tau_y/np.max(np.absolute(tau_y[:,-1]))

    plt.pcolormesh(X,Y,tau_y,cmap='RdBu_r')
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.savefig('tau_y.png',bbox_inches='tight')
    plt.close()

    return tau_y

def davis_sponge_h_timescale(X, Y):
    """Produce the sponge timescale file used by Davis et al. (2014)."""
    sponge_h_timescale = np.zeros(X.shape, dtype=np.float64)
    sponge_h_timescale[Y<480e3] = 1/(1.*30.*86400.) # six month relaxation time

    plt.pcolormesh(X,Y,sponge_h_timescale*86400.*30.)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('sponge_h_timescale.png',dpi=150)
    plt.close()

    return sponge_h_timescale

def davis_sponge_h(X, Y):
    """Produce the sponge file used by Davis et al. (2014)."""
    sponge_h = 400.*np.ones(X.shape, dtype=np.float64)

    plt.pcolormesh(X,Y,sponge_h)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('sponge_h.png',dpi=150)
    plt.close()

    return sponge_h


def davis_wind_time_series(nTimeSteps,dt):
    wind_time_series = 0.02375*np.ones(nTimeSteps,dtype=np.float64)
    time = np.arange(nTimeSteps)*dt
    wind_time_series[(np.mod(time,12.*30.*86400.)>8.*30.*86400.)] = 0.0125

    plt.plot(time/86400./30/12, wind_time_series)
    plt.savefig('wind_time_series.png')
    plt.close()

    return wind_time_series

def run_davis_2014_control(nx,ny,layers,nTimeSteps,dt,simulation=None):

    #assert layers == 1
    xlen = 1530e3
    ylen = 2730e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    with opt.working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
        drv.simulate(initHfile=[400.],
                zonalWindFile=davis_wind_x, meridionalWindFile=davis_wind_y,
                wind_mag_time_series_file=[davis_wind_time_series],
                wetMaskFile=davis_wetmask,
                spongeHTimeScaleFile=[davis_sponge_h_timescale],
                spongeHFile=[davis_sponge_h],
                nx=nx, ny=ny, dx=dx, dy=dy, 
                exe='aronnax_core', 
                dt=dt, dumpFreq=int(dt*nTimeSteps/40), nTimeSteps=nTimeSteps)


def run_davis_control_final_five(nx,ny,layers,nTimeSteps,dt,simulation=None):

    #assert layers == 1
    xlen = 1530e3
    ylen = 2730e3
    dx = xlen / nx
    dy = ylen / ny

    grid = aro.Grid(nx, ny, layers, dx, dy)

    with opt.working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
        drv.simulate(initHfile='../final_state_of_control/final.h.0001244161',
                initUfile='../final_state_of_control/final.u.0001244161',
                initVfile='../final_state_of_control/final.v.0001244161',
                zonalWindFile=davis_wind_x, meridionalWindFile=davis_wind_y,
                wind_mag_time_series_file=[davis_wind_time_series],
                wetMaskFile=davis_wetmask,
                spongeHTimeScaleFile=[davis_sponge_h_timescale],
                spongeHFile=[davis_sponge_h],
                nx=nx, ny=ny, dx=dx, dy=dy, 
                exe='aronnax_core', 
                dt=dt, dumpFreq=int(86400.*5.), nTimeSteps=nTimeSteps)


if __name__ == '__main__':
    #run_davis_2014_control(102, 182, 1, 1244160, 1000, 'control')
    run_davis_control_final_five(102, 182, 1, 155520, 1000, 'control_final_five') # need to copy final outputs from the previous run
