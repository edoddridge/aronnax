import os.path as p

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Davis_et_al_2014'))

import aronnax.driver as drv
from aronnax.utils import working_directory


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

    plt.figure()
    plt.pcolormesh(X/1e3, Y/1e3, wetmask, cmap='Greys_r')
    #plt.colorbar()
    plt.xlim(0,1500)
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.savefig('wetmask.png', dpi=150, bbox_inches='tight')
    plt.close()

    return wetmask

def davis_wind_x(X, Y):
    L = 1500e3

    r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
    theta = np.arctan2(Y-1965e3, X-765e3)

    # From Pete's code
    # This wind forcing gave an upper layer that was slightly too thin.
    # base_wind = ((2.*L)/(np.pi*r))*(
    #                     np.pi*r*np.sin(np.pi*r/L)/(8.*L) - 
    #                     (np.sin(np.pi*r/(2.*L))**2)/4. + 
    #                     np.pi**2 * r**2 / (16.*L))

    # base_wind[Y-1965e3<-L] = np.max(base_wind[r<L])
    # base_wind = base_wind/np.mean(np.fabs(base_wind[r<L]))
    # tau_x = base_wind*np.sin(theta)


    # my attempt
    # norm = (L*(2+np.pi)/2.)/(4*np.pi);
    
    tau_x = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
    tau_x = tau_x/np.mean(np.fabs(tau_x[r<L]))
    tau_x = tau_x*np.sin(theta)

    tau_x[Y<1200e3] = 0

    #-2*L*((np.pi-2)/(np.pi+2))*(-Y[Y<1200e3]/r[Y<1200e3]**2)
    #tau_x = np.sin(theta) * (r**2 * np.pi / (8. * L) + r * np.sin(r * np.pi / L)/4. + L * np.cos(r * np.pi / L)/(4. * r**2 * np.pi))
    #np.sin(theta)*(r/(4.*L) + np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
    #tau_x[Y<1200e3] = (tau_x[82,50]*(1950e3-Y[Y<1200e3]))/(r[Y<1200e3])**2
    #tau_x = tau_x/np.max(np.absolute(tau_x[:,:]))

    plt.figure()
    plt.pcolormesh(X/1e3,Y/1e3,tau_x,cmap='RdBu_r')
    CB = plt.colorbar()
    CB.set_label('Normalised pattern for x component of wind stress')
    plt.xlim(0,1500)
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.savefig('tau_x.png',bbox_inches='tight')
    plt.close()

    return tau_x

def davis_wind_y(X, Y):
    L = 1500e3

    r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
    theta = np.arctan2(Y-1965e3, X-765e3)

    # From Pete's code
    # This wind forcing gave an upper layer that was slightly too thin.
    # base_wind = ((2.*L)/(np.pi*r))*(
    #                     np.pi*r*np.sin(np.pi*r/L)/(8.*L) - 
    #                     (np.sin(np.pi*r/(2.*L))**2)/4. + 
    #                     np.pi**2 * r**2 / (16.*L))


    # base_wind[Y-1965e3<-L] = np.max(base_wind[r<L])
    # base_wind = base_wind/np.mean(np.fabs(base_wind[r<L]))

    # tau_y = -base_wind*np.cos(theta)

    
    tau_y = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
    tau_y = tau_y/np.mean(np.fabs(tau_y[r<L]))
    tau_y = -tau_y*np.cos(theta)

    tau_y[Y<1200e3] = 0



    plt.figure()
    plt.pcolormesh(X/1e3,Y/1e3,tau_y,cmap='RdBu_r')
    CB = plt.colorbar()
    CB.set_label('Normalised pattern for y component of wind stress')
    plt.xlim(0,1500)
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.savefig('tau_y.png',bbox_inches='tight')
    plt.close()

    return tau_y

def davis_sponge_h_timescale(X, Y):
    """Produce the sponge timescale file used by Davis et al. (2014)."""
    sponge_h_timescale = np.zeros(X.shape, dtype=np.float64)
    sponge_h_timescale[Y<480e3] = 1/(1.*30.*86400.) # six month relaxation time

    plt.figure()
    plt.pcolormesh(X,Y,sponge_h_timescale*86400.*30.)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('sponge_h_timescale.png',dpi=150)
    plt.close()

    return sponge_h_timescale

def davis_sponge_h(X, Y):
    """Produce the sponge file used by Davis et al. (2014)."""
    sponge_h = 400.*np.ones(X.shape, dtype=np.float64)

    plt.figure()
    plt.pcolormesh(X,Y,sponge_h)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('sponge_h.png', dpi=150)
    plt.close()

    return sponge_h


def davis_wind_time_series(n_time_steps, dt):
    wind_time_series = 0.02375*np.ones(n_time_steps, dtype=np.float64)
    time = np.arange(n_time_steps)*dt
    wind_time_series[(np.mod(time,12.*30.*86400.)>8.*30.*86400.)] = 0.0125

    plt.figure()
    plt.plot(time/86400./30/12, wind_time_series)
    plt.title('Wind stress time series')
    plt.xlabel('Time (years)')
    plt.ylabel('Wind stress (N / m^2)')
    plt.savefig('wind_time_series.png')
    plt.close()

    return wind_time_series

def run_davis_2014_control(nx, ny, layers, n_time_steps, dt, simulation=None):

    #assert layers == 1
    xlen = 1530e3
    ylen = 2730e3
    dx = xlen / nx
    dy = ylen / ny

    # grid = aro.Grid(nx, ny, layers, dx, dy)

    with working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
        drv.simulate(init_h_file=[400.],
                zonal_wind_file=[davis_wind_x], meridional_wind_file=[davis_wind_y],
                wind_mag_time_series_file=[davis_wind_time_series],
                wet_mask_file=[davis_wetmask],
                sponge_h_time_scale_file=[davis_sponge_h_timescale],
                sponge_h_file=[davis_sponge_h],
                nx=nx, ny=ny, dx=dx, dy=dy, 
                exe='aronnax_core', 
                dt=dt, dump_freq=int(dt*n_time_steps/40), n_time_steps=n_time_steps)


def run_davis_control_final_five(nx, ny, layers, n_time_steps, dt, simulation=None):

    #assert layers == 1
    xlen = 1530e3
    ylen = 2730e3
    dx = xlen / nx
    dy = ylen / ny

    # grid = aro.Grid(nx, ny, layers, dx, dy)

    with working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
        drv.simulate(init_h_file='../final_state_of_control/final.h.0001244161',
                init_u_file='../final_state_of_control/final.u.0001244161',
                init_v_file='../final_state_of_control/final.v.0001244161',
                zonal_wind_file=[davis_wind_x], meridional_wind_file=[davis_wind_y],
                wind_mag_time_series_file=[davis_wind_time_series],
                wet_mask_file=[davis_wetmask],
                sponge_h_time_scale_file=[davis_sponge_h_timescale],
                sponge_h_file=[davis_sponge_h],
                nx=nx, ny=ny, dx=dx, dy=dy, 
                exe='aronnax_core', 
                dt=dt, dump_freq=int(86400.*5.), n_time_steps=n_time_steps)


if __name__ == '__main__':
    #run_davis_2014_control(102, 182, 1, 1244160, 1000, 'control')
    run_davis_control_final_five(102, 182, 1, 155520, 1000, 'control_final_five') # need to copy final outputs from the previous run
