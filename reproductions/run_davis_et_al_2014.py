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

    assert layers == 1
    xlen = 1530e3
    ylen = 2730e3

    grid = mim.Grid(nx, ny, layers, xlen / nx, ylen / ny)


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
        f.write_record(wetmask.astype(np.float64))

        plt.pcolormesh(X,Y,wetmask)
        plt.colorbar()
        plt.savefig('wetmask.png',dpi=150)
        plt.close()

    def davis_wind_x(X, Y):
        L = 750e3
        
        r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
        theta = np.arctan2(Y-1965e3,X-765e3)

        norm = (L*(2+np.pi)/2.)/(4*np.pi);
        
        tau_x = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
        tau_x = tau_x*np.sin(theta)/norm;

        tau_x[Y<1150e3] = 0
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

    def davis_wind_y(X, Y):
        L = 750e3

        r = np.sqrt((Y-1965e3)**2 + (X-765e3)**2)
        theta = np.arctan2(Y-1965e3,X-765e3)

        norm = (L*(2+np.pi)/2.)/(4*np.pi);
        
        tau_y = (L/(2*np.pi))*np.sin(np.pi*r/L)+(r/4)+(L/(4*np.pi))*(np.cos(np.pi*r/L)-1);
        tau_y = -tau_y*np.cos(theta)/norm;

        tau_y[Y<1150e3] = 0

        # tau_y = -np.cos(theta)*(np.pi*r/(8.*L) + 
        #     np.sin(np.pi*r/(2.*L))/4. + np.cos(np.pi*r/(2.*L))/(r*8.))
        # tau_y[Y<1250e3] = -tau_x[82,50]*(750e3 - X[Y<1250e3])/(r[Y<1250e3])**2
        # tau_y = tau_y/np.max(np.absolute(tau_y[:,-1]))
        f.write_record(tau_y.astype(np.float64))

        plt.pcolormesh(X,Y,tau_y,cmap='RdBu_r')
        plt.colorbar()
        plt.savefig('tau_y.png',bbox_inches='tight')
        plt.close()

    def davis_sponge_h_timescale(X, Y):
        """Produce the sponge timescale file used by Davis et al. (2014)."""
        sponge_h_timescale = np.zeros(X.shape, dtype=np.float64)
        sponge_h_timescale[Y<300e3] = 1/(3.*30.*86400.) # six month relaxation time
        f.write_record(sponge_h_timescale)

        plt.pcolormesh(X,Y,sponge_h_timescale*86400.*30.)
        plt.colorbar()
        plt.savefig('sponge_h_timescale.png',dpi=150)
        plt.close()

    def davis_sponge_h(X, Y):
        """Produce the sponge timescale file used by Davis et al. (2014)."""
        sponge_h = 400.*np.ones(X.shape, dtype=np.float64)
        f.write_record(sponge_h)

        plt.pcolormesh(X,Y,sponge_h)
        plt.colorbar()
        plt.savefig('sponge_h.png',dpi=150)
        plt.close()


    def wind_time_series(nTimeSteps,dt):
        wind_time_series = 0.02375*np.ones(nTimeSteps,dtype=np.float64)
        time = np.arange(nTimeSteps)*dt
        wind_time_series[(np.mod(time,12.*30.*86400.)>8.*30.*86400.)] = 0.0125
        f.write_record(wind_time_series)
        plt.plot(time/86400./30/12, wind_time_series)
        plt.savefig('wind_time_series.png')
        plt.close()

    with opt.working_directory(p.join(self_path, 
        "Davis_et_al_2014/{0}".format(simulation))):
        drv.simulate(initHfile=[400.],
                zonalWindFile=wind_x, meridionalWindFile=wind_y,
                fUfile = [14.5842318e-5], fVfile = [14.5842318e-5],
                nx=nx, ny=ny, dx=dx, dy=dy, 
                exe=aro_exec, 
                dt=dt, dumpFreq=int(dt*nTimeSteps/50), nTimeSteps=nTimeSteps)



if __name__ == '__main__':
    run_davis_et_al_2014(102, 182, 1, 1244160, 1000, 'control')
    run_davis_et_al_2014(102, 182, 1, 155520, 1000, 'control_final_five') # need to copy final outputs from the previous run
