import os.path as p

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import aronnax.driver as drv
from aronnax.utils import working_directory


def Yang_et_al_spin_up():

    nx = 200
    ny = 200
    layers = 2
    xlen = 1000e3
    ylen = 1000e3
    dx = xlen / nx
    dy = ylen / ny

    def wind_x(X, Y):

        L = 500e3

        r = np.sqrt((Y-500e3)**2 + (X-500e3)**2)
        theta = np.arctan2(Y-500e3, X-500e3)

        tau_x = np.sin(0.5*np.pi*r/L)
        tau_x = tau_x*np.sin(theta)*0.025
        tau_x[r>L] = 0

        plt.pcolormesh(X, Y, tau_x)
        plt.colorbar()
        plt.savefig('tau_x.pdf')
        plt.close()

        return tau_x

    def wind_y(X, Y):
        L = 500e3

        r = np.sqrt((Y-500e3)**2 + (X-500e3)**2)
        theta = np.arctan2(Y-500e3, X-500e3)

        tau_y = np.sin(0.5* np.pi*r/L)
        tau_y = -tau_y*np.cos(theta)*0.025
        tau_y[r>L] = 0

        plt.pcolormesh(X, Y, tau_y)
        plt.colorbar()
        plt.savefig('tau_y.pdf')
        plt.close()

        return tau_y

    def wetmask(X, Y):
        """The wet mask."""

        # start with land everywhere and carve out space for water
        wetmask = np.zeros(X.shape, dtype=np.float64)
        # circular gyre region
        wetmask[((Y-500e3)**2 + (X-500e3)**2) < 500e3**2] = 1

        # clean up the edges
        wetmask[ 0, :] = 0
        wetmask[-1, :] = 0
        wetmask[ :, 0] = 0
        wetmask[ :,-1] = 0

        return wetmask

    def bathymetry(X,Y):
        mask = wetmask(X, Y)
        r = np.sqrt((Y-500e3)**2 + (X-500e3)**2)

        # creating slope near southern boundary
        depth = np.minimum(4000*np.ones(X.shape), 15250 - 0.03*r)
        # set minimum depth to 250 m (also the model won't run if depthFile
        #   contains negative numbers)
        depth = np.maximum(depth, 250)

        plt.pcolormesh(X, Y, np.ma.masked_where(mask==0, depth))
        plt.colorbar()
        plt.savefig('depth.pdf')
        plt.close()

        return depth

    def coriolis(X,Y):
        r = np.sqrt((Y-500e3)**2 + (X-500e3)**2)

        f0 = 14.5842318e-5 # at north pole
        beta = f0*np.cos(np.pi*85/180)/6371e3

        f = f0 - beta*r

        plt.pcolormesh(X, Y, f)
        plt.colorbar()
        plt.savefig('coriolis.pdf')
        plt.close()
        return f


    with working_directory(p.join(self_path, 'spin_up')):
        drv.simulate(
                # give it flat layers and let it squash them to fit
                initHfile=[400, 3600],
                zonalWindFile=[wind_x], meridionalWindFile=[wind_y],
                wind_depth=30,
                wetMaskFile=[wetmask],
                depthFile=[bathymetry],
                fUfile=[coriolis],
                fVfile=[coriolis],
                nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                exe='aronnax_external_solver', 
                dt=50, nTimeSteps=12441600)

if __name__ == '__main__':
    Yang_et_al_spin_up()
