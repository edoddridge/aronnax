import os
import os.path as p

import glob

import numpy as np
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Davis_et_al_2014'))

import subprocess as sub

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

xlen = 1530e3
ylen = 2730e3
nx = 102
ny = 182
grid = aro.Grid(nx, ny, 1, xlen / nx, ylen / ny)

def plt_h_cross_section(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    h_final = aro.interpret_raw_file(h_files[-1],102,182,1,)

    plt.figure()
    for i in xrange(100,160):
        plt.plot(h_final[:,i,0])
    plt.savefig('{0}figures/h_final.png'.format(simulation))

def plt_state(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    u_files = sorted(glob.glob("{0}output/snap.u.*".format(simulation)))

    h_max = np.zeros(len(h_files))

    for i in xrange(len(u_files)):
        h = aro.interpret_raw_file(h_files[i],102,182,1,)
        h_max[i] = np.max(h)

        u = aro.interpret_raw_file(u_files[i],102,182,1,)

        X,Y = np.meshgrid(grid.x,grid.y)

        CS = plt.contour(X,Y,np.transpose(h[:,:,0]),colors='k')
        plt.clabel(CS, inline=1, fontsize=10)
        X,Y = np.meshgrid(grid.xp1,grid.y)

        plt.pcolormesh(X,Y,np.transpose(u[:,:,0]),cmap='RdBu_r',
            vmin = -0.02, vmax = 0.02)
        plt.colorbar()

        plt.savefig('{0}figures/state_{1}.png'.format(simulation,u_files[i][-10:]),dpi=150)
        plt.close()

    plt.plot(h_max)
    plt.savefig('{0}figures/h_max.png'.format(simulation),dpi=150)
    plt.close()



if __name__ == '__main__':
    os.chdir('Davis_et_al_2014')
    sub.check_call(["rm", "-rf", "control/figures/"])
    sub.check_call(["mkdir", "-p", "control/figures/"])

    plt_h_cross_section('control/')
    plt_state('control/')