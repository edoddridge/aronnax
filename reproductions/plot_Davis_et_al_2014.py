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

import output_preservation_test as opt

import MIMutils as mim


xlen = 1530e3
ylen = 2730e3
nx = 102
ny = 182
grid = mim.Grid(nx, ny, xlen / nx, ylen / ny)

def plt_h_cross_section():
    h_files = sorted(glob.glob("output/snap.h.*"))
    h_final = opt.interpret_mim_raw_file(h_files[-1],102,182,1,)

    plt.figure()
    for i in xrange(100,160):
        plt.plot(h_final[:,i,0])
    plt.savefig('figures/h_final.png')

def plt_state():
    h_files = sorted(glob.glob("output/snap.h.*"))
    u_files = sorted(glob.glob("output/snap.u.*"))

    h_max = np.zeros(len(h_files))

    for i in xrange(len(u_files)):
        h = opt.interpret_mim_raw_file(h_files[i],102,182,1,)
        h_max[i] = np.max(h)

        u = opt.interpret_mim_raw_file(u_files[i],102,182,1,)

        X,Y = np.meshgrid(grid.x,grid.y)

        CS = plt.contour(X,Y,np.transpose(h[:,:,0]),colors='k')
        plt.clabel(CS, inline=1, fontsize=10)
        X,Y = np.meshgrid(grid.xp1,grid.y)

        plt.pcolormesh(X,Y,np.transpose(u[:,:,0]),cmap='RdBu_r',
            vmin = -0.012, vmax = 0.012)
        plt.colorbar()

        plt.savefig('figures/state_{0}.png'.format(u_files[i][-10:]),dpi=150)
        plt.close()

    plt.plot(h_max)
    plt.savefig('figures/h_max.png',dpi=150)
    plt.close()



if __name__ == '__main__':
    os.chdir('Davis_et_al_2014')
    sub.check_call(["rm", "-rf", "figures/"])
    sub.check_call(["mkdir", "-p", "figures/"])

    plt_h_cross_section()
    plt_state()