import os
import os.path as p

import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
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
    plt.close()

def plt_state(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    v_files = sorted(glob.glob("{0}output/snap.v.*".format(simulation)))

    h_max = np.zeros(len(h_files))

    for i in xrange(len(v_files)):
        h = aro.interpret_raw_file(h_files[i],102,182,1,)
        h_max[i] = np.max(h)

    # plot the final state of the run
    i = len(h_files) - 1
    v = aro.interpret_raw_file(v_files[i],102,182,1,)

    X,Y = np.meshgrid(grid.x/1e3,grid.y/1e3)

    CS = plt.contour(X,Y,np.transpose(h[:,:,0]),colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
    X,Y = np.meshgrid(grid.x/1e3,grid.yp1/1e3)

    plt.pcolormesh(X,Y,np.transpose(v[:,:,0])*100.,cmap='RdBu_r'
        ,vmin = -2, vmax = 2)
    CB = plt.colorbar()
    CB.set_label('y component of velocity (cm / s)')
    plt.xlim(0,1500)
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')

    plt.savefig('{0}figures/state_{1}.png'.format(simulation,v_files[i][-10:]),dpi=150,
        bbox_inches='tight')
    plt.close()

    if simulation == 'control_final_five/':
        plt.plot(np.arange(360.)*5./360.,h_max)
        plt.xlabel('Time (years)')
    else:
        plt.plot(h_max)
    plt.ylabel('Depth and centre of gyre (m)')
    plt.savefig('{0}figures/h_max.png'.format(simulation),dpi=150)
    plt.close()

def plot_channel_transport(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    v_files = sorted(glob.glob("{0}output/snap.v.*".format(simulation)))

    transport = np.zeros(len(h_files))

    for i in xrange(len(v_files)):
        h = aro.interpret_raw_file(h_files[i],102,182,1,)
        v = aro.interpret_raw_file(v_files[i], 102, 182, 1)


        transport[i] = np.sum(h[:,70,0]*(v[:,53,0]+v[:,54,0])/2.)*grid.dx/1e6


    if simulation == 'control_final_five/':
        plt.plot(np.arange(360.)*5./360.,transport)
        plt.hlines(0,0,5)
        plt.xlabel('Time (years)')
    else:
        plt.plot(transport)
    plt.ylabel('Transport trhough the channel (Sv)')
    plt.savefig('{0}figures/transport.png'.format(simulation),dpi=150)
    plt.close()


if __name__ == '__main__':
    os.chdir('Davis_et_al_2014')
    sub.check_call(["rm", "-rf", "control/figures/"])
    sub.check_call(["mkdir", "-p", "control/figures/"])
    sub.check_call(["rm", "-rf", "control_final_five/figures/"])
    sub.check_call(["mkdir", "-p", "control_final_five/figures/"])

    plt_h_cross_section('control/')
    plt_state('control/')
    plot_channel_transport('control/')

    plt_h_cross_section('control_final_five/')
    plt_state('control_final_five/')
    plot_channel_transport('control_final_five/')