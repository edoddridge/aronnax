import os
import os.path as p

import glob
from builtins import range

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Yang_et_al_2016'))

import subprocess as sub

import aronnax as aro


## Pretty plots
plt.rcParams['figure.figsize'] = (12, 12) # set default figure size to 12x12 inches
plt.rc('text',usetex=True)
#font = {'family':'serif','size':16}
font = {'family':'serif','size':16, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':14})
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']


def wetmask(X, Y):
    """The wet mask."""

    # start with land everywhere and carve out space for water
    wetmask = np.zeros(X.shape, dtype=np.float64)
    # circular gyre region
    wetmask[((Y-500)**2 + (X-500)**2) < 500**2] = 1

    # clean up the edges
    wetmask[ 0, :] = 0
    wetmask[-1, :] = 0
    wetmask[ :, 0] = 0
    wetmask[ :,-1] = 0

    return wetmask

nx = 200
ny = 200
layers = 2
xlen = 1000e3
ylen = 1000e3
dx = xlen / nx
dy = ylen / ny

grid = aro.Grid(nx, ny, layers, dx, dy)

X_h, Y_h = np.meshgrid(grid.x/1e3, grid.y/1e3)
mask_h = wetmask(X_h, Y_h)
X_v, Y_v = np.meshgrid(grid.x/1e3, grid.yp1/1e3)
mask_v = wetmask(X_v, Y_v)

def plt_h_cross_section(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    h_final = aro.interpret_raw_file(h_files[-1], nx, ny, layers)

    plt.figure()
    for i in range(100,160):
        plt.plot(-h_final[0,i,:])
    plt.title('Timestep {}: Day {:05d}'.format(h_files[-1][-10:], len(h_files)))
    plt.savefig('{0}figures/h_final.png'.format(simulation))
    plt.close()

    plt.pcolormesh(X_h, Y_h, h_final[0,:,:])
    plt.colorbar()
    plt.contour(X_h, Y_h, h_final[0,:,:], np.arange(0,800,50), colors='k')
    plt.title('Timestep {}: Day {:05d}'.format(h_files[-1][-10:], len(h_files)))
    plt.savefig('{0}figures/h_final_pcolor.png'.format(simulation))
    plt.close()


def plt_state(simulation=None):
    h_files = sorted(glob.glob("{0}output/snap.h.*".format(simulation)))
    v_files = sorted(glob.glob("{0}output/snap.v.*".format(simulation)))

    h_max = np.zeros(len(h_files))

    for i in range(len(v_files)):
        h = aro.interpret_raw_file(h_files[i], nx, ny, layers)
        h_max[i] = h[0,100,100] #np.max(h[0,:,:])
    
    years = np.arange(len(h_max))/360
    plt.figure()
    plt.plot(years, h_max)
    plt.ylabel('Depth at centre of gyre (m)')
    plt.xlabel('Time (years)')
    plt.savefig('{0}figures/h_max.png'.format(simulation), dpi=300)
    plt.close()

    for i in range(0, len(v_files)):
        h = aro.interpret_raw_file(h_files[i], nx, ny, layers)

        v = aro.interpret_raw_file(v_files[i], nx, ny, layers)

        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
        ax1.pcolormesh(X_h, Y_h, np.ma.masked_where(mask_h==1, mask_h), cmap='YlOrBr_r',
            vmin=0, vmax=1)
        CS = ax1.contour(X_h, Y_h, h[0,:,:], np.arange(0,800,50), colors='k')

        im = ax1.pcolormesh(X_v, Y_v, np.ma.masked_where(mask_v==0, v[0,:,:]*100.),
            cmap='RdBu_r', vmin = -75, vmax = 75)
        CB = plt.colorbar(im, ax=ax1, orientation='horizontal')
        CB.set_label('y component of velocity (cm / s)')

        ax1.plot(500, 500, 'ro', ms=8)
        ax1.set_xlim(0,1000)
        ax1.set_aspect('equal')
        ax1.set_xlabel('x coordinate (km)')
        ax1.set_ylabel('y coordinate (km)')

        ax2.plot(years, h_max)
        ax2.plot(years[i], h_max[i], 'ro', ms=8)
        ax2.set_ylabel('Depth at centre of gyre (m)')
        ax2.set_xlabel('Time (years)')

        f.suptitle('Year {:02d} Day {:03d} (Timestep {})'.format(
            int(np.floor(i/360)), int(np.mod(i, 360)), v_files[i][-10:]))

        f.savefig('{}figures/state_L0_{:06d}.png'.format(simulation,i), dpi=300,
            bbox_inches='tight')
        plt.close('all')


        plt.figure()
        plt.pcolormesh(X_h, Y_h, np.ma.masked_where(mask_h==1, mask_h), cmap='YlOrBr_r',
            vmin=0, vmax=1)
        CS = plt.contour(X_h, Y_h, h[1,:,:], np.arange(0,5000,500),
            colors='k')

        plt.pcolormesh(X_v, Y_v, np.ma.masked_where(mask_v==0, v[1,:,:]*100.),
            cmap='RdBu_r', vmin = -15, vmax = 15)
        CB = plt.colorbar()
        CB.set_label('y component of velocity (cm / s)')
        plt.xlim(0,1500)
        plt.axes().set_aspect('equal')
        plt.xlabel('x coordinate (km)')
        plt.ylabel('y coordinate (km)')
        plt.title('Timestep {}: Day {:05d}'.format(v_files[i][-10:], i))

        plt.savefig('{}figures/state_L1_{:06d}.png'.format(simulation,i), dpi=300,
            bbox_inches='tight')
        plt.close()

    try:
        sub.check_call(["ffmpeg", "-pattern_type", "glob", "-i",
            "{0}/figures/state_L0_*.png".format(simulation),
            "{0}/{1}.mp4".format(simulation, 'anim_L0')])
    except:
        print("failed to make L0 animation")

    try:
        sub.check_call(["ffmpeg", "-pattern_type", "glob", "-i",
            "{0}/figures/state_L1_*.png".format(simulation),
            "{0}/{1}.mp4".format(simulation, 'anim_L1')])
    except:
        print("failed to make L1 animation")



if __name__ == '__main__':
    sub.check_call(["mkdir", "-p", "spin_up_upwind/figures/"])

    plt_h_cross_section('spin_up_upwind/')
    plt_state('spin_up_upwind/')
