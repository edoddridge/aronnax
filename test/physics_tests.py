import os.path as p

import glob

import numpy as np
import matplotlib.pyplot as plt

import MIMutils as mim

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt


def write_f_plane_red_grav_wind_input(nx,ny,layers):
    opt.write_f_plane(nx, ny, 10e-4)
    opt.write_rectangular_pool(nx, ny)
    with opt.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)
    with opt.fortran_file('wind_x.bin','w') as f:
        wind_x = np.zeros((ny,nx+1),dtype=np.float64)
        wind_x[50,50] = 0.2
        f.write_record(wind_x)

        plt.figure()
        plt.pcolormesh(wind_x)
        plt.colorbar()
        plt.savefig('wind_x.png')

def write_f_plane_red_grav_init_u_input(nx,ny,layers):
    opt.write_f_plane(nx, ny, 10e-4)
    opt.write_rectangular_pool(nx, ny)
    with opt.fortran_file('initH.bin', 'w') as f:
        f.write_record(np.ones((nx, ny,layers), dtype=np.float64) * 400)
    with opt.fortran_file('init_u.bin','w') as f:
        init_u = np.zeros((ny,nx+1),dtype=np.float64)
        init_u[50,50] = 0.2
        f.write_record(init_u)

        plt.figure()
        plt.pcolormesh(init_u)
        plt.colorbar()
        plt.savefig('init_u.png')