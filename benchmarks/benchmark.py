import os
import os.path as p
import subprocess as sub
import time

import numpy as np
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test') )
import output_preservation_test as opt


grid_points = np.array([10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 500])

with opt.working_directory(root_path):
    sub.check_call(["make", "MIM_test"])
    sub.check_call(["make", "MIM"])

def run_experiment(write_input, nx, ny, layers, mim_exec, valgrind=False):
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])

    with opt.working_directory("input"):
        write_input(nx, ny, layers)
    opt.tweak_parameters(nx, ny, layers)
    then = time.time()
    if valgrind or 'MIM_TEST_VALGRIND_ALL' in os.environ:
        sub.check_call(["valgrind", mim_exec])
    else:
        sub.check_call([mim_exec])
    run_time = time.time() - then
    print "MIM execution took", run_time

    return run_time


def benchmark_gaussian_bump_red_grav():
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))

    with opt.working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        mim_exec = p.join(root_path, "MIM_test")
        for counter, nx in enumerate(grid_points):
            run_time_O1[counter] = run_experiment(
                  opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, mim_exec)

        mim_exec = p.join(root_path, "MIM")
        for counter, nx in enumerate(grid_points):
            run_time_Ofast[counter] = run_experiment(
                  opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, mim_exec)
        
        plt.figure()
        plt.plot(grid_points, run_time_O1, '-*', label='MIM run time -O1')
        plt.plot(grid_points, run_time_Ofast, '-*', label='MIM run time -Ofast')
        plt.plot(grid_points, (run_time_O1[-7]/(grid_points[-7]**2))*grid_points**2, '-*', label='O(nx*nx)')
        plt.legend()
        plt.xlabel('nx')
        plt.ylabel('run time (s)')
        plt.savefig('beta_plane_bump_red_grav scaling.png', dpi=150)


def benchmark_gaussian_bump():
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))

    with opt.working_directory(p.join(self_path, "beta_plane_bump")):
        mim_exec = p.join(root_path, "MIM_test")
        for counter, nx in enumerate(grid_points[:7]):
            run_time_O1[counter] = run_experiment(
                  opt.write_input_beta_plane_bump, nx, nx, 2, mim_exec)
        
        mim_exec = p.join(root_path, "MIM")
        for counter, nx in enumerate(grid_points[:7]):
            run_time_Ofast[counter] = run_experiment(
                  opt.write_input_beta_plane_bump, nx, nx, 2, mim_exec)
        
        plt.figure()
        plt.plot(grid_points[:7], run_time_O1[:7], '-*', label='MIM run time -O1')
        plt.plot(grid_points[:7], run_time_Ofast[:7], '-*',
            label='MIM run time -Ofast')
        plt.plot(grid_points[:7],
            (run_time_O1[-6]/(grid_points[-6]**2))*grid_points[:7]**2,
            '-*', label='O(nx**2)')
        plt.plot(grid_points[:7],
            (run_time_O1[-6]/(grid_points[-6]**3))*grid_points[:7]**3,
            '-*', label='O(nx**3)')
        plt.legend()
        plt.xlabel('nx')
        plt.ylabel('run time (s)')
        plt.savefig('beta_plane_bump scaling.png', dpi=150)

if __name__ == '__main__':
    benchmark_gaussian_bump_red_grav()
    benchmark_gaussian_bump()
