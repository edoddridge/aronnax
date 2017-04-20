import cPickle as pkl
import os.path as p
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aronnax.utils import working_directory
import aronnax.driver as aro

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

n_time_steps = 502.0
scale_factor = 1000 / n_time_steps # Show times in ms

def benchmark_gaussian_bump_red_grav_save(grid_points):
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))
    def bump(X, Y):
        return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        aro_exec = "aronnax_test"
        for counter, nx in enumerate(grid_points):
            run_time_O1[counter] = aro.simulate(
                exe=aro_exec, initHfile=[bump], nx=nx, ny=nx)

        aro_exec = "aronnax_core"
        for counter, nx in enumerate(grid_points):
            run_time_Ofast[counter] = aro.simulate(
                exe=aro_exec, initHfile=[bump], nx=nx, ny=nx)

        with open("times.pkl", "w") as f:
            pkl.dump((grid_points, run_time_O1, run_time_Ofast), f)

def benchmark_gaussian_bump_red_grav_plot():
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        with open("times.pkl", "r") as f:
            (grid_points, run_time_O1, run_time_Ofast) = pkl.load(f)

        plt.figure()
        plt.loglog(grid_points, run_time_O1*scale_factor,
            '-*', label='Aronnax run time -O1')
        plt.loglog(grid_points, run_time_Ofast*scale_factor,
            '-*', label='Aronnax run time -Ofast')
        scale = scale_factor * run_time_O1[-7]/(grid_points[-7]**2)
        plt.loglog(grid_points, scale*grid_points**2,
                   ':', label='O(nx**2)', color='black', linewidth=0.5)
        plt.legend()
        plt.xlabel('Resolution (grid cells on one side)')
        plt.ylabel('Avg time per integration step (ms)')
        plt.title('Runtime scaling of a 1.5-layer Aronnax simulation on a square grid')
        plt.savefig('beta_plane_bump_red_grav scaling.png', dpi=150)
        filename = p.join(root_path, 'docs/beta_plane_bump_red_grav_scaling.png')
        plt.savefig(filename, dpi=150)

def benchmark_gaussian_bump_red_grav(grid_points):
    benchmark_gaussian_bump_red_grav_save(grid_points)
    benchmark_gaussian_bump_red_grav_plot()


def benchmark_gaussian_bump_save(grid_points):
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))
    def bump(X, Y):
        return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

    with working_directory(p.join(self_path, "beta_plane_bump")):
        aro_exec = "aronnax_test"
        for counter, nx in enumerate(grid_points):
            run_time_O1[counter] = aro.simulate(
                exe=aro_exec, initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)], nx=nx, ny=nx)
        aro_exec = "aronnax_core"
        for counter, nx in enumerate(grid_points):
            run_time_Ofast[counter] = aro.simulate(
                exe=aro_exec, initHfile=[bump, lambda X, Y: 2000. - bump(X, Y)], nx=nx, ny=nx)
        with open("times.pkl", "w") as f:
            pkl.dump((grid_points, run_time_O1, run_time_Ofast), f)

def benchmark_gaussian_bump_plot():
    with working_directory(p.join(self_path, "beta_plane_bump")):
        with open("times.pkl", "r") as f:
            (grid_points, run_time_O1, run_time_Ofast) = pkl.load(f)

        plt.figure()
        plt.loglog(grid_points, run_time_O1*scale_factor,
            '-*', label='Aronnax run time -O1')
        plt.loglog(grid_points, run_time_Ofast*scale_factor,
            '-*', label='Aronnax run time -Ofast')
        scale = scale_factor * run_time_O1[3]/(grid_points[3]**3)
        plt.loglog(grid_points, scale*grid_points**3,
            ':', label='O(nx**3)', color='black', linewidth=0.5)
        plt.legend()
        plt.xlabel('Resolution (grid cells on one side)')
        plt.ylabel('Avg time per integration step (ms)')
        plt.title('Runtime scaling of a 2-layer Aronnax simulation\nwith bathymetry on a square grid')
        plt.savefig('beta_plane_bump scaling.png', dpi=150)
        filename = p.join(root_path, 'docs/beta_plane_bump_scaling.png')
        plt.savefig(filename, dpi=150)

def benchmark_gaussian_bump(grid_points):
    benchmark_gaussian_bump_save(grid_points)
    benchmark_gaussian_bump_plot()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == "save":
            benchmark_gaussian_bump_red_grav_save(np.array([10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 500]))
            benchmark_gaussian_bump_save(np.array([10, 20, 40, 60, 80, 100]))
        else:
            benchmark_gaussian_bump_red_grav_plot()
            benchmark_gaussian_bump_plot()
    else:
        benchmark_gaussian_bump_red_grav(np.array([10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 500]))
        benchmark_gaussian_bump(np.array([10, 20, 40, 60, 80, 100]))