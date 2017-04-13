import cPickle as pkl
import os.path as p
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt


grid_points = np.array([10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 500])

def benchmark_gaussian_bump_red_grav_save():
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))

    with opt.working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        aro_exec = "aronnax_test"
        for counter, nx in enumerate(grid_points):
            run_time_O1[counter] = opt.run_experiment(
                opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, aro_exec)

        aro_exec = "aronnax_core"
        for counter, nx in enumerate(grid_points):
            run_time_Ofast[counter] = opt.run_experiment(
                opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, aro_exec)

        with open("times.pkl", "w") as f:
            pkl.dump((run_time_O1, run_time_Ofast), f)

def benchmark_gaussian_bump_red_grav_plot():
    with opt.working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        with open("times.pkl", "r") as f:
            (run_time_O1, run_time_Ofast) = pkl.load(f)

        plt.figure()
        plt.plot(grid_points, run_time_O1, '-*', label='Aronnax run time -O1')
        plt.plot(grid_points, run_time_Ofast,
            '-*', label='Aronnax run time -Ofast')
        plt.plot(grid_points,
            (run_time_O1[-7]/(grid_points[-7]**2))*grid_points**2,
            '-*', label='O(nx**2)')
        plt.legend()
        plt.xlabel('nx')
        plt.ylabel('run time (s)')
        plt.savefig('beta_plane_bump_red_grav scaling.png', dpi=150)
        filename = p.join(root_path, 'docs/beta_plane_bump_red_grav_scaling.png')
        plt.savefig(filename, dpi=150)

def benchmark_gaussian_bump_red_grav():
    benchmark_gaussian_bump_red_grav_save()
    benchmark_gaussian_bump_red_grav_plot()


def benchmark_gaussian_bump_save():
    run_time_O1 = np.zeros(len(grid_points))
    run_time_Ofast = np.zeros(len(grid_points))

    with opt.working_directory(p.join(self_path, "beta_plane_bump")):
        aro_exec = "aronnax_test"
        for counter, nx in enumerate(grid_points[:6]):
            run_time_O1[counter] = opt.run_experiment(
                  opt.write_input_beta_plane_bump, nx, nx, 2, aro_exec)
        aro_exec = "aronnax_core"
        for counter, nx in enumerate(grid_points[:6]):
            run_time_Ofast[counter] = opt.run_experiment(
                  opt.write_input_beta_plane_bump, nx, nx, 2, aro_exec)
        with open("times.pkl", "w") as f:
            pkl.dump((run_time_O1, run_time_Ofast), f)

def benchmark_gaussian_bump_plot():
    with opt.working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        with open("times.pkl", "r") as f:
            (run_time_O1, run_time_Ofast) = pkl.load(f)

        plt.figure()
        plt.plot(grid_points[:6], run_time_O1[:6],
            '-*', label='Aronnax run time -O1')
        plt.plot(grid_points[:6], run_time_Ofast[:6], '-*',
            label='Aronnax run time -Ofast')
        plt.plot(grid_points[:6],
            (run_time_O1[3]/(grid_points[3]**2))*grid_points[:6]**2,
            '-*', label='O(nx**2)')
        plt.plot(grid_points[:6],
            (run_time_O1[3]/(grid_points[3]**3))*grid_points[:6]**3,
            '-*', label='O(nx**3)')
        plt.legend()
        plt.xlabel('nx')
        plt.ylabel('run time (s)')
        plt.savefig('beta_plane_bump scaling.png', dpi=150)
        filename = p.join(root_path, 'docs/beta_plane_bump_scaling.png')
        plt.savefig(filename, dpi=150)

def benchmark_gaussian_bump():
    benchmark_gaussian_bump_save()
    benchmark_gaussian_bump_plot()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == "save":
            benchmark_gaussian_bump_red_grav_save()
            benchmark_gaussian_bump_save()
        else:
            benchmark_gaussian_bump_red_grav_plot()
            benchmark_gaussian_bump_plot()
    else:
        benchmark_gaussian_bump_red_grav()
        benchmark_gaussian_bump()
