from __future__ import print_function

import os.path as p
import sys

import numpy as np

from aronnax.utils import working_directory
import aronnax.driver as aro

self_path = p.dirname(p.abspath(__file__))

def bump(X, Y):
    return 500. + 20*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

def do_red_grav(nx, aro_build, perf):
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        aro.simulate(exe=aro_build, init_h_file=[bump], nx=nx, ny=nx, perf=perf)

def do_n_layer(nx, aro_build, perf):
    with working_directory(p.join(self_path, "beta_plane_bump")):
        aro.simulate(exe=aro_build, nx=nx, ny=nx, perf=perf,
                     init_h_file=[bump, lambda X, Y: 2000. - bump(X, Y)])

def main():
    aro_build = "aronnax_prof"
    nx = 100
    red_grav = True
    perf = False
    for arg in sys.argv[1:]:
        if arg.lower() == "red-grav":
            red_grav = True
        elif arg.lower() == "n-layer":
            red_grav = False
        elif arg.lower() == "perf":
            perf = True
            aro_build = "aronnax_core"
        else:
            nx = int(arg)
    if red_grav:
        print("Profiling a reduced gravity configuration" \
                    + " of size %dx%dx1 by 502 time steps." % (nx, nx))
        do_red_grav(nx, aro_build, perf)
        outpath = p.join(self_path, "beta_plane_bump_red_grav", "gmon.out")
    else:
        print("Profiling an n-layer configuration" \
                    + " of size %dx%dx2 by 502 time steps." % (nx, nx))
        do_n_layer(nx, aro_build, perf)
        outpath = p.join(self_path, "beta_plane_bump", "gmon.out")
    if not perf:
        print("Inspect %s,\nfor example with gprof aronnax_prof %s" % (outpath, outpath))

if __name__ == '__main__':
    main()
