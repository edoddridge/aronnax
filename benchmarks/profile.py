import os.path as p

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt
from aronnax.utils import working_directory

def do_red_grav(nx, aro_build, perf):
    with working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        opt.run_experiment(
            opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, aro_build, perf=perf)

def do_n_layer(nx, aro_build, perf):
    with working_directory(p.join(self_path, "beta_plane_bump")):
        opt.run_experiment(
            opt.write_input_beta_plane_bump, nx, nx, 2, aro_build, perf=perf)

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
        print "Profiling a reduced gravity configuration" \
            + " of size %dx%dx1 by 502 time steps." % (nx, nx)
        do_red_grav(nx, aro_build, perf)
        outpath = p.join(self_path, "beta_plane_bump_red_grav", "gmon.out")
    else:
        print "Profiling an n-layer configuration" \
            + " of size %dx%dx2 by 502 time steps." % (nx, nx)
        do_n_layer(nx, aro_build, perf)
        outpath = p.join(self_path, "beta_plane_bump", "gmon.out")
    if not perf:
        print "Inspect %s,\nfor example with gprof aronnax_prof %s" % (outpath, outpath)

if __name__ == '__main__':
    main()
