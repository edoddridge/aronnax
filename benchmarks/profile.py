import os.path as p

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt

def do_red_grav(nx):
    with opt.working_directory(p.join(self_path, "beta_plane_bump_red_grav")):
        opt.run_experiment(
            opt.write_input_beta_plane_bump_red_grav, nx, nx, 1, "MIM_prof")

def main():
    nx = 100
    if len(sys.argv) > 1:
        nx = int(sys.argv[1])
    do_red_grav(nx)
    outpath = p.join(self_path, "beta_plane_bump_red_grav", "gmon.out")
    print "Inspect %s,\nfor example with gprof MIM_prof %s" % (outpath, outpath)

if __name__ == '__main__':
    main()
