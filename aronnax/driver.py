import ConfigParser as par
import os
import os.path as p
import subprocess as sub

from aronnax.core import fortran_file
from aronnax.core import interpret_requested_data
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

def simulate(work_dir=".", config_path="aronnax.conf", **options):
    config_file = p.join(work_dir, config_path)
    config = par.RawConfigParser()
    config.read(config_file)
    merge_config(config, options)
    with working_directory(work_dir):
        compile_core(config)
        # XXX Try to avoid overwriting the input configuration
        with open('aronnax-merged.conf', 'w') as f:
            config.write(f)
        sub.check_call(["rm", "-rf", "input/"])
        sub.check_call(["rm", "-rf", "output/"])
        sub.check_call(["mkdir", "-p", "output/"])
        with working_directory("input"):
            generate_input_data_files(config)
        generate_parameters_file(config)
        run_executable(config)
        sub.check_call(["rm", "-rf", "netcdf-output/"])
        sub.check_call(["mkdir", "-p", "netcdf-output/"])
        convert_output_to_netcdf(config)

section_map = {
    "au"                   : "numerics",
    "ah"                   : "numerics",
    "ar"                   : "numerics",
    "dt"                   : "numerics",
    "slip"                 : "numerics",
    "nTimeSteps"           : "numerics",
    "dumpFreq"             : "numerics",
    "avFreq"               : "numerics",
    "hmin"                 : "numerics",
    "maxits"               : "numerics",
    "eps"                  : "numerics",
    "freesurfFac"          : "numerics",
    "thickness_error"      : "numerics",
    "hmean"                : "model",
    "depthFile"            : "model",
    "H0"                   : "model",
    "RedGrav"              : "model",
    "spongeHTimeScaleFile" : "sponge",
    "spongeUTimeScaleFile" : "sponge",
    "spongeVTimeScaleFile" : "sponge",
    "spongeHFile"          : "sponge",
    "spongeUFile"          : "sponge",
    "spongeVFile"          : "sponge",
    "g_vec"                : "physics",
    "rho0"                 : "physics",
    "nx"                   : "grid",
    "ny"                   : "grid",
    "layers"               : "grid",
    "dx"                   : "grid",
    "dy"                   : "grid",
    "fUfile"               : "grid",
    "fVfile"               : "grid",
    "wetMaskFile"          : "grid",
    "initUfile"            : "initial_conditions",
    "initVfile"            : "initial_conditions",
    "initHfile"            : "initial_conditions",
    "initEtaFile"          : "initial_conditions",
    "zonalWindFile"        : "external_forcing",
    "meridionalWindFile"   : "external_forcing",
    "DumpWind"             : "external_forcing",
    "wind_mag_time_series_file" : "external_forcing",
    "exe"                  : "executable",
    "valgrind"             : "executable",
    "perf"                 : "executable",
}

def merge_config(config, options):
    for k, v in options.iteritems():
        if k in section_map:
            section = section_map[k]
            if not config.has_section(section):
                config.add_section(section)
            config.set(section, k, str(v))
        else:
            raise Exception("Unrecognized option", k)

def compile_core(config):
    core_name = config.get("executable", "exe")
    with working_directory(root_path):
        sub.check_call(["make", core_name])

def generate_input_data_files(config):
    for name, section in section_map.iteritems():
        if not name.endswith("File") and not name.endswith("file"):
            continue
        requested_data = config.get(section, name)
        generated_data = interpret_requested_data(requested_data, config)
        if generated_data is not None:
            with fortran_file(name + '.bin', 'w') as f:
                f.write_record(generated_data)

def generate_parameters_file(config):
    with open('parameters.in', 'w') as f:
        for section in config.sections():
            f.write(' &')
            f.write(section.toupper())
            f.write('\n')
            for (name, value) in config.items(section):
                f.write(' %s = %s,\n' % (name, value))
            f.write(' /\n')

def run_executable(config):
    core_name = config.get("executable", "exec")
    env = dict(os.environ, GFORTRAN_STDERR_UNIT="17")
    if config.getboolean("executable", "valgrind") \
       or 'ARONNAX_TEST_VALGRIND_ALL' in os.environ:
        assert not config.getboolean("executable", "perf")
        sub.check_call(["valgrind", "--error-exitcode=5", p.join(root_path, aro_exec)],
            env=env)
    elif config.getboolean("executable", "perf"):
        perf_cmds = ["perf", "stat", "-e", "r530010", # "flops", on my CPU.
            "-e", "L1-dcache-loads", "-e", "L1-dcache-load-misses",
            "-e", "L1-dcache-stores", "-e", "L1-dcache-store-misses",
            "-e", "L1-icache-loads", "-e", "L1-icache-misses",
            "-e", "L1-dcache-prefetches",
            "-e", "branch-instructions", "-e", "branch-misses"]
        sub.check_call(perf_cmds + [p.join(root_path, aro_exec)], env=env)
    else:
        sub.check_call([p.join(root_path, aro_exec)], env=env)

def convert_output_to_netcdf(config):
    # TODO Issue #30
    pass
