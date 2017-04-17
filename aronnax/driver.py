import ConfigParser as par
import os.path as p
import subprocess as sub

from aronnax.core import fortran_file
from aronnax.utils import working_directory

def simulate(work_dir=".", config_path="aronnax.conf", **options):
    config_file = p.join(work_dir, config_path)
    config = par.RawConfigParser()
    config.read(config_file)
    merge_config(config, options)
    with working_directory(config.get("...", "core_work_dir")):
        compile_core(config)
        # XXX Try to avoid overwriting the input configuration
        with open('aronnax.conf', 'w') as f:
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

section_map = dict(
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
    )

def merge_config(config, options):
    for k, v in options.iteritems():
        if k in section_map:
            config.set(section_map[k], k, str(v))
        else:
            raise Exception("Unrecognized option", k)

def generate_data_files(config):
    for name, section in section_map.iteritems():
        if not name.endswith("File") and not name.endswith("file"):
            continue
        requested_data = config.get(section, name)
        generated_data = interpret_requested_data(requested_data, config)
        if generated_data is not None:
            with fortran_file(name + '.bin', 'w') as f:
                f.write_record(generated_data)

def convert_output_to_netcdf(config):
    # TODO Issue #30
    pass
