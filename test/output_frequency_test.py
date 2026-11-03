'''Tests for validation of the output frequency options.'''

import pytest

import aronnax.driver as drv


def make_config(**options):
    config = drv.default_configuration()
    drv.merge_config(config, options)
    return config


def test_valid_output_frequencies():
    '''Zero (disabled) and multiples of dt are accepted.'''
    config = make_config(dt=600., dump_freq=6000., av_freq=0.,
                         checkpoint_freq=600., diag_freq=1.2e5)
    drv.validate_output_frequencies(config)


def test_no_frequencies_set():
    '''Options that are not set are left for the core's defaults.'''
    drv.validate_output_frequencies(make_config(dt=600.))
    drv.validate_output_frequencies(make_config())


@pytest.mark.parametrize("name", drv.output_frequency_options)
@pytest.mark.parametrize("value", [-600., 100.])
def test_invalid_output_frequencies(name, value):
    '''Negative frequencies and ones shorter than dt are rejected.'''
    config = make_config(dt=600., **{name: value})
    with pytest.raises(ValueError, match=name):
        drv.validate_output_frequencies(config)


def test_simulate_rejects_invalid_frequency_before_running(tmp_path):
    '''simulate fails fast, without compiling or running the core.'''
    with pytest.raises(ValueError, match="dump_freq"):
        drv.simulate(work_dir=str(tmp_path), dt=600., dump_freq=100.)


def test_merge_config_keeps_numbers_and_converts_booleans():
    '''0 and 1 must stay numeric (issue #132), while booleans are
    written as yes/no.'''
    config = make_config(av_freq=0., checkpoint_freq=1, red_grav=True,
                         dump_wind=False)
    assert config.getfloat("numerics", "av_freq") == 0.
    assert config.getfloat("numerics", "checkpoint_freq") == 1.
    assert config.get("model", "red_grav") == "yes"
    assert config.get("external_forcing", "dump_wind") == "no"
