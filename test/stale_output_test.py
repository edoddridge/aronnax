'''Tests that results from earlier runs are not mixed with new ones.'''

import os
import os.path as p

import aronnax.driver as drv


def touch(path):
    with open(path, 'w') as f:
        f.write('x')


def make_run_dir(tmp_path):
    names = {
        'output': ['snap.h.0000000001', 'snap.h.0000000101',
                   'snap.h.0000000201', 'av.u.0000000301',
                   'diagnostic.h.csv', 'notes.txt'],
        'checkpoints': ['h.0000000001', 'h.0000000101', 'h.0000000201'],
    }
    for directory, files in names.items():
        os.makedirs(str(tmp_path / directory))
        for name in files:
            touch(str(tmp_path / directory / name))


def remaining(tmp_path, directory):
    return sorted(os.listdir(str(tmp_path / directory)))


def test_new_run_removes_all_numbered_files(tmp_path, monkeypatch):
    make_run_dir(tmp_path)
    monkeypatch.chdir(str(tmp_path))
    drv.remove_stale_output(drv.default_configuration())
    # Files not named by iteration are left alone.
    assert remaining(tmp_path, 'output') == ['diagnostic.h.csv', 'notes.txt']
    assert remaining(tmp_path, 'checkpoints') == []


def test_restart_keeps_files_up_to_niter0(tmp_path, monkeypatch):
    make_run_dir(tmp_path)
    monkeypatch.chdir(str(tmp_path))
    config = drv.default_configuration()
    drv.merge_config(config, {'niter0': 101})
    drv.remove_stale_output(config)
    assert remaining(tmp_path, 'output') == [
        'diagnostic.h.csv', 'notes.txt',
        'snap.h.0000000001', 'snap.h.0000000101']
    assert remaining(tmp_path, 'checkpoints') == [
        'h.0000000001', 'h.0000000101']


def test_missing_directories_are_ignored(tmp_path, monkeypatch):
    monkeypatch.chdir(str(tmp_path))
    drv.remove_stale_output(drv.default_configuration())
