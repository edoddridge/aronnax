'''Tests for the helper that compares outputs with the blessed outputs.'''

import pytest

import output_preservation_test as opt


def make_dirs(tmp_path, out_names, good_names):
    for directory, names in (('output', out_names),
                             ('good-output', good_names)):
        (tmp_path / directory).mkdir()
        for name in names:
            (tmp_path / directory / name).write_bytes(b'')


def test_missing_output_file_is_detected(tmp_path, monkeypatch):
    '''A run that writes fewer files than the blessed run must fail.'''
    make_dirs(tmp_path, ['snap.h.0000000001'],
              ['snap.h.0000000001', 'snap.h.0000000101'])
    monkeypatch.chdir(str(tmp_path))
    with pytest.raises(AssertionError, match='snap.h.0000000101'):
        opt.assert_outputs_close(10, 10, 1, 1e-15)


def test_unexpected_output_file_is_detected(tmp_path, monkeypatch):
    '''A run that writes extra files must fail too.'''
    make_dirs(tmp_path, ['snap.h.0000000001', 'snap.u.0000000001'],
              ['snap.h.0000000001'])
    monkeypatch.chdir(str(tmp_path))
    with pytest.raises(AssertionError, match='snap.u.0000000001'):
        opt.assert_outputs_close(10, 10, 1, 1e-15)
