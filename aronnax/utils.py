from contextlib import contextmanager
import os
import subprocess as sub

@contextmanager
def working_directory(path):
    old_path = os.getcwd()
    sub.check_call(["mkdir", "-p", path])
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_path)

