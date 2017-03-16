import os.path as p

import numpy as np
import matplotlib.pyplot as plt

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
sys.path.append(p.join(root_path, 'reproductions/Davis_et_al_2014'))

import output_preservation_test as opt

import subprocess as sub

import MIMutils as mim

import time

mim_exec = p.join(root_path, "MIM")


def run_davis_et_al_2014(nx,ny,layers,nTimeSteps,dt):
    with opt.working_directory(p.join(self_path, "Davis_et_al_2014")):
        run_experiment(write_input_davis_et_al_2014, nx, ny, layers, nTimeSteps, dt)


# Generic function to run an experiment
def run_experiment(write_input, nx, ny, layers, nTimeSteps,dt):
    """ Run MIM"""
    sub.check_call(["rm", "-rf", "input/"])
    sub.check_call(["rm", "-rf", "output/"])
    sub.check_call(["mkdir", "-p", "output/"])
    sub.check_call(["rm", "-rf", "run_finished.txt"])

    with opt.working_directory(root_path):
        sub.check_call(["make", "MIM"])
    with opt.working_directory("input"):
        write_input(nx, ny, layers,nTimeSteps,dt)
    tweak_parameters(nx, ny, layers,nTimeSteps, dt)

    then = time.time()
    sub.check_call([mim_exec])
    print "MIM execution took", time.time() - then
