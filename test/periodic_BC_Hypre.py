import os
import os.path as p
import subprocess as sub
import time

import glob

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

import output_preservation_test as opt

self_path = p.dirname(p.abspath(__file__))

def test_periodic_BC_Hypre():
    nx = 50
    ny = 20
    layers = 2

    dx = 5e4
    dy = 5e4

    grid = aro.Grid(nx, ny, layers, dx, dy)

    rho0 = 1035

    def wetmask(X, Y):
        mask = np.ones(X.shape, dtype=np.float64)
        return mask

    def eta(X, Y):
        return 0. + 1.*np.exp(-((6e5-X)**2 + (5e5-Y)**2)/(2*1e5**2))

    with working_directory(p.join(self_path, "periodic_BC")):
        drv.simulate(initHfile=[500.,500.],
                     nx=nx, ny=ny, layers=layers, dx=dx, dy=dy,
                     exe="aronnax_external_solver_test",
                     wetMaskFile=wetmask, 
                     fUfile=-1e-4,
                     fVfile=-1e-4,
                     initEtaFile=eta,
                     depthFile=1000., nTimeSteps=801,
                     dumpFreq=10000)
        opt.assert_outputs_close(nx, ny, layers, 3e-12)
        opt.assert_volume_conservation(nx, ny, layers, 3e-5)

if __name__ == '__main__':
    test_periodic_BC_Hypre()