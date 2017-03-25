import os.path as p

import glob

import numpy as np
import matplotlib.pyplot as plt

import MIMutils as mim

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

import sys
sys.path.append(p.join(root_path, 'test'))
import output_preservation_test as opt