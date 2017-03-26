"""!
The main file with the class definitions

Core
==============

This file contains all of the classes for the module.
"""

import numpy as np

class Grid(object):
    """!Make a grid object containing all of the axes."""

    def __init__(self,nx,ny,dx,dy,x0=0,y0=0):
        """!Instantiate a grid object for MIM."""

        # axes for vorticity points
        self.xp1 = np.linspace(x0,nx*dx+x0,nx+1)
        self.yp1 = np.linspace(y0,ny*dy+y0,ny+1)

        # Axes for tracer points.
        self.x = (self.xp1[1:] + self.xp1[:-1])/2.
        self.y = (self.yp1[1:] + self.yp1[:-1])/2.
