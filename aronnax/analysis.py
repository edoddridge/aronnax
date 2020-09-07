"""
Analysis functions for Aronnax.

"""

import xarray as xr

def calc_zeta (u,v, grid, xgcm_grid):
    """
    Calculate relative vorticity from given velocity fields.
    Located at the vorticity point.
    Velocity fields must be generated with the 'open_mfdataarray' function.
    """

    zeta = (xgcm_grid.diff(v, 'X', boundary='extend')/grid.dx
        - xgcm_grid.diff(u, 'Y', boundary='extend')/grid.dy)

    zeta.name = 'zeta'

    return zeta


def calc_KE (u,v,h, rho0=1028):
    """
    Calculate kinetic energy at each grid point. Located at the tracer point.
    Velocity fields must be generated with the 'open_mfdataarray' function.
    """

    u_sq = (u[:,:,:,1:].rename({'xp1':'x'}).assign_coords(x=h['x'])**2 +
           u[:,:,:,:-1].rename({'xp1':'x'}).assign_coords(x=h['x'])**2)/2

    v_sq = (v[:,:,1:,:].rename({'yp1':'y'}).assign_coords(y=h['y'])**2 +
            v[:,:,:-1,:].rename({'yp1':'y'}).assign_coords(y=h['y'])**2)/2

    KE = (u_sq + v_sq)*h*rho0/2

    KE.name = 'KE'

    return KE

def calc_streamfunction(h, u, xgcm_grid):
    """
    Calculate the streamfunction. Located at the vorticity point.
    """

    h_on_u = xgcm_grid.interp(h, 'X', boundary='extend')

    psi = xgcm_grid.cumint(h_on_u*u, 'Y', to='outer', boundary='extend')


    psi.name = 'psi'
    return psi

def calc_speed (u,v):
    """
    Calculate flow speed at each grid point.
    Velocity fields must be generated with the 'open_mfdataarray' function.
    """

    u_sq = (u[:,:,:,1:].rename({'xp1':'x'}).assign_coords(x=v['x'])**2 +
           u[:,:,:,:-1].rename({'xp1':'x'}).assign_coords(x=v['x'])**2)/2

    v_sq = (v[:,:,1:,:].rename({'yp1':'y'}).assign_coords(y=u['y'])**2 +
            v[:,:,:-1,:].rename({'yp1':'y'}).assign_coords(y=u['y'])**2)/2

    speed = xr.ufuncs.sqrt(u_sq + v_sq)

    speed.name = 'speed'

    return speed
