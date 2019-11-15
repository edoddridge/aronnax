"""
Analysis functions for Aronnax.

"""

def calc_zeta (u,v):
    """
    Calculate relative vorticity from given velocity fields.
    Velocity fields must be generated with the 'open_mfdataarray' function.
    """
    dvdy = v.diff(dim='x')[:,:,1:-1,:]/(v.x[2] - v.x[1])
    dvdy = dvdy.rename({'x': 'xp1'})
    dvdy = dvdy.assign_coords(xp1= u['xp1'][1:-1])

    dudx = u.diff(dim='y')[:,:,:,1:-1]
    dudx = dudx.rename({'y': 'yp1'})
    dudx = dudx.assign_coords(yp1= v['yp1'][1:-1])/(u.y[2] - u.y[1])

    zeta = dvdy - dudx

    zeta.name = 'zeta'

    return zeta


def calc_KE (u,v,h, rho0=1028):
    """
    Calculate kinetic energy at each grid point.
    Velocity fields must be generated with the 'open_mfdataarray' function.
    """

    u_sq = (u[:,:,:,1:].rename({'xp1':'x'}).assign_coords(x=h['x'])**2 +
           u[:,:,:,:-1].rename({'xp1':'x'}).assign_coords(x=h['x'])**2)/2

    v_sq = (v[:,:,1:,:].rename({'yp1':'y'}).assign_coords(y=h['y'])**2 +
            v[:,:,:-1,:].rename({'yp1':'y'}).assign_coords(y=h['y'])**2)/2

    KE = (u_sq + v_sq)*h*rho0/2

    return KE


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

    return speed
