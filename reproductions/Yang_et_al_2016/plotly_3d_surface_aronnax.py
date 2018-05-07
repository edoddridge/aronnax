import numpy as np
import pandas as pd
import plotly.offline as pyo
import plotly.graph_objs as go

import os.path as p
import glob

import aronnax as aro


def wetmask(X, Y):
    """The wet mask."""

    # start with land everywhere and carve out space for water
    wetmask = np.zeros(X.shape, dtype=np.float64)
    # circular gyre region
    wetmask[((Y-500e3)**2 + (X-500e3)**2) < 500e3**2] = 1

    # clean up the edges
    wetmask[ 0, :] = 0
    wetmask[-1, :] = 0
    wetmask[ :, 0] = 0
    wetmask[ :,-1] = 0

    return wetmask

def bathymetry(X,Y):
    mask = wetmask(X, Y)
    r = np.sqrt((Y-500e3)**2 + (X-500e3)**2)

    # creating slope near southern boundary
    depth = np.minimum(4000*np.ones(X.shape), 15250 - 0.03*r)
    # set minimum depth to 250 m (also the model won't run if depthFile
    #   contains negative numbers)
    depth = np.maximum(depth, 250)
    depth[r>=499e3] = 0
    depth = np.ma.masked_where(r>507e3, depth)


    return depth

nx = 200
ny = 200
layers = 2
xlen = 1000e3
ylen = 1000e3
dx = xlen / nx
dy = ylen / ny

grid = aro.Grid(nx, ny, layers, dx, dy)

X_h, Y_h = np.meshgrid(grid.x, grid.y)
mask_h = wetmask(X_h, Y_h)
X_v, Y_v = np.meshgrid(grid.x, grid.yp1)
mask_v = wetmask(X_v, Y_v)

depth = bathymetry(X_h, Y_h)

file_path = "/Users/doddridge/Documents/Edward/Code/aronnax/reproductions/Yang_et_al_2016/spin_up/output/"

h_files = sorted(glob.glob(p.join(file_path, 'snap.h.*')))
eta_files = sorted(glob.glob(p.join(file_path, 'snap.eta.*')))
u_files = sorted(glob.glob(p.join(file_path, 'snap.u.*')))
v_files = sorted(glob.glob(p.join(file_path, 'snap.v.*')))



h_final = aro.interpret_raw_file(h_files[-1], nx, ny, layers)
h_masked = np.ma.masked_where(mask_h==0, h_final[0,:,:])
u_final = aro.interpret_raw_file(u_files[-1], nx, ny, layers)
v_final = aro.interpret_raw_file(v_files[-1], nx, ny, layers)

speed = np.sqrt((u_final[0,:,:-1]+u_final[0,:,1:])**2 + (v_final[:,1:,0] + v_final[:,:-1,0])**2)
speed[mask_h==0] = np.nan

eta_final = aro.interpret_raw_file(eta_files[-1], nx, ny, 1)
eta_masked = np.ma.masked_where(mask_h==0, eta_final[0,:,:])*1e2
eta_masked = eta_masked - np.min(eta_masked)

data = [go.Surface(x=grid.x/1e3, y=grid.y/1e3, z=-h_masked,
                    surfacecolor=speed,
                    colorscale='Viridis',
                    colorbar=dict(title='Speed (m/s)', titlefont={'size':30},
                        tickfont={'size':18}, x=1),
                    name='Depth of interface'),
        go.Surface(x=grid.x[:100]/1e3, y=grid.y/1e3, z=-depth[:,:100], surfacecolor=-depth[:,:100], 
                    colorscale='inferno',
                    # colorbar=dict(title='Bathymetry (m)', titlefont={'size':30},
                    #     tickfont={'size':18}, x=-0.1),
                    name='Bathymetry', showscale = False)]

layout = go.Layout(font=dict(family='Times New Roman'),
                    xaxis = dict(
                        title='x axis', titlefont={'size':20}),
                    yaxis=dict(title='Y axis'),
                    title='Layer interface depth after 14 years',
                    titlefont=dict(size=40, family='Times New Roman'),
                    )

fig = go.Figure(data=data, layout=layout)
pyo.plot(fig, filename='spin_up/final_state.html')
