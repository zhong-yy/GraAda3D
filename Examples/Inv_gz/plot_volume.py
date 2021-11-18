#https://plotly.com/python/3d-volume-plots/#reference
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import RBFInterpolator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import plotly.graph_objects as go

x,y,z,v=np.loadtxt("gz_result.txt",skiprows=2,usecols=[6,7,8,9],unpack=True)

xyz=np.hstack((x[:,np.newaxis],y[:,np.newaxis],z[:,np.newaxis]))

xi=np.linspace(0,2000,41)
yi=np.linspace(0,2000,41)
zi=np.linspace(0,1000,21)
XI,YI,ZI=np.meshgrid(xi,yi,zi)
#VI=griddata((x,y,z),v,(XI,YI,ZI),method='linear')

rbf=RBFInterpolator(xyz,v)

XI_flat=XI.flatten()
YI_flat=YI.flatten()
ZI_flat=ZI.flatten()

XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
VI_flat=rbf(XYZI_flat)

size=XI.shape

VI=np.reshape(VI_flat,size)

#fig = go.Figure(data=go.Volume(
#    x=XI.flatten(),
#    y=YI.flatten(),
#    z=ZI.flatten(),
#    value=VI.flatten(),
#    isomin=-300,
#    isomax=300,
#    opacity=0.5, # needs to be small to see through all surfaces
#    surface_count=11, # needs to be a large number for good volume rendering
#    slices_z=dict(show=True, locations=[350]),
#    ))
#fig.show()

fig = go.Figure(data=go.Volume(
    x=XI.flatten(),
    y=YI.flatten(),
    z=ZI.flatten(),
    value=VI.flatten(),
    slices_z=dict(show=True, locations=[350]),
    slices_y=dict(show=True, locations=[700,1450]),
    caps= dict(x_show=False, y_show=False, z_show=False),
    colorscale='jet'
    ))
fig.update_scenes(zaxis_autorange="reversed")
fig.show()
