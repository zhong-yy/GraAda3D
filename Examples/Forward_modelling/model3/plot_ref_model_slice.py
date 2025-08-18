from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
#from xyz2XYZ import xyz2XYZ
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

mydata=Dataset('model/reference.readable.nc', "r")# r mean read
print(mydata)
var=mydata.variables
x=var['x']
y=var['y']
z=var['z']
density=var['density']

def node_to_pixel(x):
    dx=(0.5*(x[1]-x[0]))
    grid_line=x-dx
    grid_line=np.concatenate((grid_line,x[-1]+dx),axis=None)
    return grid_line

def get_slice_index(xs_pixel,x0):
    index_slice=-1
    for i in range(0,xs_pixel.size-1):
        if (x0>xs_pixel[i] or np.abs(x0-xs_pixel[i])<1e-15) and (x0<xs_pixel[i+1] or np.abs(x0-xs_pixel[i+1])<1e-15):
            index_slice=i
            print(xs_pixel[i:i+2])
            break
    return index_slice

xs=node_to_pixel(x)
ys=node_to_pixel(y)
zs=node_to_pixel(z)


z_id=get_slice_index(zs,400)
z_slice=density[z_id,:,:]

y_id=get_slice_index(ys,1000)
y_slice1=density[:,y_id,:]

min_value_shown=-500
max_value_shown=500
fts=12
clrmap='jet'
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(5,8),constrained_layout=True,height_ratios=[1,2])
# plt.subplots_adjust(hspace=0.1)

X,Z=np.meshgrid(xs,zs)
pc=ax1.pcolormesh(X,Z,y_slice1,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown)
ax1.set_aspect('equal')#相当于matlab中的axis equal
#ax1.set_title(r'Y=700m')
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_yaxis()
#clb=fig.colorbar(pc,ax=ax1,shrink=0.8,orientation='horizontal')
#clb_title=clb.ax.set_title(r'kg/m$^3$')
ax1.set_xlabel("x (m)",fontsize=fts)
ax1.set_ylabel("z (m)",fontsize=fts)
ax1.text(0.01,0.88,"(a) Y=1000 m",fontsize=fts, transform=ax1.transAxes)


X,Y=np.meshgrid(xs,ys)
pc=ax2.pcolormesh(X,Y,z_slice,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown)
ax2.set_aspect('equal')#相当于matlab中的axis equal
#ax2.set_title(r'Z=350m')
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.invert_yaxis()

clb=fig.colorbar(pc,ax=ax2,shrink=0.9,orientation='horizontal')
clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts)
clb.ax.tick_params(labelsize=12.5)

ax2.set_xlabel("x (m)",fontsize=fts)
ax2.set_ylabel("y (m)",fontsize=fts)
ax2.text(0.01,0.92,"(b) Z=400 m",fontsize=fts, transform=ax2.transAxes)
#plt.savefig('inversion_model_Z_slices.jpg',dpi=300,bbox_inches='tight')
for ax in fig.axes:
    ax.tick_params(which='both',labelsize=fts)
print(np.max(z_slice))
print(np.min(z_slice))
plt.savefig('reference_model.jpg',dpi=300,bbox_inches='tight')


