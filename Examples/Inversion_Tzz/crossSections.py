from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
#from xyz2XYZ import xyz2XYZ
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

mydata=Dataset('Tzz_result.nc', "r")# r mean read
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


z_id=get_slice_index(xs,350)
z_slice=density[z_id,:,:]

y_id=get_slice_index(ys,700)
y_slice1=density[:,y_id,:]

y_id=get_slice_index(ys,1450)
y_slice2=density[:,y_id,:]


fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,6))

Y,Z=np.meshgrid(ys,zs)
pc=ax1.pcolormesh(Y,Z,y_slice1,cmap='jet')
ax1.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax1.set_title(r'Y=700m')
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_yaxis()
clb=fig.colorbar(pc,ax=ax1,shrink=0.8,orientation='horizontal')
clb_title=clb.ax.set_title(r'kg/m$^3$')
ax1.set_xlabel("y (m)")
ax1.set_ylabel("z (m)")

pc=ax2.pcolormesh(Y,Z,y_slice2,cmap='jet')
ax2.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax2.set_title(r'Y=1450m')
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.invert_yaxis()
clb=fig.colorbar(pc,ax=ax2,shrink=0.8,orientation='horizontal')
clb_title=clb.ax.set_title(r'kg/m$^3$')
ax2.set_xlabel("y (m)")
ax2.set_ylabel("z (m)")
plt.savefig('inversion_model_Y_slices.jpg',dpi=300,bbox_inches='tight')


fig,ax3=plt.subplots(1,1,figsize=(6,4))
X,Y=np.meshgrid(xs,ys)
pc=ax3.pcolormesh(X,Y,z_slice,cmap='jet')
ax3.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax3.set_title(r'Z=350m')
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.invert_yaxis()
clb=fig.colorbar(pc,ax=ax3,shrink=1)
clb_title=clb.ax.set_title(r'kg/m$^3$')
ax2.set_xlabel("x (m)")
ax2.set_ylabel("y (m)")
plt.savefig('inversion_model_Z_slices.jpg',dpi=300,bbox_inches='tight')

