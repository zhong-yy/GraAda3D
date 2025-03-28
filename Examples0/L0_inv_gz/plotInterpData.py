import sys
import os
include_path=os.path.abspath('../../PythonPlot')
sys.path.append(include_path)
#own module
import interpData

#plot interpolated data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

#own module
import interpData

x,y,z,density=np.loadtxt('gz_result.txt',skiprows=2,usecols=[6,7,8,9],unpack=True)
print(np.min(density),np.max(density))

fig=plt.figure(figsize=(12,6))
gs0 = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1])

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0],hspace=0.4)

min_value_shown=-200
max_value_shown=200
fts=12
clrmap='jet'

model='gz_result.txt'
#Y=700
Xi1,Zi1,Vi1=interpData.interp_yslice(model,y0=700,
                         start_x=0,stop_x=2000,num_x=81,
                         start_z=0,stop_z=1000,num_z=41,neighbors=100,kernel='linear')
                         
ax1=fig.add_subplot(gs00[0])
pc=ax1.pcolormesh(Xi1,Zi1,Vi1,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax1.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_yaxis()
ax1.set_xlabel("x (m)",fontsize=fts)
ax1.set_ylabel("z (m)",fontsize=fts)
ax1.text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=ax1.transAxes)

#Y=1450
Xi2,Zi2,Vi2=interpData.interp_yslice(model,y0=1450,
                         start_x=0,stop_x=2000,num_x=81,
                         start_z=0,stop_z=1000,num_z=41,neighbors=100,kernel='linear')
ax2=fig.add_subplot(gs00[1])
pc=ax2.pcolormesh(Xi2,Zi2,Vi2,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax2.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.invert_yaxis()
ax2.set_xlabel("x (m)",fontsize=fts)
ax2.set_ylabel("z (m)",fontsize=fts)
ax2.text(0.01,0.9,"(b) Y=1450 m",fontsize=fts, transform=ax2.transAxes)


#Z=350
Xi3,Yi3,Vi3=interpData.interp_zslice(model,z0=350,
                         start_x=0,stop_x=2000,num_x=81,
                         start_y=0,stop_y=2000,num_y=81,neighbors=100,kernel='linear')
gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],height_ratios=[1,7,1])
ax3=fig.add_subplot(gs01[1,0])
pc=ax3.pcolormesh(Xi3,Yi3,Vi3,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax3.set_aspect('equal', 'box')#相当于matlab中的axis equal
#ax3.set_title(r'Z=350m')
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.invert_yaxis()
clb=fig.colorbar(pc,ax=ax3,shrink=1)
clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts)
clb.ax.tick_params(labelsize=12.5)
ax3.set_xlabel("x (m)",fontsize=fts)
ax3.set_ylabel("y (m)",fontsize=fts)
ax3.text(0.01,0.92,"(c) Z=350 m",fontsize=fts, transform=ax3.transAxes)
#plt.savefig('inversion_model_Z_slices.jpg',dpi=300,bbox_inches='tight')
for ax in fig.axes:
    ax.tick_params(which='both',labelsize=fts)

plt.savefig('Inversion_Slices_interp.jpg',dpi=300,bbox_inches='tight')


