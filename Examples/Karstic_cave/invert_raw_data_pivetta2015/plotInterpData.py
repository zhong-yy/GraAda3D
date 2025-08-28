import sys
import os
#include_path=os.path.abspath('../../PythonPlot')
#sys.path.append(include_path)
#my module
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

fig=plt.figure(figsize=(10.5,6))
gs0 = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1])

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0],hspace=0.5,wspace=0.5)

min_value_shown=-1000
max_value_shown=1000
fts=12
clrmap='jet'

model='gz_result.txt'
#Y=700
Xi1,Zi1,Vi1=interpData.interp_yslice(model,y0=403840,
                         start_x=5062255,stop_x=5062575,num_x=100,
                         start_z=-315,stop_z=-100,num_z=100,neighbors=50,kernel='linear')
                         
ax1=fig.add_subplot(gs00[0])
pc=ax1.pcolormesh(Xi1,Zi1,Vi1,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax1.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_yaxis()
ax1.set_xlabel("x/Northing (m)",fontsize=fts)
ax1.set_ylabel("z/Depth (m)",fontsize=fts)
ax1.text(0.01,1.02,"(a) y=403840 m",fontsize=fts, transform=ax1.transAxes)

#Y=1450
Yi2,Zi2,Vi2=interpData.interp_xslice(model,x0=5062370,
                         start_y=403710,stop_y=403970,num_y=100,
                         start_z=-315,stop_z=-100,num_z=100,neighbors=50,kernel='linear')
ax2=fig.add_subplot(gs00[1])
pc=ax2.pcolormesh(Yi2,Zi2,Vi2,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax2.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.invert_yaxis()
ax2.set_xlabel("y/Easting (m)",fontsize=fts)
ax2.set_ylabel("z/Depth (m)",fontsize=fts)
ax2.text(0.01,1.02,"(b) x=5062360 m",fontsize=fts, transform=ax2.transAxes)


#Z=350
Xi3,Yi3,Vi3=interpData.interp_zslice(model,z0=-260,
                         start_x=5062255,stop_x=5062575,num_x=100,
                         start_y=403710,stop_y=403970,num_y=100,neighbors=50,kernel='linear')
gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],height_ratios=[1,7,1])
ax3=fig.add_subplot(gs01[1,0])
pc=ax3.pcolormesh(Yi3.T,Xi3.T,Vi3.T,cmap=clrmap,vmin=min_value_shown,vmax=max_value_shown,shading='gouraud')
ax3.set_aspect('equal', 'box')#相当于matlab中的axis equal
#ax3.set_title(r'Z=350m')
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.xaxis.set_minor_locator(AutoMinorLocator())
clb=fig.colorbar(pc,ax=ax3,shrink=1)
clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts)
clb.ax.tick_params(labelsize=12.5)
ax3.set_xlabel("y/Easting (m)",fontsize=fts)
ax3.set_ylabel("x/Northing (m)",fontsize=fts)
ax3.text(0.01,1.02,"(c) z=-260 m",fontsize=fts, transform=ax3.transAxes)
#plt.savefig('inversion_model_Z_slices.jpg',dpi=300,bbox_inches='tight')
for ax in fig.axes:
    ax.tick_params(which='both',labelsize=fts)
    ax.ticklabel_format(style="plain", scilimits=(0, 10), useOffset=False)
    
plt.savefig('Inversion_Slices_interp.jpg',dpi=300,bbox_inches='tight')


