
#plot interpolated data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec


#own module
import interpData
import similarity

def normalized_to_1_1(x):
    result=(x-np.min(x))/(np.max(x)-np.min(x))
    return result
    
x,y,z,density=np.loadtxt('gz_result.txt',skiprows=2,usecols=[6,7,8,9],unpack=True)
print(np.min(density),np.max(density))



#gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0],hspace=0.4)

min_value_shown1=-300
max_value_shown1=300
min_value_shown2=-300
max_value_shown2=300
min_value_shown3=-1
max_value_shown3=1
fts=12
clrmap1='jet'
clrmap='Blues_r'

model1='gz_result.txt'
model2='gz_result_with_ref.txt'
comparison_method=similarity.normalized_cross_correlation
win_size=(30,30)

#fig=plt.figure(figsize=(12,6))
#gs0 = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1])

fig, axs=plt.subplots(3,3,figsize=(12,8.))
plt.subplots_adjust(hspace=0.25,wspace=0.35)
#Y=700
Xi1,Zi1,slice1_model1=interpData.interp_yslice(model1,y0=700,
                         start_x=0,stop_x=2000,num_x=201,
                         start_z=0,stop_z=1000,num_z=101,neighbors=10,kernel='linear')
Xi1,Zi1,slice1_model2=interpData.interp_yslice(model2,y0=700,
                         start_x=0,stop_x=2000,num_x=201,
                         start_z=0,stop_z=1000,num_z=101,neighbors=10,kernel='linear')
norm_slice1_model1=normalized_to_1_1(slice1_model1)
norm_slice1_model2=normalized_to_1_1(slice1_model2)
pc=axs[0,0].pcolormesh(Xi1,Zi1,slice1_model1,cmap=clrmap1,vmin=min_value_shown1,vmax=max_value_shown1,shading='gouraud')
axs[0,0].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].invert_yaxis()
axs[0,0].set_xlabel("x (m)",fontsize=fts)
axs[0,0].set_ylabel("z (m)",fontsize=fts)
axs[0,0].set_title("Result 1")
axs[0,0].text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=axs[0,0].transAxes)
print(np.min(slice1_model2),np.max(slice1_model1))    
pc=axs[0,1].pcolormesh(Xi1,Zi1,slice1_model2,cmap=clrmap1,vmin=min_value_shown2,vmax=max_value_shown2,shading='gouraud')
axs[0,1].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].invert_yaxis()
axs[0,1].set_xlabel("x (m)",fontsize=fts)
axs[0,1].set_title("Result 2")
#axs[0,1].set_ylabel("z (m)",fontsize=fts)
#axs[0,1].text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=axs[0,1].transAxes)

     

Vi1=similarity.compute_sim(slice1_model1,slice1_model2,win_size,method=comparison_method)      
print(np.min(Vi1),np.max(Vi1))      
pc=axs[0,2].pcolormesh(Xi1,Zi1,Vi1,cmap=clrmap,vmin=min_value_shown3,vmax=max_value_shown3,shading='gouraud')
axs[0,2].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[0,2].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,2].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,2].invert_yaxis()
axs[0,2].set_xlabel("y (m)",fontsize=fts)
axs[0,2].set_title("Normalized cross correlation")
#axs[0,2].set_ylabel("z (m)",fontsize=fts)
#axs[0,2].text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=axs[0,2].transAxes)

#Y=1450
Xi2,Zi2,slice2_model1=interpData.interp_yslice(model1,y0=1450,
                         start_x=0,stop_x=2000,num_x=201,
                         start_z=0,stop_z=1000,num_z=101,neighbors=10,kernel='linear')
Xi2,Zi2,slice2_model2=interpData.interp_yslice(model2,y0=1450,
                         start_x=0,stop_x=2000,num_x=201,
                         start_z=0,stop_z=1000,num_z=101,neighbors=10,kernel='linear')
norm_slice2_model1=normalized_to_1_1(slice2_model1)
norm_slice2_model2=normalized_to_1_1(slice2_model2)
pc=axs[1,0].pcolormesh(Xi2,Zi2,slice2_model1,cmap=clrmap1,vmin=min_value_shown1,vmax=max_value_shown1,shading='gouraud')
axs[1,0].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].invert_yaxis()
axs[1,0].set_xlabel("x (m)",fontsize=fts)
axs[1,0].set_ylabel("z (m)",fontsize=fts)
axs[1,0].text(0.01,0.9,"(b) Y=1450 m",fontsize=fts, transform=axs[1,0].transAxes)

pc=axs[1,1].pcolormesh(Xi2,Zi2,slice2_model2,cmap=clrmap1,vmin=min_value_shown2,vmax=max_value_shown2,shading='gouraud')
axs[1,1].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].invert_yaxis()
axs[1,1].set_xlabel("x (m)",fontsize=fts)
#axs[1,1].set_ylabel("z (m)",fontsize=fts)
#axs[1,1].text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=ax1.transAxes)


Vi2=similarity.compute_sim(slice2_model1,slice2_model2,win_size,method=comparison_method)
#ax2=fig.add_subplot(gs00[1])
print(np.min(Vi2),np.max(Vi2))
pc=axs[1,2].pcolormesh(Xi2,Zi2,Vi2,cmap=clrmap,vmin=min_value_shown3,vmax=max_value_shown3,shading='gouraud')
axs[1,2].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[1,2].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,2].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,2].invert_yaxis()
axs[1,2].set_xlabel("x (m)",fontsize=fts)
#axs[1,2].set_ylabel("z (m)",fontsize=fts)
#axs[1,2].text(0.01,0.9,"(b) Y=1450 m",fontsize=fts, transform=ax2.transAxes)


#Z=350
Xi3,Yi3,slice3_model1=interpData.interp_zslice(model1,z0=350,
                         start_x=0,stop_x=2000,num_x=201,
                         start_y=0,stop_y=2000,num_y=201,neighbors=10,kernel='linear')
Xi3,Yi3,slice3_model2=interpData.interp_zslice(model2,z0=350,
                         start_x=0,stop_x=2000,num_x=201,
                         start_y=0,stop_y=2000,num_y=201,neighbors=10,kernel='linear')
norm_slice3_model1=normalized_to_1_1(slice3_model1)
norm_slice3_model2=normalized_to_1_1(slice3_model2)
pc=axs[2,0].pcolormesh(Xi3,Yi3,slice3_model1,cmap=clrmap1,vmin=min_value_shown1,vmax=max_value_shown1,shading='gouraud')
axs[2,0].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[2,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[2,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[2,0].invert_yaxis()
axs[2,0].set_xlabel("x (m)",fontsize=fts)
axs[2,0].set_ylabel("y (m)",fontsize=fts)
axs[2,0].text(0.01,0.92,"(c) Z=350 m",fontsize=fts, transform=axs[2,0].transAxes)
clb=fig.colorbar(pc,ax=axs[2,0],shrink=0.85)
clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts,loc='left')
clb.ax.tick_params(labelsize=12.5)

pc=axs[2,1].pcolormesh(Xi3,Yi3,slice3_model2,cmap=clrmap1,vmin=min_value_shown2,vmax=max_value_shown2,shading='gouraud')
axs[2,1].set_aspect('equal', 'box')#相当于matlab中的axis equal
axs[2,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[2,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[2,1].invert_yaxis()
axs[2,1].set_xlabel("x (m)",fontsize=fts)
#axs[2,1].set_ylabel("z (m)",fontsize=fts)
#axs[2,1].text(0.01,0.9,"(a) Y=700 m",fontsize=fts, transform=ax1.transAxes)
clb=fig.colorbar(pc,ax=axs[2,1],shrink=0.85)
clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts,loc='left')
clb.ax.tick_params(labelsize=12.5)



Vi3=similarity.compute_sim(slice3_model1,slice3_model2,win_size,method=comparison_method)
print(np.min(Vi3),np.max(Vi3))  
#gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],height_ratios=[1,7,1])
#ax3=fig.add_subplot(gs01[1,0])
pc=axs[2,2].pcolormesh(Xi3,Yi3,Vi3,cmap=clrmap,vmin=min_value_shown3,vmax=max_value_shown3,shading='gouraud')
axs[2,2].set_aspect('equal', 'box')#相当于matlab中的axis equal
#ax3.set_title(r'Z=350m')
axs[2,2].yaxis.set_minor_locator(AutoMinorLocator())
axs[2,2].xaxis.set_minor_locator(AutoMinorLocator())
axs[2,2].invert_yaxis()
clb=fig.colorbar(pc,ax=axs[2,2],shrink=0.85)
#clb_title=clb.ax.set_title(r'kg/m$^3$',fontsize=fts)
clb.ax.tick_params(labelsize=12.5)
axs[2,2].set_xlabel("x (m)",fontsize=fts)
#axs[2,2].set_ylabel("y (m)",fontsize=fts)
#axs[2,2].text(0.01,0.92,"(c) Z=350 m",fontsize=fts, transform=ax3.transAxes)
#plt.savefig('inversion_model_Z_slices.jpg',dpi=300,bbox_inches='tight')
for ax in fig.axes:
    ax.tick_params(which='both',labelsize=fts)

plt.savefig('sim.jpg',dpi=300,bbox_inches='tight')


