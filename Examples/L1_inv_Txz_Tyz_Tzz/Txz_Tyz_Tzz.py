import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors


clrmap=mcolors.LinearSegmentedColormap.from_list("mycmap", ["magenta","blueviolet","royalblue","aqua","springgreen","lawngreen","yellow","orangered","red","white"])

x,y,gxz=np.loadtxt('./dobs_T_xz',unpack=True)
x,y,gyz=np.loadtxt('./dobs_T_yz',unpack=True)
x,y,gzz=np.loadtxt('./dobs_T_zz',unpack=True)


x,y,gxz_pre=np.loadtxt('./dpredicted_T_xz',unpack=True)
x,y,gyz_pre=np.loadtxt('./dpredicted_T_yz',unpack=True)
x,y,gzz_pre=np.loadtxt('./dpredicted_T_zz',unpack=True)

X=np.reshape(x,(41,41),order='C')
Y=np.reshape(y,(41,41),order='C')
GZZ=np.reshape(gzz,(41,41),order='C')
GXZ=np.reshape(gxz,(41,41),order='C')
GYZ=np.reshape(gyz,(41,41),order='C')

GZZ_predicted=np.reshape(gzz_pre,(41,41),order='C')
GXZ_predicted=np.reshape(gxz_pre,(41,41),order='C')
GYZ_predicted=np.reshape(gyz_pre,(41,41),order='C')


#fig,ax=plt.subplots(1,1,figsize=(6,4))

fig=plt.figure(figsize=(10,9))

plt.subplots_adjust(hspace=0.25,wspace=0.35)
#Txz
ax=plt.subplot(3,3,1)

pc=ax.pcolormesh(X,Y,GXZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'$T_{xz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,4)

pc=ax.pcolormesh(X,Y,GXZ_predicted,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Calculated $T_{xz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,7)

pc=ax.pcolormesh(X,Y,GXZ_predicted-GXZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Residuals of $T_{xz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

#Tyz
ax=plt.subplot(3,3,2)
pc=ax.pcolormesh(X,Y,GYZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'$T_{yz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,5)
pc=ax.pcolormesh(X,Y,GYZ_predicted,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Calculated $T_{yz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,8)
pc=ax.pcolormesh(X,Y,GYZ_predicted-GYZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Residuals of $T_{yz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")



#Tzz
ax=plt.subplot(3,3,3)

pc=ax.pcolormesh(X,Y,GZZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'$T_{zz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,6)

pc=ax.pcolormesh(X,Y,GZZ_predicted,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Calculated $T_{zz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

ax=plt.subplot(3,3,9)

pc=ax.pcolormesh(X,Y,GZZ_predicted-GZZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax.set_title(r'Residuals of $T_{zz}$')
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.invert_yaxis()
clb=fig.colorbar(pc,ax=ax,shrink=0.8)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'E')
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

plt.savefig('Txz_Tyz_Tzz.jpg',dpi=300,bbox_inches='tight')
