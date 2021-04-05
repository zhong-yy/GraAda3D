import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
#from xyz2XYZ import xyz2XYZ
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors


clrmap=mcolors.LinearSegmentedColormap.from_list("mycmap", ["magenta","blueviolet","royalblue","aqua","springgreen","lawngreen","yellow","orangered","red","white"])
#clrmap="jet"

x,y,gz=np.loadtxt('dobs_g_z',unpack=True)
X=np.reshape(x,(41,41),order='C')
Y=np.reshape(y,(41,41),order='C')
GZ=np.reshape(gz,(41,41),order='C')


fig,ax1=plt.subplots(1,1,figsize=(6,4))    

pc=ax1.pcolormesh(X,Y,GZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
ax1.set_aspect('equal', 'box')#相当于matlab中的axis equal
ax1.set_title(r'$g_z$')
#ax1.yaxis.set_major_locator(MultipleLocator(1))
#ax1.xaxis.set_major_locator(MultipleLocator(1))
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_yaxis()
clb=fig.colorbar(pc,ax=ax1,shrink=1)
clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
#clb.ax.ticklabel_format(style='scientific',scilimits=(0,0),useMathText=True)
clb_title=clb.ax.set_title(r'mGal',loc="left",pad=3)
#clb_title.set_position((1.2,1))
#clb.ax.yaxis.get_offset_text().set_position((1,1))
ax1.set_xlabel("x (m)")
ax1.set_ylabel("y (m)")
plt.savefig('gz.jpg',dpi=300,bbox_inches='tight')
