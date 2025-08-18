import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
#from xyz2XYZ import xyz2XYZ
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors


clrmap=mcolors.LinearSegmentedColormap.from_list("mycmap", ["magenta","blueviolet","royalblue","aqua","springgreen","lawngreen","yellow","orangered","red","white"])
clrmap="rainbow"
for fname,outfile_remark in zip(['result/data_no_noise','result/data'],["","_with_noise"]):
    x,y,z,gz,Tzz,Txz,Tyz,Txx,Txy,Tyy=np.loadtxt(fname,unpack=True)
    # _,_,gz2=np.loadtxt('../../../Examples0/dobs_no_noise_g_z',unpack=True)
    # _,_,Tzz2=np.loadtxt('../../../Examples0/dobs_no_noise_T_zz',unpack=True)
    # _,_,Txz2=np.loadtxt('../../../Examples0/dobs_no_noise_T_xz',unpack=True)
    # _,_,Tyz2=np.loadtxt('../../../Examples0/dobs_no_noise_T_yz',unpack=True)
    # gz=gz-gz2
    # Tzz=Tzz-Tzz2
    # Txz=Txz-Txz2
    # Tyz=Tyz-Tyz2
    X=np.reshape(x,(41,41),order='C')
    Y=np.reshape(y,(41,41),order='C')
    GZ=np.reshape(gz,(41,41),order='C')

    TXX=np.reshape(Txx,(41,41),order='C')
    TXY=np.reshape(Txy,(41,41),order='C')
    TXZ=np.reshape(Txz,(41,41),order='C')
    TYY=np.reshape(Tyy,(41,41),order='C')
    TYZ=np.reshape(Tyz,(41,41),order='C')
    TZZ=np.reshape(Tzz,(41,41),order='C')


    fig,ax1=plt.subplots(1,1,figsize=(6,4))    

    pc=ax1.pcolormesh(X,Y,GZ,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
    ax1.contour(X,Y,GZ,levels=np.arange(-1,1,0.1),colors='k')
    ax1.set_aspect('equal', 'box')#相当于matlab中的axis equal
    ax1.set_title(r'$g_z$')
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
    plt.savefig(f'gz{outfile_remark}.jpg',dpi=300,bbox_inches='tight')

    

    ## Tzz
    # fig=plt.figure(figsize=(15,4),layout="constrained")
    fig,axs=plt.subplots(3,3,figsize=(8,6),layout="constrained")
    # plt.subplots_adjust(wspace=0.25,hspace=0.25)
    # fig,axs=plt.subplot_mosaic([["Txx","Txy","Txz"],[".","Tyy","Tyz"],[".",".","Tzz"]],figsize=(15,4),layout="constrained")

    for i,(name,field) in enumerate(zip([r'$T_{xx}$',r'$T_{xy}$',r'$T_{xz}$',None,r'$T_{yy}$',r'$T_{yz}$',None,None,r'$T_{zz}$'],[TXX,TXY,TXZ,None,TYY,TYZ,None,None,TZZ])):
        # ax=plt.subplot(3,3,i+1)
        ax=axs.flat[i]
        if field is not None:
            # ax=axs[i][2]
            pc=ax.pcolormesh(X,Y,field,shading='gouraud',cmap=clrmap)#vmin,vmax 设置色标的上下限
            ax.contour(X,Y,field,levels=np.arange(-30,30,5),colors='k')
            ax.set_aspect('equal', 'box')#相当于matlab中的axis equal
            
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.invert_yaxis()
            clb=fig.colorbar(pc,ax=ax,shrink=1)
            clb.ax.yaxis.set_minor_locator(AutoMinorLocator())
            clb_title=clb.ax.set_title(r'E',loc="left",pad=3)
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            ax.set_title(name)
        else:
            fig.delaxes(ax)
            # ax.set_visible(False)
    # axs[0][0].set_title(r'$T_{xx}$')
    # axs[0][1].set_title(r'$T_{xy}$')
    # axs[0][2].set_title(r'$T_{xz}$')

    # axs[0][2].set_title(r'$T_{xz}$')
    # axs[1][2].set_title(r'$T_{yz}$')
    # axs[2][2].set_title(r'$T_{zz}$')
    plt.savefig(f'T{outfile_remark}.jpg',dpi=300,bbox_inches='tight')
