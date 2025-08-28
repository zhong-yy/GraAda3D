import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RBFInterpolator

x,y,dobs=np.loadtxt("./dobs_g_z",unpack=True)
x,y,dpre=np.loadtxt("./dpredicted_g_z",unpack=True)
residual=dobs-dpre


# nx=121
# ny=94

# X=np.reshape(x,(nx,ny),order='C')
# Y=np.reshape(y,(nx,ny),order='C')
# Dobs=np.reshape(dobs,(nx,ny),order='C')
# Dpre=np.reshape(dpre,(nx,ny),order='C')
# Residual=np.reshape(residual,(nx,ny),order='C')

cm=1/2.54
fig,axs=plt.subplots(2,2,figsize=(18*cm,16*cm),height_ratios=[1,1],width_ratios=[1,1])
plt.subplots_adjust(wspace=0.35,hspace=0.35)
fts=8

def plot_field(ax,xs,ys,values,subplot_title,fts=12,cbar_title='mGal',cmap='rainbow',invert_y=False):
    # c=ax.contourf(X,Y,FIELD,cmap=cmap)
    # clb=fig.colorbar(c,ax=ax,shrink=0.9,orientation="vertical")
    # clb_title=clb.ax.set_title(cbar_title,fontsize=fts)
    # clb.ax.tick_params(labelsize=fts)
    # ax.set_aspect('equal', 'box')
    # if invert_y:
    #     ax.invert_yaxis()
    # ax.tick_params(which='both',labelsize=fts)
    # ax.set_xlabel("y (m)",fontsize=fts)
    # ax.set_ylabel("x (m)",fontsize=fts)
    # ax.set_title(subplot_title,fontsize=fts)
# fig, ax = plt.subplots(1, 1)




    xy = np.hstack((xs[:, np.newaxis], ys[:, np.newaxis]))
    rbf = RBFInterpolator(xy, values, neighbors=50)

    min_x = np.min(xs)
    max_x = np.max(xs)
    min_y = np.min(ys)
    max_y = np.max(ys)

    Yi, Xi = np.meshgrid(np.linspace(min_y, max_y, 100), np.linspace(min_x, max_x, 100))
    xy_interp = np.stack((Xi.flatten(), Yi.flatten())).T
    data_i = rbf(xy_interp)
    data_i = np.reshape(data_i, Xi.shape)
    
    vmin = np.min(data_i)
    vmax = np.max(data_i)

    c = ax.contour(
        Yi,
        Xi,
        data_i,
        levels=20,
        linewidths=0.5,
        cmap=cmap,
    )
    c = ax.scatter(
        ys,
        xs,
        # s=20,
        c=values,
        cmap=cmap,
        edgecolors="black",
        marker="o",
        vmin=vmin,
        vmax=vmax,
    )
    clb = fig.colorbar(c, ax=ax, extend="both")
    clb_title=clb.ax.set_title(cbar_title,fontsize=fts)
    clb.ax.tick_params(labelsize=fts)    
    ax.ticklabel_format(style="plain", scilimits=(0, 10), useOffset=False)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("y/east (m)")
    ax.set_ylabel("x/north (m)")
    ax.set_title(subplot_title,fontsize=fts)
    if invert_y:
        ax.invert_yaxis()
    ax.tick_params(which='both',labelsize=fts)
    ax.set_xlabel("y (m)",fontsize=fts)
    ax.set_ylabel("x (m)",fontsize=fts)    


plot_field(axs[0][0],x,y,dobs,"Observed data",fts=fts)
plot_field(axs[0][1],x,y,dpre,"Predicted data",fts=fts)
plot_field(axs[1][0],x,y,residual,"Residuals",fts=fts)


axs[1][1].hist(residual,bins=20,color="lightgray",edgecolor="black")
axs[1][1].set_title("Histogram of residuals",fontsize=fts)
axs[1][1].set_xlabel(r"$d_{obs}-d_{pre} (mGal)$")
axs[1][1].tick_params(which='both',labelsize=fts)

plt.savefig('data_fitting.jpg',dpi=300,bbox_inches='tight')
