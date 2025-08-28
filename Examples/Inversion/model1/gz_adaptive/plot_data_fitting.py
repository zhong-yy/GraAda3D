import matplotlib.pyplot as plt
import numpy as np

x,y,z,dobs=np.loadtxt("./dobs_g_z",unpack=True)
x,y,z,dpre=np.loadtxt("./dpredicted_g_z",unpack=True)
residual=dobs-dpre


nx=41
ny=41

X=np.reshape(x,(ny,nx),order='C')
Y=np.reshape(y,(ny,nx),order='C')
Dobs=np.reshape(dobs,(ny,nx),order='C')
Dpre=np.reshape(dpre,(ny,nx),order='C')
Residual=np.reshape(residual,(ny,nx),order='C')

cm=1/2.54
fig,axs=plt.subplots(2,2,figsize=(16*cm,15*cm),height_ratios=[1,1],width_ratios=[1,1])
plt.subplots_adjust(wspace=0.35,hspace=0.35)
fts=8

def plot_field(ax,X,Y,FIELD,subplot_title,fts=12,cbar_title='mGal',cmap='rainbow',invert_y=False):
    c=ax.contourf(X,Y,FIELD,cmap=cmap,levels=20)
    clb=fig.colorbar(c,ax=ax,shrink=0.9,orientation="vertical")
    clb_title=clb.ax.set_title(cbar_title,fontsize=fts)
    clb.ax.tick_params(labelsize=fts)
    ax.set_aspect('equal', 'box')
    if invert_y:
        ax.invert_yaxis()
    ax.tick_params(which='both',labelsize=fts)
    ax.set_xlabel("x (m)",fontsize=fts)
    ax.set_ylabel("y (m)",fontsize=fts)
    ax.set_title(subplot_title,fontsize=fts)


plot_field(axs[0][0],X,Y,Dobs,"Observed data",fts=fts,invert_y=True)
plot_field(axs[0][1],X,Y,Dpre,"Predicted data",fts=fts,invert_y=True)
plot_field(axs[1][0],X,Y,Residual,"Residuals",fts=fts,invert_y=True)


axs[1][1].hist(residual,bins=20,color="lightgray",edgecolor="black")
axs[1][1].set_title("Histogram of residuals",fontsize=fts)
axs[1][1].set_xlabel(r"$d_{obs}-d_{pre} (mGal)$")
axs[1][1].tick_params(which='both',labelsize=fts)

plt.savefig('data_fitting.jpg',dpi=300,bbox_inches='tight')
