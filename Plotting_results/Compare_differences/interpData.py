import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import RBFInterpolator

def interp_xslice(filename,x0,start_y,stop_y,num_y,start_z,stop_z,num_z,neighbors=None,kernel='thin_plate_spline'):
    #Load data
    x,y,z,v=np.loadtxt(filename,skiprows=2,usecols=[6,7,8,9],unpack=True)
    xyz=np.hstack((x[:,np.newaxis],y[:,np.newaxis],z[:,np.newaxis]))
    
    #Construct interpolation grid points
    xi=np.array([x0])
    yi=np.linspace(start_y,stop_y,num_y)
    zi=np.linspace(start_z,stop_z,num_z)
    XI,YI,ZI=np.meshgrid(xi,yi,zi)
    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()
    
    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)
     
    YI_2d=YI[:,0,:].T # shape of YI is (ny,nx,nz)
    ZI_2d=ZI[:,0,:].T # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,0,:].T
    return YI_2d,ZI_2d,VI_2d
    
def interp_yslice(filename,y0,start_x,stop_x,num_x,start_z,stop_z,num_z,neighbors=None,kernel='thin_plate_spline'):
    #Load data
    x,y,z,v=np.loadtxt(filename,skiprows=2,usecols=[6,7,8,9],unpack=True)
    xyz=np.hstack((x[:,np.newaxis],y[:,np.newaxis],z[:,np.newaxis]))
    
    #Construct interpolation grid points
    yi=np.array([y0])
    xi=np.linspace(start_x,stop_x,num_x)
    zi=np.linspace(start_z,stop_z,num_z)
    XI,YI,ZI=np.meshgrid(xi,yi,zi)
    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()
    
    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)
     
    XI_2d=XI[0,:,:].T # shape of XI is (ny,nx,nz)
    ZI_2d=ZI[0,:,:].T # shape of ZI is (ny,nx,nz)
    VI_2d=VI[0,:,:].T
    return XI_2d,ZI_2d,VI_2d
    
def interp_zslice(filename,z0,start_x,stop_x,num_x,start_y,stop_y,num_y,neighbors=None,kernel='thin_plate_spline'):
    #Load data
    x,y,z,v=np.loadtxt(filename,skiprows=2,usecols=[6,7,8,9],unpack=True)
    xyz=np.hstack((x[:,np.newaxis],y[:,np.newaxis],z[:,np.newaxis]))
    
    #Construct interpolation grid points
    zi=np.array([z0])
    xi=np.linspace(start_x,stop_x,num_x)
    yi=np.linspace(start_y,stop_y,num_y)
    XI,YI,ZI=np.meshgrid(xi,yi,zi)
    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()
    
    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)
     
    XI_2d=XI[:,:,0] # shape of XI is (ny,nx,nz)
    YI_2d=YI[:,:,0] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,:,0]
    return XI_2d,YI_2d,VI_2d
    

