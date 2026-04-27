import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import scipy
from scipy.stats import norm
from netCDF4 import MFDataset
import xarray as xr


#This function calculates the dynamical heating according to equations (2) and (3) in Friedel et al., 2022a
#Needs w* and theta on pressure levels as input

def calculate_dyn_heating(nc_fid_t, nc_fid_w, nc_fid_theta, temp):
    
    #"nc_fid_t": file containing temperature variable
    #"nc_fid_w": file containing vertical velocity w*
    #"nc_fid_theta": file containing potential temperature
    #"temp": name of temperature variable
    
    # get vertical residual velocity
    w = nc_fid_w['vert_res_circ']
    plev_w = nc_fid_w['lev']
    lats = nc_fid_w['lat']
    time_w = nc_fid_w['time']

    w=w.sortby('time')
    w=w.resample(time="1MS").mean(dim="time") #calculate monthly means
    
    # get temperature
    temp = nc_fid_t[temp]
    temp=temp.sortby('time')
    temp=temp.resample(time="1MS").mean(dim="time") #calculate monthly means
    
    plev_t = nc_fid_t.variables['plev'][:]
    time = nc_fid_t.variables['time'][:]
    
    #get potential temperature
    theta=nc_fid_theta['THETA']
                
    theta=theta.sortby('time')
    theta=theta.resample(time="1MS").mean(dim="time") #calculate monthly means
    
    time_theta=theta.time

    theta=np.array(theta)
    
    # get pressure level in hPa if in Pa 
    if plev_w[len(plev_w)-1] > 1000:
        plev_w=plev_w/100
        
    # get pressure level in hPa if in Pa 
    if plev_t[len(plev_t)-1] > 1000:
        plev_t=plev_t/100    
    
    
    temp = temp.mean(dim='lon') # longitudinal mean
    temp = temp.mean(dim='plev') # calculate mean temperature in atmospheric column --> has to be weighted by plev thickness!!!!!!!!!
    
    temp = np.array(temp)
    
    temp_col=np.zeros((len(temp[:,0]), len(plev_w), len(temp[0,:])))

    # set temperature on each level to mean temperature
    for lev in range(len(plev_w)):
        temp_col[:,lev,:]=temp
    
    hlev = np.zeros((len(plev_w)))    # altitude in height instead of pressure
        
    for level in range(len(plev_w)):
        hlev[level] = -7000 * math.log(plev_w[level]/1013.25)
    
    theta=np.squeeze(theta)

    theta_grad = np.gradient(theta,hlev, axis=1)

    # calculate the dynamical heating rate in K/day
    
    w=np.reshape(w, (len(time_theta), len(plev_w), len(lats))) 
    
    dyn_hr = (temp_col) / theta * theta_grad  * w  * 24*60*60* (-1) 
    
    dyn_hr = xr.DataArray(dyn_hr, coords=[time,plev_w,lats], dims=['time','plev','lat'],name='dyn_hr')

    return dyn_hr
    