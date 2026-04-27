import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import scipy
from scipy.stats import norm
from scipy import signal
from find_FW import *
from eofs.standard import Eof
import xarray as xr 
import pandas as pd
from bootstrapping import *
from sklearn.linear_model import LinearRegression
import random 
import time

 
def interpolate_pole(time, lons, lats, array):
    
    #interpolate to 90°N and 360° longitude
    # extent latitude and longitude vector
    lats_90=np.zeros((len(lats)+2))
    
    for lat in range(1,len(lats)+1):
        lats_90[lat]=lats[lat-1]
        
    if lats[0]>0:    
        lats_90[0]=90    
        lats_90[len(lats)+1]=-90
        
    if lats[0]<0:      
        lats_90[0]=-90    
        lats_90[len(lats)+1]=90
        
    
    lons_360=np.zeros((len(lons)+1))
    
    for lon in range(len(lons)):
        lons_360[lon]=lons[lon]
        
    lons_360[len(lons)]=360
    
    
    #interpolate array to new grid points
    
    array_interpolated=np.zeros((len(lats)+2,len(lons)+1))
    
    array_interpolated[1:len(lats)+1,0:len(lons)]=array
    
    array_interpolated[1:len(lats)+1,len(lons)]=array[:,0]
  
    array_interpolated[0,:]=np.mean(array[0,:])
    array_interpolated[len(lats)+1,:]=np.mean(array[len(lats)-1,:])
        
    return lats_90, lons_360, array_interpolated  


#______________________________________________________________________________________


def interpolate_lons(time, lons, lats, array):
    
    lats_90 = lats
    
    # extent longitude vector 
    
    lons_360=np.zeros((len(lons)+1))
    
    for lon in range(len(lons)):
        lons_360[lon]=lons[lon]
        
    lons_360[len(lons)]=360
    
    #interpolate array to new grid point
    
    array_interpolated=np.zeros((len(lats),len(lons)+1))
    
    array_interpolated[:,0:len(lons)]=array
    
    array_interpolated[:,len(lons)]=array[:,0]
        
    
    return lats_90, lons_360, array_interpolated    
    
#______________________________________________________________________________________
  

def find_early_late_FW(FW_dates, years):
    
    #find years with the earliest and latest FW dates

    extreme_years=int(years/2) #divide years 50/50 (50% latest and 50% earliest FWs)
    
    FW_dates=np.array(FW_dates)
    
    FW_latest=FW_dates.argsort()[-extreme_years:][::-1]
    FW_earliest=FW_dates.argsort()[0:extreme_years]

    FW_latest_index=np.zeros((extreme_years))
    FW_earliest_index=np.zeros((extreme_years))

    
    for i in range(extreme_years):
        FW_earliest_index[i]=int(FW_dates[FW_earliest[i]])
        FW_latest_index[i]=int(FW_dates[FW_latest[i]])
        
        
    return FW_earliest, FW_latest, FW_earliest_index, FW_latest_index    

#______________________________________________________________________________________


def find_SSW_FW_dates(FW_dates, SSW_years):
    
    #this function selects the indices contained in the variable SSW_years
    #very badly programmed --> FW_SSW_index=FW_dates[SSW_dates] should give the same result
    
    FW_SSW_index=np.zeros((len(SSW_years)))
       
    for i in range(len(SSW_years)):
        
        if math.isnan(FW_dates[int(SSW_years[i])]) == True:
            
            FW_SSW_index[i]=np.nan
            
        else: 
             FW_SSW_index[i]=int(FW_dates[int(SSW_years[i])])
                
    return FW_SSW_index

#______________________________________________________________________________________



def find_ozone_extremes_FW(nc_fid, var, lev, years, extreme_years, model, FW_date, runmean):

    # This function finds years with the highest and lowest ozone values in March and April 
    # based on daily mean ozone values.
    
    #This function is slightly different from "find_ozone_extremes" in the module 
    #"ozone_extremes_leap", as it only consideres ozone values until the FW date 
    #to count towards high and low ozone years (in "find_ozone_extremes", all ozone values
    #in March and April are being considered, irrespective of the FW date. 
    
    # The module requires the input of a daily averaged zm ozone file on pressure levels.
    #"nc_fid": file containing ozone data
    #"var": name of ozone variable
    # "lev": name of pressure level variable
    # "years": number of years in input file 
    # "extreme_years": number of desired ozone extreme years.
    # "model": SOCOL, WACCM or MERRA
    # "FW_date": array of FW dates in each year
    #"runmean": if "True", a 5-day running mean is calculated for the ozone data before finding high/low ozone years

    
    partial_column=True # if True: calculates ozone extremes based on 30-70 hPa partial ozone column; 
                        # if False: calculates ozone extremes based on 70 hPa O3 mixing ratio

    O3=nc_fid[var]
    plev=nc_fid[lev]
    time=nc_fid['time']
    
    # interpolate to 2.5°
    
    O3=O3.interp(lat=np.linspace(-90,90,73))

    #calculate partial ozone column from 30 to 70 hPa over the polar cap

    delta_p = np.zeros((len(O3.time), len(plev)))

    m_air = 28.964/(6.022E23)
    g = 980.6

    if plev[len(plev)-1] <= 1000 and model!='MERRA': # for pressure levels in hPa

        for level in range(1,len(plev)):
            delta_p[:,level].fill( plev[level] - plev[level-1])

        #average over polar cap

        O3=O3.sel(lat=slice(60,90)) 

        weights = np.cos(np.deg2rad(O3.lat)) # latitudinal weights
        O3 = O3.weighted(weights)     
        O3=O3.mean(dim='lat')

        weights_p = xr.DataArray(delta_p*100, dims=['time','plev'], coords=[time,plev]) # difference between pressure levels in Pa

        O3 = O3 * weights_p * 10/ (g * m_air)

        O3=O3.sel(plev=slice(30,70)) 
        O3 = O3.sum(dim='plev')

        O3 = O3/2.687E16  #calculate DU

        plev=plev.sel(plev=slice(30,70))


    if plev[len(plev)-1] > 1000: # for pressure levels in hPa

        for level in range(1,len(plev)):
            delta_p[:,level].fill( plev[level] - plev[level-1])

        #average over polar cap    

        O3=O3.sel(lat=slice(60,90))

        weights = np.cos(np.deg2rad(O3.lat)) # latitudinal weights
        O3 = O3.weighted(weights)     

        O3=O3.mean(dim='lat')

        weights_p = xr.DataArray(delta_p, dims=['time','plev'], coords=[time,plev])

        O3 = O3 * weights_p * 10/ (g * m_air)

        O3=O3.sel(plev=slice(3000,7000))
        O3 = O3.sum(dim='plev')

        O3 = O3/2.687E16 #calculate DU

        plev=plev.sel(plev=slice(3000,7000))
        plev=plev/100    



    if model=='MERRA':

        for level in range(0,len(plev)-1):
            delta_p[:,level].fill( plev[level] - plev[level+1])

        #average over polar cap   

        O3=O3.sel(lat=slice(60,90))*28.970/47.9982

        weights = np.cos(np.deg2rad(O3.lat)) #latitudinal weights
        O3 = O3.weighted(weights)     
        O3=O3.mean(dim='lat')

        weights_p = xr.DataArray(delta_p*100, dims=['time','lev'], coords=[time,plev]) # difference between pressure levels in Pa

        O3 = O3 * weights_p * 10/ (g * m_air)

        O3=O3.sel(lev=slice(70,30)) 
        O3 = O3.sum(dim='lev')

        O3 = O3/2.687E16  #calculate DU

        plev=plev.sel(lev=slice(70,30))

    
    # select March and April
    
    var=O3
    
    if runmean==True:
        var=var.rolling(time=15, center=True).mean()
    
    O3_clim=O3.groupby("time.month").mean() #calculate monthly mean climatology
    
    #select values within March and April
    
    O3=O3.sel(time=nc_fid.time.dt.month.isin([3,4]))
    time=time.sel(time=nc_fid.time.dt.month.isin([3,4]))
    
    
        # select 41 random years to make sample comparable to MERRA2 (subsampling)
    
# =============================================================================
  #  O3_years_long=np.array(O3.time.dt.year)
  #  O3_years_long=set(O3_years_long.tolist())
#  
   # random_years = random.sample(O3_years_long, k=41)
#        
   # time=time.sel(time=O3.time.dt.year.isin(random_years))
   # O3=O3.sel(time=O3.time.dt.year.isin(random_years))
#     
   # years = 41
# =============================================================================
    
    
    # select highest and lowest ozone values in each year

    O3=O3.groupby("time.year")
    
    
    O3_highest_indices=np.empty((years))
    O3_highest_values=np.empty((years))
    
    O3_lowest_indices=np.empty((years))
    O3_lowest_values=np.empty((years))
    
    i=0
    
    for year, group in O3:
        
        if FW_date[i]-58 <= len(group) and FW_date[i]-58 > 0:
            
           # group_new=group[0:FW_date[i]-58+7]
        
            group_new=group[0:FW_date[i]]

            O3_highest_values[i]=int(group_new.max())
            O3_highest_indices[i]=int(group_new.argmax())

            O3_lowest_values[i]=int(group_new.min())
            O3_lowest_indices[i]=int(group_new.argmin())    
            
        else:

            O3_highest_values[i]=int(group.max())
            O3_highest_indices[i]=int(group.argmax())

            O3_lowest_values[i]=int(group.min())
            O3_lowest_indices[i]=int(group.argmin())
            
        i=i+1    
    
    time=time.groupby("time.year") 
    
    
    # find dates with the highest and lowest ozone values for each year
    
    O3_lowest_dates=[] # finds dates when ozone maximizes/minimizes each year
    O3_highest_dates=[]
    
    O3_years=np.zeros((years)) # list of all the years in the data

    for i,(year, group) in enumerate(time):

        O3_lowest_dates.append(np.array(group[int(O3_lowest_indices[i])]))
        O3_highest_dates.append(np.array(group[int(O3_highest_indices[i])]))
        O3_years[i]=year
    
    
    O3_highest_values=np.array(O3_highest_values)
    O3_lowest_values=np.array(O3_lowest_values)
    
    O3_highest_values=np.reshape(O3_highest_values, (years,))
    O3_lowest_values=np.reshape(O3_lowest_values, (years,))
    
    # find 50 highest and 50 lowest March and April ozone years based on daily values (this return the index of the respective year ozone extreme years)
    
    O3_highest=O3_highest_values.argsort()[-extreme_years:][::-1]
    O3_lowest=O3_lowest_values.argsort()[0:extreme_years]

    # find ozone value and date of the 50 extreme ozone years

    O3_lowest_index=np.zeros((extreme_years)) # this saves the index when the ozone extreme maximizes in each respective year
    O3_highest_index=np.zeros((extreme_years))
    
    O3_lowest_date=[] # this saves the date vector when the ozone extreme maximizes in each respective year
    O3_highest_date=[]

    
    for i in range(extreme_years):
        O3_lowest_index[i]=int(O3_lowest_indices[O3_lowest[i]])
        O3_highest_index[i]=int(O3_highest_indices[O3_highest[i]])
        
        O3_lowest_date.append(O3_lowest_dates[O3_lowest[i]])
        O3_highest_date.append(O3_highest_dates[O3_highest[i]])
        
    
    O3_lowest_index=O3_lowest_index.astype(int)
    O3_highest_index=O3_highest_index.astype(int)
    
    
    print('Mean minimum ozone day: ' + str(np.mean(O3_lowest_index)) + ' ± ' + str(np.std(O3_lowest_index)))
    print('Mean maximum ozone day: ' + str(np.mean(O3_highest_index)) + ' ± ' + str(np.std(O3_highest_index)))
    
    
# =============================================================================
   # for i in range(extreme_years): # comment this for FW calculations (this saves the actual year instead of the year index)
   #     O3_lowest[i]=O3_years[O3_lowest[i]]
   #     O3_highest[i]=O3_years[O3_highest[i]]
# =============================================================================


    O3_intersect=len(np.intersect1d(O3_highest, O3_lowest)) #this counts the number of years where there is both an ozone maximum and minimum


    
    print('Numer of intersections: ' + str(O3_intersect))
  

    return O3_highest, O3_lowest, np.reshape(O3_lowest_date, (extreme_years,)), np.reshape(O3_highest_date, (extreme_years,)), np.reshape(O3_lowest_index, (extreme_years,)), np.reshape(O3_highest_index, (extreme_years,)), O3_intersect, O3_years
  
    
#______________________________________________________________________________________




def FW_surf_pattern(nc_fid, var, years, FW_dates, FW_dates_index, model, O3_years):
    
    #calculate the surface pattern following the FW date in specific years (e.g. high/low ozone years or early/late FW years)
    
    #"nc_fid": file containing surface variable
    #"var": name of surface variable
    #"years": number of years in the record
    #"FW_dates": years in which the surface response of FWs should be calcuated (e.g. high/low ozone years)
    #"FW_dates_index": timing (day of the year) of the FW in the specific years
    #"model": "SOCOL", "WACCM", or "MERRA"
    #"O3_years": all years in the record
    
    var=nc_fid[var]
    
    lons=nc_fid['lon']
    lons=np.array(lons)
    lats=nc_fid['lat']
    lats=np.array(lats)
    time=nc_fid['time']

    extreme_years=int(years/2)

    # Calculate anomalies
  
    if model=='MERRA': #For MERRA2, exclude year 2020 from calculation of anomalies
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time")
        years=40
    else:
        var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")
    

    var_anomalies_xr=var_anomalies
    var_anomalies_np=np.array(var_anomalies)
    var_anomalies=var_anomalies.sel(time=nc_fid.time.dt.month.isin([1,2,3,4,5,6,7,8])) #select months Jan-Aug
    time=time.sel(time=nc_fid.time.dt.month.isin([1,2,3,4,5,6,7,8]))


    var_anomalies=var_anomalies.groupby("time.year")

    var_FW=np.empty((len(FW_dates),30,len(lats),len(lons)))

    i=0
    j=0
    
    for year, group in var_anomalies: #loop through all years in the record
        
        if i in FW_dates: #if the year is contained in the years of interest (e.g. high/low ozone years) 

            index=np.where(FW_dates==i)  #Finds location of the specific early/late FW within the array
        
            var_FW[j,:,:,:]=group[int(FW_dates_index[index]):int(FW_dates_index[index]+30),:,:] #select the 30 days after the FW
        
            j=j+1
        
        i=i+1    
        
        
    var_FW = np.nanmean(var_FW[:,:,:,:], axis=1)
    
    var_FW = np.nanmean(var_FW, axis=0)
    

    #calculate significance according to bootstrapping
    if model=='MERRA' or model=='SOCOL':
        significance = bootstrapping_leap_FW(nc_fid, var_anomalies_xr, FW_dates_index, var_FW, O3_years)

    if model=='WACCM':
        significance = bootstrapping_FW(np.reshape(np.array(var_anomalies_xr), (200,365, len(lats), len(lons))), FW_dates_index, var_FW,200)



    if 90 and -90 not in lats:
        lats_90, lons_360, var_FW = interpolate_pole(0, lons, lats, var_FW)
        lats_90, lons_360, significance = interpolate_pole(0, lons, lats, significance)

    else:
        lats_90, lons_360, var_FW  = interpolate_lons(0, lons, lats, var_FW)
        lats_90, lons_360, significance  = interpolate_lons(0, lons, lats, significance)

    #new coordinates for plotting
    m = Basemap(projection='ortho', lat_0=90,lon_0=0,resolution='l')
    lons_m,lats_m=np.meshgrid(lons_360,lats_90)
    xpt,ypt = m(lons_m,lats_m)
    
    return xpt, ypt, var_FW, significance
    


#______________________________________________________________________________________
    
def calculate_AO_FW(nc_fid, var, lev,  FW_dates, FW_dates_index, model ):

    
    # This function calculates the Arctic oscillation index in the -30/+60 days time period around the ozone extreme event.
    
    #"nc_fid": file containing geopotential height variable
    #"var": name of geopotential height variable
    #"lev": name of pressure level variable
    #"FW_dates": years of interest (e.g. high/low ozone years)
    #"FW_dates_index": indices of FWs in years of interest
    #"model": "WACCM", "SOCOL", or "MERRA"
    
    
    manual=False
    
    var=nc_fid[var]
    plev=nc_fid[lev]
    lats=nc_fid['lat']
    
    #get orientation of latitude dimension
    if lats[0]>0:
        lats_orientation = 'negative'
    else:
        lats_orientation = 'positive'

        
  #  plev=np.array(plev)/100
    
    if plev[len(plev)-1] > 1000:
        var=var.interp(plev=100000)
        
    if plev[len(plev)-1] <= 1000 and model!='MERRA':   
        var=var.interp(plev=1000)
        
    if model=='MERRA':
        var=var.interp(lev=1000)  
       
    
    var=var.interp(lat=np.linspace(-90,90,73))
    lats=lats.interp(lat=np.linspace(-90,90,73))
    
    lats=lats.sel(lat=slice(20,90))
    lats=np.array(lats)

    #calculate anomalies
    if model=='MERRA':
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time")
    else:
          var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")
    
 
    var_anomalies=var_anomalies.sel(lat=slice(20,90)) 


    time=nc_fid['time']
    time=np.array(time)
    
    coslat = np.cos(np.deg2rad(lats).clip(0., 1.))
    wgts = np.sqrt(coslat)[..., np.newaxis]
    
    # get geopotential height of polar cap 

    gh_layer = np.array(var_anomalies)
 
    if len(np.shape(gh_layer)) == 3:
        gh_layer = np.reshape(gh_layer, (np.shape(gh_layer)[0], len(lats)))

    #calculate EOF loading pattern
    solver = Eof(gh_layer, weights=wgts, center=True)
    
    #calculate PC time series
    AO=np.reshape(solver.pcs(npcs=1, pcscaling=1), (np.shape(gh_layer)[0],))
    
    AO = xr.DataArray(AO, coords=[time], dims=['time'],name='AO')
    AO_xr = xr.DataArray(AO, coords=[time], dims=['time'])  #save as xarray
    
    AO = AO.groupby('time.year')
    
    AO_FW=np.zeros((len(FW_dates),60))
    

    i=0
    j=0
    
    
    for year, group in AO: #loop through all years in AO variable
        
        if i in FW_dates: #if year is contained in years of interest:
    
            index=int(np.array(np.where(FW_dates==i)))  #get FW index in year of interest
        
            AO_FW[j,:]=group[int(FW_dates_index[index]-30):int(FW_dates_index[index]+30)]
        
            j=j+1
        
        i=i+1    
                    
    nof_lats = len(lats)
  
    max_AO = AO.argmax()
    
    #check if sign of AO is correct
    if lats_orientation=='negative':
               
             if np.nanmean(gh_layer[max_AO,0:int(nof_lats/2)-1])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+1:nof_lats-1])<0:
                     AO_FW[:]=-AO_FW[:]
               
    else:
              
               if np.nanmean(gh_layer[max_AO,0:int(nof_lats/2)-3])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+3:nof_lats-1])<0:
                       AO_FW[:]=-AO_FW[:]
            
    print('Mean AO index difference 30 days after max. ozone depletion for low ozone years:  '   + str(np.mean(AO_FW[:,30:59])) +' ± ' + str(np.std(-AO_FW[:,30:59])))  

    
    return AO_FW


      
#______________________________________________________________________________________
    
    
    
def FW_stratosphere_seas(nc_fid, var, lev, years,  FW_dates, model, anomalies, lats, late, runmean):
    
    # This function creates a time - altitude plot of the seasonal anomalies in specific years ("FW_dates") 
    # of a certain variable.
    # It needs input of a zonal mean daily mean variable on pressure levels
    
    #"nc_fid": file containing the variable
    #"var": name of the variable of interest
    #"lev": name of pressure level variable
    #"years": number of years in the record
    #"FW_dates": years of interest (e.g high/low ozone years)
    #"model": WACCM, SOCOL, or MERRA
    #"anomalies": if "True": calculates anomalies, if "False": calculates absolute values
    #"lats/late": start and end latitude for latitudinal averaging
    #"runmean": if "True": calculates 30 day running means (for noisy variables like EPF divergence)

  
    factor=1
    var=nc_fid[var]*factor
 
    time=nc_fid['time']
    plev=nc_fid[lev]

    
    # get pressure level in hPa if in Pa and get pressure levels between 1 and 1000 hPa

    var=var.interp(lat=np.linspace(-90,90,73))
    
    if plev[len(plev)-1] <= 1000 and model!='MERRA':
              var=var.sel(lat=slice(lats,late), plev=slice(1,1000)) 
            #  var=var.where(var < 30)
              plev=plev.sel(plev=slice(1,1000))    

    if plev[len(plev)-1] > 1000:
              var=var.sel(lat=slice(lats,late), plev=slice(100,100000))
              plev=plev.sel(plev=slice(100,100000))
              plev=np.array(plev)/100
           #   var=var.where(var < 30)

    if model=='MERRA':
              var=var.sel(lat=slice(lats,late), lev=slice(1000,1))
              years=40
              plev=plev.sel(lev=slice(1000,1))
            
    
    weights = np.cos(np.deg2rad(var.lat))
    weights.name = "weights"
    var = var.weighted(weights)    
    
    var=var.mean(dim='lat')
    
    
    if anomalies is False:
        var_anomalies=var
                            
    if anomalies is True:
        
        var_anomalies=var.groupby('time.dayofyear')-var.groupby('time.dayofyear').mean()
     
    if runmean is True: 

        var_anomalies=var_anomalies.rolling(time=30, center=True).mean()
    
    
    var_anomalies=var_anomalies.sel(time=nc_fid.time.dt.month.isin([3,4,5])) #select springtime

    var_anomalies_out=var_anomalies
    
    time=time.sel(time=nc_fid.time.dt.month.isin([3,4,5]))
    
    time_days=var_anomalies.sel(plev=100).groupby('time.dayofyear').mean()

    var_anomalies=var_anomalies.groupby("time.year")

    if model == 'SOCOL': var_FW=np.empty((len(FW_dates),len(time_days)-1,len(plev)))
    else: var_FW=np.empty((len(FW_dates),len(time_days),len(plev)))
    
    i=0
    j=0
    
    for year, group in var_anomalies:
        
        if i in FW_dates:
            if model== 'SOCOL': var_FW[j,:,:]=np.squeeze(np.array(group[0:len(time_days)-1,:]))
            else: var_FW[j,:,:]=np.squeeze(np.array(group[0:len(time_days),:]))
            j=j+1
        
        i=i+1    
    
     
    var_FW=np.where(var_FW ==0, np.nan, var_FW)
    if anomalies is True: var_FW=np.where(np.abs(var_FW)>30, 0, var_FW)
    var_FW = np.nanmean(var_FW, axis=0)

    if anomalies==True:
            levels = np.linspace(-20,20, 21)
        
    if anomalies==False:    
            levels = np.linspace(-15,15, 21)
    
    return var_FW, var_anomalies_out


#______________________________________________________________________________________




def FW_vertical_ozone_years(file, nc_fid, lev, O3_lowest, O3_highest, extreme_years, Data):

    #reads in all FW dates across altitudes for all years and selects the FW profile
    #for high/low ozone years
    
    #"file": file containing all FW dates across altitudes
    #"nc_fid": file containing pressure variable
    #"lev": name of pressure variable
    #"O3_lowest/O3_highest": high/low ozone years
    #"extreme_years": number of extreme years considered (number of high/low ozone years)
    #"Data": WACCM, SOCOL, or MERRA

    FW_dates_lev = np.load(file)
    plev=nc_fid[lev]
    
    if Data=='MERRA':
        plev=plev.sel(lev=slice(1,0.1))
    else:    
        if plev[len(plev)-1] > 1000:
            plev=plev.sel(plev=slice(10,100))
            plev=plev/100
        else:
           plev=plev.sel(plev=slice(0.1,1))

    #save anomalies of FW dates in years of interest
    high_FW_vertical = np.empty((len(FW_dates_lev[:,0]), extreme_years))
    low_FW_vertical = np.empty((len(FW_dates_lev[:,0]), extreme_years))
    
    #save absolute FW dates
    high_FW_vertical_full = np.empty((len(FW_dates_lev[:,0]), extreme_years))
    low_FW_vertical_full = np.empty((len(FW_dates_lev[:,0]), extreme_years))
    
    for year in range(extreme_years):
        low_FW_vertical_full[:,year] = FW_dates_lev[:,O3_lowest[year]]
        low_FW_vertical[:,year] = FW_dates_lev[:,O3_lowest[year]] - np.mean(FW_dates_lev, axis=1)
        
    for year in range(extreme_years):
        
        high_FW_vertical[:,year] = FW_dates_lev[:,O3_highest[year]] - np.mean(FW_dates_lev, axis=1)
        high_FW_vertical_full[:,year] = FW_dates_lev[:,O3_highest[year]]
        

    print("mean FW date: " + str(np.mean(FW_dates_lev, axis=1)))
    print("mean FW date in low ozone years: " + str(np.mean(low_FW_vertical_full, axis=1)))
    print("mean FW date in high ozone years: " + str(np.mean(high_FW_vertical_full, axis=1)))        
        
    return low_FW_vertical, high_FW_vertical


#______________________________________________________________________________________



"""
BELOW, THE FOLLOWING PART CONTAINS THE BOOTSTRAPPING FUNCTIONS USED FOR FRIEDEL ET AL. (2022).
NOTE THAT THEY ARE VERY SIMILAR TO THE BOOTSTRAPPING FUNCTIONS CONTAINED IN THE MODULE bootstrapping.py.
PLEASE REFER TO THIS MODULE FOR DETAILED DESCRIPTION OF THE FUNCTIONS.
"""


def bootstrapping_alt_time_seas(array1, array2, diff, O3_years1, O3_years2, nc_fid_1, nc_fid_2, extreme_years, model):
    
    # input arrays must be reduced to March-May

    time_days=array1.sel(plev=100).groupby('time.dayofyear').mean()
    
    bootstrap=500
    
    array1=np.squeeze(array1)
    array2=np.squeeze(array2)
  
    print('starting bootstrapping....')
    
    if model =='SOCOL': bootstrap_composite=np.zeros((bootstrap,len(time_days)-1,len(array1[0,:])))
    else: bootstrap_composite=np.zeros((bootstrap,len(time_days),len(array1[0,:])))

    i=0
    while i<bootstrap:

        # for each event, find date in random year and calculate composite of -60 to + 60 days after each SSW and append to array var_all_SSWs
 
        if model=='SOCOL':
            var_all_1= np.zeros((extreme_years,len(time_days)-1,len(array1[0,:])))
            var_all_2= np.zeros((extreme_years,len(time_days)-1,len(array2[0,:])))  
        if model=='WACCM':
            var_all_1= np.zeros((extreme_years,len(time_days),len(array1[0,:])))
            var_all_2= np.zeros((extreme_years,len(time_days),len(array2[0,:])))  
        
        for index in range(extreme_years):

            #randomly select one year of the first sample

            random_number_1 = random.randint(0,len(O3_years1)-1)
            random_year_1 = int(O3_years1[random_number_1])

            array_random_year_1 = array1.sel(time=array1.time.dt.year.isin([random_year_1]))
            
            if model=='SOCOL': var_all_1[index,:,:]=np.squeeze(array_random_year_1[0:len(time_days)-1,:])
            else: var_all_1[index,:,:]=np.squeeze(array_random_year_1[0:len(time_days),:])


            #randomly select one year of the second sample

            random_number_2 = random.randint(0,len(O3_years2)-1)
            random_year_2 = int(O3_years2[random_number_2])
            array_random_year_2 = array2.sel(time=array2.time.dt.year.isin([random_year_2]))  
            
            
            if model=='SOCOL': var_all_2[index,:,:]=np.squeeze(array_random_year_2[0:len(time_days)-1,:])
            else:  var_all_2[index,:,:]=np.squeeze(array_random_year_2[0:len(time_days),:])

                 # calculate composite over all SSWs for one bootstrap

        
        bootstrap_composite[i,:,:] = np.mean(var_all_1, axis=0) - np.mean(var_all_2, axis=0)
        i=i+1
        

    # calculate mean and significance of 500 sample composite
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - diff)

    # calculate significane: anomaly significant, if it differs more than 2 standard deviations from the composite mean value
    significance=np.greater(np.abs(diff_bootstrap), 2*std_bootstrap)
    
    print('boostrapping done...')

    return significance


#______________________________________________________________________________________


def bootstrapping_alt_time_seas_np(array1, array2, diff, O3_years1, O3_years2, nc_fid_1, nc_fid_2, extreme_years, model):
    
    # input arrays must be reduced to March-May

    
    time_days=len(array1[0,:,0])
    
    bootstrap=500
    
    array1=np.squeeze(array1)
    array2=np.squeeze(array2)
  
    print('starting bootstrapping....')
    
    if model =='SOCOL': bootstrap_composite=np.zeros((bootstrap,time_days,len(array1[0,0,:])))
    else: bootstrap_composite=np.zeros((bootstrap,time_days,len(array1[0,0,:])))

    i=0
    while i<bootstrap:

        # for each event, find date in random year and calculate composite of -60 to + 60 days after each SSW and append to array var_all_SSWs
 
        if model=='SOCOL':
            var_all_1= np.zeros((extreme_years,time_days,len(array1[0,0,:])))
            var_all_2= np.zeros((extreme_years,time_days,len(array2[0,0,:])))  
        if model=='WACCM':
            var_all_1= np.zeros((extreme_years,time_days,len(array1[0,0,:])))
            var_all_2= np.zeros((extreme_years,time_days,len(array2[0,0,:])))  
        
        for index in range(extreme_years):

            #randomly select one year of the first sample

            random_number_1 = random.randint(0,len(O3_years1)-1)

            array_random_year_1 = array1[random_number_1,:,:]
            
            if model=='SOCOL': var_all_1[index,:,:]=np.squeeze(array_random_year_1[0:time_days,:])
            else: var_all_1[index,:,:]=np.squeeze(array_random_year_1[0:time_days,:])


            #randomly select one year of the second sample

            random_number_2 = random.randint(0,len(O3_years2)-1)
            array_random_year_2 = array2[random_number_2,:]
            
            
            if model=='SOCOL': var_all_2[index,:,:]=np.squeeze(array_random_year_2[0:time_days,:])
            else:  var_all_2[index,:,:]=np.squeeze(array_random_year_2[0:time_days,:])

                 # calculate composite over all SSWs for one bootstrap
        
        bootstrap_composite[i,:,:] = np.mean(var_all_1, axis=0) - np.mean(var_all_2, axis=0)
        i=i+1
        

    # calculate mean and significance of 500 sample composite
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - diff)

    # calculate significane: anomaly significant, if it differs more than 2 standard deviations from the composite mean value
    significance=np.greater(np.abs(diff_bootstrap), 2*std_bootstrap)
    
    print('boostrapping done...')

    return significance



#______________________________________________________________________________________


def bootstrapping_leap_FW(nc_fid, array, indices, composite, O3_years):

    bootstrap=500

    # perform bootstrapping 5000 times
    bootstrap_composite=np.zeros((bootstrap,len(array[0,:,0]),len(array[0,0,:])))
    
    
    print("starting bootstrapping...")
    i=0
    while i<bootstrap:
        
        # for each event, find date in random year and calculate composite of 30 days after each SSW and append to array var_all_SSWs
        #print(i)
        var_all= np.zeros((len(indices)*30,len(array[0,:,0]),len(array[0,0,:])))
        
        
        for index in range(len(indices)):
            
            indice = int(indices[index]) # as index counting starts in March only
         

            random_number = random.randint(0,len(O3_years)-1)
            random_year = O3_years[random_number]
            
            array_random_year = array.sel(time=nc_fid.time.dt.year.isin([random_year]))
            
            var_all[index*30:index*30+30,:,:]=array_random_year[indice:indice+30, :, :]
  
        # calculate composite over all SSWs for one bootstrap
        
        bootstrap_composite[i,:,:] = np.mean(var_all, axis=0)
    
        i=i+1
        
        
    print("bootstrapping finished.")
    # calculate mean and significance of 500 sample composite
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    
    
    # calculate significane: anomaly significant, if it differs more than 2 standard deviations from the composite mean value
    significance=np.greater(np.abs(diff_bootstrap), 2*std_bootstrap)
    
    return significance



#______________________________________________________________________________________




def bootstrapping_FW(array, indices, composite, years):

    bootstrap=500
    
    # perform bootstrapping 5000 times
    bootstrap_composite=np.zeros((bootstrap,len(array[0,0,:,0]),len(array[0,0,0,:])))
    
    print("starting bootstrapping...")
    
    i=0
    while i<bootstrap:
        
        # for each SSW, find date in random year and calculate composite of 30 days after each SSW and append to array var_all_SSWs
      #  print(i)
        var_all= np.zeros((len(indices)*30,len(array[0,0,:,0]),len(array[0,0,0,:])))
        
        for index in range(len(indices)):
            indice = int(indices[index]) +59
            random_year = random.randint(0, years-1)
            var_all[index*30:index*30+30,:,:]=array[random_year, indice:indice+30, :, :]
  
        # calculate composite over all SSWs for one bootstrap
        
        
        bootstrap_composite[i,:,:] = np.mean(var_all, axis=0)
    
        i=i+1
        
        
    print("bootstrapping completed.")
    # calculate mean and significance of 500 sample composite
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    
    
    # calculate significane: anomaly significant, if it differs more than 2 standard deviations from the composite mean value
    significance=np.greater(np.abs(diff_bootstrap), 2*std_bootstrap)
    
    return significance


#______________________________________________________________________________________


