import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import scipy
from scipy.stats import norm
from scipy import signal
from find_FW import *   # scripts for the calculation of the Final Stratospheric Warming date
from eofs.standard import Eof
import xarray as xr 
import pandas as pd
from bootstrapping import * # scripts for the bootstrapping calculation
from sklearn.linear_model import LinearRegression
import random
import time



def interpolate_pole(time, lons, lats, array):
    
    #This function interpolates an array to 90°N and 90°S latitude and 360° longitude
    
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
  

#_______________________________________________________________________________






def find_ozone_extremes_mm(nc_fid, var, lev, years, extreme_years, model, nc_fid_SURF):

    # This function has been used to produce supplementary Fig. S9
    
    # This function finds the 25% of strongest ozone minima in March (monthly mean) and calculates
    # SLP anomaly composites in March and April of the same years 
    
    # "years": number of years in input file 
    # "extreme_years": number of desired ozone extreme years.
    # "model": SOCOL, WACCM or MERRA
    # "var": name of ozone variable in input file
    # "lev": name of pressure level variable
    
    
    partial_column=True # if True: calculates ozone extremes based on 30-70 hPa partial ozone column; 
                        # if False: calculates ozone extremes based on 70 hPa O3 mixing ratio
        
        
    O3=nc_fid[var] # O3 variable
    plev=nc_fid[lev] # xarray containing pressure levels
    time=nc_fid['time'] # time xarray
    

    # interpolate to 2.5°
    
    O3=O3.interp(lat=np.linspace(-90,90,73)) #interpolate to a 2.5° latitudinal grid
    
    if partial_column ==True:
    
        #calculate partial ozone column from 30 to 70 hPa over the polar cap
    
        delta_p = np.zeros((len(O3.time), len(plev)))
        m_air = 28.964/(6.022E23)
        g = 980.6
            
        if plev[len(plev)-1] <= 1000 and model!='MERRA': # for pressure levels in hPa
            
            for level in range(1,len(plev)):
                delta_p[:,level].fill( plev[level] - plev[level-1])

            O3=O3.sel(lat=slice(60,90)) #select polar cap
    
            weights = np.cos(np.deg2rad(O3.lat)) # latitudinal weights 
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')
            
            weights_p = xr.DataArray(delta_p*100, dims=['time','plev'], coords=[time,plev]) # difference between pressure levels in Pa

            O3 = O3 * weights_p * 10/ (g * m_air)
            
            O3=O3.sel(plev=slice(30,70)) 
            O3 = O3.sum(dim='plev')
           
            O3 = O3/2.687E16  # convert in DU

            plev=plev.sel(plev=slice(30,70))
            
            
        if plev[len(plev)-1] > 1000: # for pressure levels in Pa
            
            for level in range(1,len(plev)):
                delta_p[:,level].fill( plev[level] - plev[level-1])
                
            O3=O3.sel(lat=slice(60,90)) # average over polar cap    

            weights = np.cos(np.deg2rad(O3.lat)) # latitudinal weights 
            O3 = O3.weighted(weights)     
   
            O3=O3.mean(dim='lat')
  
            weights_p = xr.DataArray(delta_p, dims=['time','plev'], coords=[time,plev])
            
            O3 = O3 * weights_p * 10/ (g * m_air)
            
            O3=O3.sel(plev=slice(3000,7000))
            O3 = O3.sum(dim='plev')
           
            O3 = O3/2.687E16 # convert in DU

            plev=plev.sel(plev=slice(3000,7000))
            plev=plev/100    # convert to hPa

    
        if model=='MERRA':
            
            for level in range(0,len(plev)-1):
                delta_p[:,level].fill( plev[level] - plev[level+1])

            O3=O3.sel(lat=slice(60,90))*28.970/47.9982 # convert in mol/mol

            weights = np.cos(np.deg2rad(O3.lat)) # latitudinal weights
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')

            weights_p = xr.DataArray(delta_p*100, dims=['time','lev'], coords=[time,plev]) # difference between pressure levels in Pa
            
            O3 = O3 * weights_p * 10/ (g * m_air)
 
            O3=O3.sel(lev=slice(70,30)) 
            O3 = O3.sum(dim='lev')
           
            O3 = O3/2.687E16  # convert in DU

            plev=plev.sel(lev=slice(70,30))

        
    O3=O3.sel(time=nc_fid.time.dt.month.isin([3])) # select March
    O3=O3.groupby("time.year").mean() # calculate March monthly means
    
   
    # find 50 highest and 50 lowest March and April ozone years based on March monthly mean values 
    
    O3=np.reshape(np.array(O3), (years,))
    
    O3_highest=O3.argsort()[-extreme_years:][::-1] #50 highest ozone values
    O3_lowest=O3.argsort()[0:extreme_years] #50 lowest ozone values
    
    
    # calculate SLP anomalies in March and April of ozone minimum years
    
    SLP=nc_fid_SURF['SLP'] # SLP array
    time=nc_fid_SURF['time']
    
    lats=nc_fid_SURF['lat']
    lons=nc_fid_SURF['lon'] 
    
    SLP = SLP.groupby("time.dayofyear") - SLP.groupby("time.dayofyear").mean("time") # calculate daily anomalies
    
    SLP=SLP.sel(time=nc_fid_SURF.time.dt.month.isin([3,4])) #select March and April 
    time=time.sel(time=nc_fid_SURF.time.dt.month.isin([3,4]))
    SLP=SLP.groupby("time.year").mean() # calculate mean over March and April of each year
    
    SLP=np.array(SLP)
    
    SLP_low=np.zeros((extreme_years,len(lats),len(lons)))
    
    #select ozone minimum years
    
    for year in range(extreme_years):
        SLP_low[year,:,:]=SLP[O3_lowest[year],:,:]
        
    t_array, p_array = scipy.stats.ttest_1samp(SLP_low,0, axis=0)   #calculate significance based on t-test
    SLP_low=np.mean(SLP_low,axis=0)    # calculate composite mean
    
    
    # interpolate to latitude of -90 and 90° and longitude of 360°
     
    if 90 and -90 not in lats:
        lats_90, lons_360, SLP_low  = interpolate_pole(0, lons, lats, SLP_low)
        lats_90, lons_360, p_array  = interpolate_pole(0, lons, lats, p_array)
        
    else:
        lats_90, lons_360, SLP_low = interpolate_lons(0, lons, lats, SLP_low)
        lats_90, lons_360, p_array  = interpolate_lons(0, lons, lats, p_array)
    
    
    return SLP_low
    
   
#_______________________________________________________________________________



def find_ozone_extremes(nc_fid, var, lev, years, extreme_years, model):
    
    
    # This function finds years with the highest and lowest ozone values in March and April 
    # based on daily mean ozone values.
    
    # The module requires the input of a daily averaged zm ozone file on pressure levels.
    
    # "nc_fid": input file including ozone variable
    # "years": number of years in input file 
    # "extreme_years": number of desired ozone extreme years.
    # "model": SOCOL, WACCM or MERRA
    # "var": name of ozone variable in input file
    # "lev": name of pressure level variable
    
    partial_column=True # if True: calculates ozone extremes based on 30-70 hPa partial ozone column; 
                        # if False: calculates ozone extremes based on 70 hPa O3 mixing ratio

    O3=nc_fid[var]
    plev=nc_fid[lev]
    time=nc_fid['time']
    
    O3=O3.interp(lat=np.linspace(-90,90,73))     # interpolate to 2.5°
    
    
    if partial_column ==True:
        
        #calculate partial ozone column from 30 to 70 hPa over the polar cap
    
        delta_p = np.zeros((len(O3.time), len(plev)))
    
        m_air = 28.964/(6.022E23)
        g = 980.6
            
        if plev[len(plev)-1] <= 1000 and model!='MERRA': # for pressure levels in hPa
            
            for level in range(1,len(plev)):
                delta_p[:,level].fill( plev[level] - plev[level-1])

            O3=O3.sel(lat=slice(60,90)) #average over polar cap
            
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
                
            O3=O3.sel(lat=slice(60,90))  #average over polar cap  
            
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
            
            O3=O3.sel(lat=slice(60,90))*28.970/47.9982 #average over polar cap and conversion to mol/mol
     
            weights = np.cos(np.deg2rad(O3.lat)) #latitudinal weights
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')
            
            weights_p = xr.DataArray(delta_p*100, dims=['time','lev'], coords=[time,plev]) # difference between pressure levels in Pa
            
            O3 = O3 * weights_p * 10/ (g * m_air)
            
            O3=O3.sel(lev=slice(70,30)) 
            O3 = O3.sum(dim='lev')
           
            O3 = O3/2.687E16  #calculate DU

            plev=plev.sel(lev=slice(70,30))

        
    if partial_column==False:

        # average over polar cap. Choose values at 70hPa
        if plev[len(plev)-1] > 1000:
            O3=O3.sel(lat=slice(60,90),plev=7000)
            
        if plev[len(plev)-1] <= 1000 and model!='MERRA':
            O3=O3.sel(lat=slice(60,90),plev=70) 
            
        if model=='MERRA':
            O3=O3.sel(lat=slice(60,90),lev=70)
     
        weights = np.cos(np.deg2rad(O3.lat))
        weights.name = "weights"
        
        O3 = O3.weighted(weights)     
        
        O3=O3.mean(dim='lat')*1000000 # calculate ppmv
        
    
    # select March and April
    
   # O3=O3.groupby("time.dayofyear")-O3.groupby("time.dayofyear").mean() #detection based on anomalies
    var=O3
 
    
    O3_clim=O3.groupby("time.month").mean() #calculate monthly mean climatology
    
    #select values within March and April
    
    O3=O3.sel(time=nc_fid.time.dt.month.isin([3,4]))
    time=time.sel(time=nc_fid.time.dt.month.isin([3,4]))


    # select 41 random years to make sample comparable to MERRA2 (subsampling)
    
# =============================================================================
#     O3_years_long=np.array(O3.time.dt.year)
#     O3_years_long=set(O3_years_long.tolist())
#  
#     random_years = random.sample(O3_years_long, k=41)
#        
#     time=time.sel(time=O3.time.dt.year.isin(random_years))
#     O3=O3.sel(time=O3.time.dt.year.isin(random_years))
#     
#     years = 41
# =============================================================================

    
    # select highest and lowest ozone values in each year
    #因为版本不兼容所以把armax改了
    O3_highest_values=O3.groupby("time.year").max() #select maximum ozone value in each spring
   # O3_highest_indices=O3.groupby("time.year").argmax() #select index of maximum ozone value in each spring
   # O3_highest_indices = O3.groupby("time.year").idxmax()
    O3_highest_indices = O3.groupby("time.year").map(lambda x: x.argmax(dim='time'))
    
    O3_lowest_values=O3.groupby("time.year").min() #select minimum ozone value in each spring
   # O3_lowest_indices=O3.groupby("time.year").argmin() #select index of minimum ozone value in each spring
   # O3_lowest_indices = O3.groupby("time.year").idxmin()
    O3_lowest_indices = O3.groupby("time.year").map(lambda x: x.argmin(dim='time'))
    time=time.groupby("time.year") 
    
    # find dates with the highest and lowest ozone values for each year
    
    O3_lowest_dates=[] # finds dates when ozone maximizes/minimizes each year
    O3_highest_dates=[]
    
    O3_years=np.zeros((years)) # list of all the years in the data

    for i,(year, group) in enumerate(time): #loop through every year and save date when lowest/highest ozone values occur

        O3_lowest_dates.append(np.array(group[O3_lowest_indices[i]]))
        O3_highest_dates.append(np.array(group[O3_highest_indices[i]]))
        O3_years[i]=year
    
    
    O3_highest_values=np.array(O3_highest_values)
    O3_lowest_values=np.array(O3_lowest_values)
    
    O3_highest_values=np.reshape(O3_highest_values, (years,))
    O3_lowest_values=np.reshape(O3_lowest_values, (years,))
    
    # find 50 highest and 50 lowest March and April ozone years based on daily values (this return the index of the respective year ozone extreme years)
    
    O3_highest=O3_highest_values.argsort()[-extreme_years:][::-1] #find 50 highest ozone values
    O3_lowest=O3_lowest_values.argsort()[0:extreme_years] #find 50 lowest ozone values

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
    

    O3_intersect=len(np.intersect1d(O3_highest, O3_lowest)) #this counts the number of years where there is both an ozone maximum and minimum

    #This function returns:
    
    #"O3_highest": indices of the maximum ozone years
    #"O3_lowest": indices of the minimum ozone years
    #np.reshape(O3_lowest_date, (extreme_years,)): dates of occurrence of the 50 strongest ozone minima
    #np.reshape(O3_highest_date, (extreme_years,)): dates of occurrence of the 50 strongest ozone maxima
    #np.reshape(O3_lowest_index, (extreme_years,)): index (days after March 1st) when strongest 50 ozone minima occur
    #np.reshape(O3_highest_index, (extreme_years,)): index (days after March 1st) when strongest 50 ozone maxima occur
    #"O3_intersect": number of years when both an ozone minimum and maximum occurs
    #"O3_years": year labels as in the netCDF input file

    return O3_highest, O3_lowest, np.reshape(O3_lowest_date, (extreme_years,)), np.reshape(O3_highest_date, (extreme_years,)), np.reshape(O3_lowest_index, (extreme_years,)), np.reshape(O3_highest_index, (extreme_years,)), O3_intersect, O3_years
  
    

#______________________________________________________________________________________



def FW_diff_ozone_extremes(nc_fid, var, lev, years, extreme_years, O3_lowest, O3_highest, O3_lowest_index, model):
    
    
    # This function calculates the final warming date of low and high ozone years and 
    # plots it in a histogram and as density functions.
    
    # This function needs input in form of daily zonal mean wind on pressure levels
    # as well as the number of years in the data file, the number of extreme ozone years
    # and two arrays containing a list with the highest and lowest ozone years. 
    
    # "nc_fid": input file containing the wind variable
    # "years": number of years in input file 
    # "extreme_years": number of desired ozone extreme years.
    # "model": SOCOL, WACCM or MERRA
    # "var": name of wind variable in input file
    # "lev": name of pressure level variable
    
    
    plev=nc_fid[lev]
    plev=np.array(plev)
    
    FW_dates=find_FW_new_leap(nc_fid, years, var, lev, model, 10) # this function calculates the FW date in each year at 10 hPa

    FW_dates_high_o3=np.zeros(extreme_years) #FSW dates in high ozone years
    FW_dates_low_o3=np.zeros(extreme_years)  #FSW dates in low ozone years
    
    for i in range(extreme_years):   
        FW_dates_high_o3[i]=FW_dates[O3_highest[i]]
        FW_dates_low_o3[i]=FW_dates[O3_lowest[i]]
        
    print('Mean final warming date: ' + str(np.mean(FW_dates)) + '  low ozone years: ' + str(np.mean(FW_dates_low_o3)) + '   high ozone years: ' +str(np.mean(FW_dates_high_o3)))
    
    date_diff=FW_dates_low_o3-O3_lowest_index
    
    print('Difference between mean low ozone and mean final warming date: ' + str(np.mean(date_diff)))
    
    plt.figure(figsize=(10,8)) 
    x=np.linspace(np.mean(FW_dates)-50,np.mean(FW_dates)+50,100)
    plt.plot(x, norm.pdf(x, np.mean(FW_dates), np.std(FW_dates)), color='k' , label='all years')
    plt.plot(x, norm.pdf(x, np.mean(FW_dates_low_o3), np.std(FW_dates_low_o3)), color='blue' , label='low ozone years')
    plt.plot(x, norm.pdf(x, np.mean(FW_dates_high_o3), np.std(FW_dates_high_o3)), color='red' , label='high ozone years')
    plt.xlabel('day of the year', fontsize=18)
    plt.ylabel('density', fontsize=18)
    plt.xticks([59,90,120,151], ('March', 'Apr', 'May','June'), fontsize=18) 
    
    plt.figure(figsize=(10,8)) 
    plt.hist(FW_dates, bins=20, color='grey', density=True)
    plt.hist(FW_dates_low_o3, bins=20, color='blue', alpha=0.4, density=True)
    plt.hist(FW_dates_high_o3, bins=20, color='red', alpha=0.4, density=True)
    plt.xlabel('day of the year', fontsize=18)
    plt.ylabel('density', fontsize=18)
    plt.xticks([59,90,120,151], ('March', 'Apr', 'May','June'), fontsize=18) 
    plt.axvline(x=np.mean(FW_dates), ymin=0, ymax=17, linestyle='--', color='grey')
    plt.axvline(x=np.mean(FW_dates_high_o3), ymin=0, ymax=17, linestyle='--', color='red')
    plt.axvline(x=np.mean(FW_dates_low_o3), ymin=0, ymax=17, linestyle='--', color='blue')    

    plt.figure(figsize=(10,8)) 
    plt.hist(date_diff, bins=20, color='yellow', density=True)
    plt.xlabel('day of the year', fontsize=18)
    plt.ylabel('density', fontsize=18)
             
    return  date_diff #differences between FSW date and ozone minimum date



#_________________________________________________________________________________________   
    
    

    
    
    
def analyse_precip(nc_fid1, nc_fid2, var1, var2, extreme_years, O3_highest_date, O3_lowest_date, O3_lowest_index, O3_years,model):

    # This function calculates and plots the total precipitation response for high and low 
    # ozone years.
    
    # It requires input in form of daily surface fiels (lat - lon) of the variable "var1" and "var2" (convective and large-scale precipitation).
    
    # "extreme_years": number of desired ozone extreme years.
    # if "surface_pressure=True": the input variable is surface pressure in Pa. The variable will be divided by 100 to get result in hPa
    # "O3_highest_date": date of the ozone maxima
    # "O3_highest_date": date of the ozone minima
    # "O3_lowest_index": index (day after March 1st) when ozone minima occur
    # "O3_years": years in which ozone minima occur
    # "model": SOCOL, WACCM or MERRA
    
    var1=nc_fid1[var1]
    var2=nc_fid1[var2]
    
    lons=nc_fid1['lon']
    lons=np.array(lons)
    lats=nc_fid1['lat']
    lats=np.array(lats)

    # Calculate anomalies
    
    var=var1+var2 #calculates total precipitation

    var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")

    # FIND INDICES OF OZONE EXTREMES   
    
    var_anomalies_xr=var_anomalies
    time=nc_fid1['time']

    time=np.array(time)
    
    ozone_low_date=np.zeros((extreme_years))
    ozone_high_date=np.zeros((extreme_years))
    
    for year in range(extreme_years):

         ozone_low_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_lowest_date[year]))), (1,))) #select variable on the date when the ozone minimum occurs
         ozone_high_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_highest_date[year]))), (1,)))

    # get data following the 30 days after an ozone extreme event
    
    var_min=np.zeros((extreme_years,30,len(lats),len(lons)))

    for i in range(extreme_years):   
  
            var_min[i,:,:,:]=var_anomalies[int(ozone_low_date[i]):int(ozone_low_date[i]+30),:,:]
           
    # average over the 30-day period
    
    var_min = np.nanmean(var_min[:,:,:,:], axis=1)
    var_min_zm = np.nanmean(var_min, axis=0)

    # calculate significance of difference of high and low ozone years

    if model=='MERRA' or model=='SOCOL':
        significance_low = bootstrapping_leap(nc_fid1, var_anomalies, O3_lowest_index, var_min_zm, O3_years )
    if model=='WACCM':
        significance_low = bootstrapping(np.reshape(np.array(var_anomalies), (200,365, len(lats), len(lons))), O3_lowest_index, var_min_zm, 200)

    # interpolate the values (to get values at the poles)
    
    if 90 and -90 not in lats:
        lats_90, lons_360, var_min_zm  = interpolate_pole(0, lons, lats, var_min_zm)
     #   lats_90, lons_360, p_array_min = interpolate_pole(0, lons, lats, p_array_min)
        lats_90, lons_360, significance_low = interpolate_pole(0, lons, lats, significance_low)
       
    else:
        lats_90, lons_360, var_min_zm  = interpolate_lons(0, lons, lats, var_min_zm)
    #    lats_90, lons_360, p_array_min = interpolate_lons(0, lons, lats, p_array_min)
        lats_90, lons_360, significance_low = interpolate_lons(0, lons, lats, significance_low)
      
    #function returns:
    
    #"var_min_zm": mean surface response in the 30 days following the ozone minima (lat-lon)
    # "lons_360" & "lats_90": coordinates used for plotting
    # "var_anomalies_xr": xarray containing the surface anomalies
    # "significance_low": significance of the surface response in low ozone years (lat-lon)
    

    return var_min_zm, lons_360, lats_90, var_anomalies_xr, significance_low
    

#_________________________________________________________________________________________________________________
    
   
    
    
def plot_time_alt(nc_fid, var, lev, extreme_years, O3_highest_date, O3_lowest_date, model, O3_lowest_index, O3_years):
    
    # This function creates a time - altitude plot of low/high ozone years
    # of a certain variable around the extreme ozone date.
    
    # It needs input of a zonal mean daily mean variable on pressure levels
    # nc_fid: file containing variable
    # "var": name of desired variable
    # "extreme_years": number of desired ozone extreme years.
    # "O3_highest_date": dates of the ozone maxima
    # "O3_lowest_date": dates of the ozone minima
    # "model": SOCOL, WACCM or MERRA
    # "O3_lowest_index": index (day of the year after March 1st) when ozone minima occur
    # "O3_years": years in which ozone minima occur
  
    factor=1 # the variable will be multiplied by this factor
    
    var=nc_fid[var]*factor
    plev=nc_fid[lev]

    var.fillna(0)
    var=var.interp(lat=np.linspace(-90,90,73))

    # average over polar cap. Choose values between 1 and 1000 hPa
      
    if plev[len(plev)-1] <= 1000 and model!='MERRA':
          var=var.sel(lat=slice(60,90), plev=slice(1,1000)) 
          plev=plev.sel(plev=slice(1,1000))
    
    if plev[len(plev)-1] > 1000:
          var=var.sel(lat=slice(60,90), plev=slice(100,100000))
          plev=plev.sel(plev=slice(100,100000))
          plev=np.array(plev)/100
          
    if model=='MERRA':
          var=var.sel(lat=slice(60,90), lev=slice(1000,1))
          plev=plev.sel(lev=slice(1000,1))
               
    weights = np.cos(np.deg2rad(var.lat)) # latitudinal weights
    weights.name = "weights"
    var = var.weighted(weights)    
    
    var=var.mean(dim='lat') # average over the polar cap
    
    var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time") # calculates anomalies
    

    # FIND INDICES OF OZONE EXTREMES   
    time=nc_fid['time']
    time=np.array(time)

    ozone_low_index=np.zeros((extreme_years))
    ozone_high_index=np.zeros((extreme_years))

    #select time index when ozone minima occur
    
    for year in range(extreme_years):
         ozone_low_index[year]=int(np.array(np.where(time==np.array(O3_lowest_date[year]))))
         ozone_high_index[year]=int(np.array(np.where(time==np.array(O3_highest_date[year]))))
    
   # get time frame -30 to +60 days after the ozone extremes 
         
    var_min=np.zeros((extreme_years,90,len(plev)))
    var_max=np.zeros((extreme_years,90,len(plev)))
    
    var_anomalies_xr=var_anomalies
    var_anomalies=np.array(var_anomalies) 
    var_anomalies[var_anomalies==0]=0
    
    for i in range(extreme_years):   
        var_min[i,:,:]=np.reshape(var_anomalies[int(ozone_low_index[i]-29):int(ozone_low_index[i]+61),:], (90,len(plev)))
        var_max[i,:,:]=np.reshape(var_anomalies[int(ozone_high_index[i]-29):int(ozone_high_index[i]+61),:], (90,len(plev)))
    
    # average over all years
    
    var_min_all=var_min 
    var_high_all=var_max
    
    var_min = np.mean(var_min, axis=0)
    var_max = np.mean(var_max, axis=0)

    #calculate significance of anomalies based on a bootstrapping test
    if model=='SOCOL' or model=='MERRA':
        significance=bootstrapping_leap_alt_time(nc_fid, var_anomalies_xr, O3_lowest_index, var_min, O3_years, model)
    if model=='WACCM': 
        significance=bootstrapping_leap_alt_time(nc_fid, np.reshape(np.array(var_anomalies_xr), (200,365,len(plev))), O3_lowest_index, var_min, O3_years, model)

    # This function returns:
    
    #"var_high_all": contains the 90-day time frame of anomalies around the ozone maxima date for ALL years
    #"var_min_all": contains the 90-day time frame of anomalies around the ozone minima date for ALL years  
    #"var_min": contains the averages 90-day time fram of anomalies around the ozone minima
    #"var_anomalies_xr": contains the variable anomalies 
    #"plev": contains the pressure levels between 1 and 1000 hPa
    #"significance": contains the significance vector for the mean 90-day time frame
    
    return var_high_all, var_min_all, var_min, var_anomalies_xr, plev, significance
    
  
#______________________________________________________________________________________



def calculate_AO(nc_fid, var, lev, extreme_years, O3_highest_date, O3_lowest_date, model, O3_lowest_index, O3_years):
    
    # This function calculates the surface Arctic oscillation index in the -60/+60 days time period around the ozone extreme event.
    
    #"nc_fid": file including variable
    # "var": name of zonal mean geopotential height variable
    # "lev": name of pressure level variable
    # "extreme_years": number of desired ozone extreme years.
    # "O3_highest_date": dates of ozone maxima
    # "O3_lowest_date": dates of ozone minima
    # "model": SOCOL, WACCM or MERRA
    # "O3_lowest_index": indices of the ozone minima
    # "O3_years": years in which ozone minima occurs
    
    manual=False # uses python function to calculate PC time series ("manual=False") or uses array multiplication ("manual=True")
    
    var=nc_fid[var]
    plev=nc_fid[lev]
    lats=nc_fid['lat']
    
    if lats[0]>0:
        lats_orientation = 'negative'
    else:
        lats_orientation = 'positive'

    # select surface variables (1000 hPa)    
    if plev[len(plev)-1] > 1000:
        var=var.interp(plev=100000)
        
    if plev[len(plev)-1] <= 1000 and model!='MERRA':   
        var=var.interp(plev=1000)
        
    if model=='MERRA':
        var=var.interp(lev=1000)  

    # interpolation to 2.5°
    var=var.interp(lat=np.linspace(-90,90,73))
    lats=lats.interp(lat=np.linspace(-90,90,73))
        
    lats=lats.sel(lat=slice(20,90)) #select latitudes between 20 and 90°N
    lats=np.array(lats)

    #calculate daily anomalies
    if model=='MERRA':
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time")
        #for MERRA, exclude year 2020 from the anomaly calculation
    else:
         var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")
 
    var_anomalies=var_anomalies.sel(lat=slice(20,90)) 

    # FIND INDICES OF OZONE EXTREMES   
    time=nc_fid['time']
    time=np.array(time)
    
    ozone_low_index=np.zeros((extreme_years))
    ozone_high_index=np.zeros((extreme_years))
    
    #select time indices when ozone extremes happen
    for year in range(extreme_years):
         ozone_low_index[year]=int(np.array(np.where(time==np.array(O3_lowest_date[year]))))
         ozone_high_index[year]=int(np.array(np.where(time==np.array(O3_highest_date[year]))))

    coslat = np.cos(np.deg2rad(lats).clip(0., 1.)) # latitudinal weights for PC calculation
    wgts = np.sqrt(coslat)[..., np.newaxis]
    
    gh_layer = np.array(var_anomalies)
 
    if len(np.shape(gh_layer)) == 3:
        gh_layer = np.reshape(gh_layer, (np.shape(gh_layer)[0], len(lats)))
        
    # calculate Arctic Oscillation index (1st EOF of all-year surface geopotential height anomalies projected onto daily anomalies)
        
    solver = Eof(gh_layer, weights=wgts, center=True)
    
    if manual==True: #if TRUE: the spatial EOF pattern is calculated using the python eofs function. 
        #The NAM time series is than calculated by projecting daily zonal mean geopotential height anomalies onto the 
        #EOF loading pattern to derive the principal component time series
        EOF=solver.eofs(neofs=1, eofscaling=2) #eofscaling=2 : EOFs are multiplied by the square-root of their eigenvalues.
    
        EOF=np.reshape((EOF), (nof_lats))
        EOF[np.isnan(EOF)]=0
        
        gh_layer[np.isnan(gh_layer)]=0
    
        weighted_EOF=EOF*coslat
        weighted_EOF[np.isnan(weighted_EOF)]=0
        
        AO = np.matmul(gh_layer, weighted_EOF)/(np.matmul(EOF.transpose(), weighted_EOF))
    
    if manual==False: #if FALSE: The PC time series (NAM indices) are derived directly via the python eofs function
        
        AO=np.reshape(solver.pcs(npcs=1, pcscaling=1), (np.shape(gh_layer)[0],))
    
    AO_xr = xr.DataArray(AO, coords=[time], dims=['time']) # AO xarray
    
    AO_xr_March=AO_xr.sel(time=nc_fid.time.dt.month.isin([3,4,5])) #select AO values in March, April, May
    AO_xr_March=AO_xr_March.groupby("time.year").mean("time")
    
    # get time frame -60 to +60 days after the ozone extremes 
    AO_min=np.zeros((extreme_years,120)) # AO index in the days-60/+60 around the ozone minima
    
    gh_layer_min=np.zeros((extreme_years,120, len(lats)))
    
    for i in range(extreme_years):   
        AO_min[i,:]=AO[int(ozone_low_index[i]-59):int(ozone_low_index[i]+61)]
        gh_layer_min[i,:,:]=gh_layer[int(ozone_low_index[i]-59):int(ozone_low_index[i]+61), :]
    
    # check, if sign of PC is correct
    nof_lats = len(lats)
    max_AO = AO.argmax()
    
    if lats_orientation=='negative':
               
            if np.nanmean(gh_layer[max_AO,0:int(nof_lats/2)-1])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+1:nof_lats-1])<0:
                    AO_min[:]=-AO_min[:]        
    else:  
              if np.nanmean(gh_layer[max_AO,0:int(nof_lats/2)-3])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+3:nof_lats-1])<0:
                      AO_min[:]=-AO_min[:]

                    
    print('Mean AO index difference 30 days after max. ozone depletion for low ozone years:  '   + str(np.mean(AO_min[:,60:90])) +' ± ' + str(np.std(-AO_min[:,60:90])))  
    print('Mean AO index difference 60 days after max. ozone depletion for low ozone years:  '   + str(np.mean(AO_min[:,60:119])) +' ± ' + str(np.std(-AO_min[:,60:119]))) 
    
    
    #This function returns: 
    # "AO_min": AO Index for all ozone minimum years in the ±60 days around the ozone minimum
    # "AO_xr": whole AO index time series as xarray
    # "AO_xr_March": March to May AO Indices of all years as xarray
    

    return AO_min, AO_xr, AO_xr_March

#___________________________________________________________________________



def calculate_NAM(nc_fid, var, lev, extreme_years, O3_highest_date, O3_lowest_date, model, O3_lowest_index, O3_years):
     
    # This function calculates the NAM index in the -30/+60 days time period around the ozone extreme event on all levels.
    
    #"nc_fid": file including variable
    # "var": name of zonal mean geopotential height variable
    # "lev": name of pressure level variable
    # "extreme_years": number of desired ozone extreme years.
    # "O3_highest_date": dates of ozone maxima
    # "O3_lowest_date": dates of ozone minima
    # "model": SOCOL, WACCM or MERRA
    # "O3_lowest_index": indices of the ozone minima
    # "O3_years": years in which ozone minima occurs
    
    manual=False # uses python function to calculate PC time series ("manual=False") or uses array multiplication ("manual=True")
    
    var=nc_fid[var]
    plev=nc_fid[lev]
    lats=nc_fid['lat']

    if lats[0]>0:
        lats_orientation = 'negative'
    else:
        lats_orientation = 'positive'
        
    # interpolate to 2.5°

    var=var.interp(lat=np.linspace(-90,90,73))
    lats=lats.interp(lat=np.linspace(-90,90,73))

    lats=lats.sel(lat=slice(20,90))
    lats=np.array(lats)
    
    # select pressure levels between 1 and 1000 hPa

    if plev[len(plev)-1] <= 1000 and model!='MERRA':
          plev=plev.sel(plev=slice(1,1000))
          var=var.sel(plev=slice(1,1000)) 
          plev=np.array(plev)
          
    
    if plev[len(plev)-1] > 1000:
          var=var.sel(plev=slice(100,100000))
          plev=plev.sel(plev=slice(100,100000))
          plev=np.array(plev)/100

    if model=='MERRA':
        plev=plev.sel(lev=slice(1000,1))
        var=var.sel(lev=slice(1000,1)) 
        plev=np.array(plev)
        
        #calculate daily anomalies
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time")
    else:
         var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")

    var_anomalies=var_anomalies.sel(lat=slice(20,90)) #select latitudes between 20 and 90°N

    # FIND INDICES OF OZONE EXTREMES   
    time=nc_fid['time']
    
    time=np.array(time)
    
    ozone_low_index=np.zeros((extreme_years))
    ozone_high_index=np.zeros((extreme_years))
    
    for year in range(extreme_years):
         ozone_low_index[year]=int(np.array(np.where(time==np.array(O3_lowest_date[year]))))
         ozone_high_index[year]=int(np.array(np.where(time==np.array(O3_highest_date[year]))))
    
    coslat = np.cos(np.deg2rad(lats).clip(0., 1.)) # latitudinal weights for PC calculation
    wgts = np.sqrt(coslat)[..., np.newaxis]

    # get geopotential height of polar cap 

    var_anomalies=np.array(var_anomalies)
    if len(np.shape(var_anomalies)) == 4:
        var_anomalies = np.reshape(var_anomalies, (np.shape(var_anomalies)[0], len(plev), len(lats)))

    AO = np.zeros((len(time), len(plev)))
    nof_lats = len(lats)

    for level in range(len(plev)): # calculate PC time series for each level
        
        gh_layer = np.array(var_anomalies[:,level,:])  
        solver = Eof(gh_layer, weights=wgts, center=True)
        
        if manual==True: #if TRUE: the spatial EOF pattern is calculated using the python eofs function. 
        #The NAM time series is than calculated by projecting daily zonal mean geopotential height anomalies onto the 
        #EOF loading pattern to derive the principal component time series
            EOF=solver.eofs(neofs=1, eofscaling=2) #eofscaling=2 : EOFs are multiplied by the square-root of their eigenvalues.
        
            EOF=np.reshape((EOF), (nof_lats))
            EOF[np.isnan(EOF)]=0
            
            gh_layer[np.isnan(gh_layer)]=0
            weighted_EOF=EOF*coslat
            weighted_EOF[np.isnan(weighted_EOF)]=0

            AO= np.matmul(gh_layer, weighted_EOF)/(np.matmul(EOF.transpose(), weighted_EOF))
        
        if manual==False:  #if FALSE: The PC time series (NAM indices) are derived directly via the python eofs function
            
            AO[:,level] =np.reshape(solver.pcs(npcs=1, pcscaling=1), (np.shape(gh_layer)[0],))
            
        max_AO = AO[:,level].argmax()
        
                    
         # check sign of PCs
        
        if model=='WACCM':

            if lats_orientation=='negative':
                    if np.nanmean(gh_layer[max_AO,1:int(nof_lats/2)-2])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+1:nof_lats-2])>0:
                            AO[:,level]=-AO[:,level]
            else:
                      if np.nanmean(gh_layer[max_AO,1:int(nof_lats/2)-2])-np.nanmean(gh_layer[max_AO,int(nof_lats/2)+1:nof_lats-2])<0:
                            AO[:,level]=-AO[:,level]
                      

    # get time frame -60 to +60 days after the ozone extremes 
         
    AO_min=np.zeros((extreme_years,90,len(plev)))
    gh_layer_min=np.zeros((extreme_years,90,len(plev),len(lats)))
     
    for i in range(extreme_years):   
        AO_min[i,:,:]=AO[int(ozone_low_index[i]-29):int(ozone_low_index[i]+61)]
        gh_layer_min[i,:,:,:]=var_anomalies[int(ozone_low_index[i]-29):int(ozone_low_index[i]+61)]

    #calculate composite mean
    AO_min=np.mean(AO_min, axis=0)
    gh_layer_min=np.mean(gh_layer_min, axis=0)
    
    
    if model=='SOCOL':
        for layer in range(len(plev)):      
       
            if lats_orientation=='negative':
                    if np.nanmean(gh_layer_min[40,layer, 0:int(nof_lats/2)-3])-np.nanmean(gh_layer_min[40, layer, int(nof_lats/2)+3:nof_lats-2])>0:
                        if AO_min[40,layer]<0:
                            AO_min[:,layer]=-AO_min[:,layer]
                        else: continue
                    
                    if np.nanmean(gh_layer_min[40,layer, 0:int(nof_lats/2)-3])-np.nanmean(gh_layer_min[40, layer, int(nof_lats/2)+3:nof_lats-2])<0:
                        if AO_min[40,layer]>0:
                            AO_min[:,layer]=-AO_min[:,layer]
                        else: continue
                  
                   
            else:
                if np.mean(gh_layer_min[40, layer, 0:int(nof_lats/2)-3])-np.mean(gh_layer_min[40, layer, int(nof_lats/2)+3:nof_lats-2])<0:
                       
                        if AO_min[40,layer]>0:
                          AO_min[:,layer]=-AO_min[:,layer]
                        else: continue
                        
                        
                if np.mean(gh_layer_min[40, layer, 0:int(nof_lats/2)-3])-np.mean(gh_layer_min[40,  layer, int(nof_lats/2)+3:nof_lats-2])>0:
                         
                        if AO_min[40,layer]<0:
                          AO_min[:,layer]=-AO_min[:,layer]
                        else: continue    
    
    AO_xr = xr.DataArray(AO, coords=[time, plev], dims=['time', 'plev'])    # save AO array as xarray

    #calculate significance based on a bootstrapping test
    if model=='SOCOL' or model=='MERRA':
        significance = bootstrapping_leap_alt_time(nc_fid, AO_xr, O3_lowest_index, AO_min, O3_years,model)
    if model=='WACCM':
        significance=bootstrapping_leap_alt_time(nc_fid, np.reshape( np.array(AO_xr), (200,365,len(plev))), O3_lowest_index, AO_min, O3_years,model)
    
    
    print('Mean AO index difference 30 days after max. ozone depletion for low ozone years at 1000 hPa:  '   + str(np.mean(AO_min[30:60, np.argwhere(plev==1000)])))  
    print('Mean AO index difference 30 days after max. ozone depletion for low ozone years at 500 hPa:  '   + str(np.mean(AO_min[30:60, np.argwhere(plev==500)])))  
    
    #This function returns: 
    # "AO_min": NAM Index for all ozone minimum years at all levels in the ±60 days around the ozone minimum
    # "AO_xr": whole NAM index time series as xarray
    # "AO_xr_March": March to May NAM Indices of all years as xarray

    
    return AO_min, AO_xr, significance



#__________________________________________________________________________________________________________


def plot_2D_anomalies(nc_fid, var, lev, extreme_years, O3_highest_date, O3_lowest_date, O3_lowest_index, O3_years,model):

   
    #This function calculates 2D anomalies (lat-lon) in the 30 days after the ozone extreme at 10 hPa averaged over all high/low ozone years
    
    #"nc_fid": file including variable
    # "var": name of zonal mean geopotential height variable
    # "lev": name of pressure level variable
    # "extreme_years": number of desired ozone extreme years.
    # "O3_highest_date": dates of ozone maxima
    # "O3_lowest_date": dates of ozone minima
    # "O3_lowest_index": indices of the ozone minima
    # "O3_years": years in which ozone minima occurs
    # "model": SOCOL, WACCM or MERRA
    
    factor=1000000 # the variables will be multiplied by this factor
    
    var=nc_fid[var]*factor
    
    lats=nc_fid['lat']
    lats=np.array(lats)
    
    lons=nc_fid['lon']
    lons=np.array(lons)

    time=nc_fid['time']
    time=np.array(time)
    plev=nc_fid[lev]
    
    #select 10 hPa level

    if plev[len(plev)-1] > 1000:
        var=var.interp(plev=1000)
        
    if plev[len(plev)-1] <= 1000 and model!='MERRA':   
        var=var.interp(lev=10)
        
    if model=='MERRA':
        var=var.interp(lev=10)  

    # Calculate anomalies    
    if model=='MERRA':
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time")
        #for MERRA, exclude year 2020 from the anomaly calculation
    else:
        var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")

    # FIND INDICES OF OZONE EXTREMES   

    ozone_low_date=np.empty((extreme_years))
    ozone_high_date=np.empty((extreme_years))

    for year in range(extreme_years):

        if np.array(O3_lowest_date[year]) in time:
           
            ozone_low_date[year]=int(np.reshape(np.array((np.where(time==np.array(O3_lowest_date[year])))), (1,)))
            ozone_high_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_highest_date[year]))), (1,)))

    # get data following the 30 days after an ozone extreme event
    
    var_min=np.empty((extreme_years,30,len(lats),len(lons)))
    var_max=np.empty((extreme_years,30,len(lats),len(lons)))
     
    for i in range(extreme_years):   
        if np.array(O3_lowest_date[i]) in time:
                var_min[i,:,:,:]=var_anomalies[int(ozone_low_date[i]):int(ozone_low_date[i]+30),:,:]
                var_max[i,:,:,:]=var_anomalies[int(ozone_high_date[i]):int(ozone_high_date[i]+30),:,:]
         
    # average over the 30-day period
    var_min = np.nanmean(var_min[:,:,:,:], axis=1)
    var_max = np.nanmean(var_max[:,:,:,:], axis=1)
   
    #average over all years
    var_min_zm = np.nanmean(var_min, axis=0)
    var_max_zm = np.nanmean(var_max, axis=0)
    
    #calculate significance of anomalies based on a bootstrapping test
    #significance_low = bootstrapping_leap(nc_fid, var_anomalies, O3_lowest_index, var_min_zm, O3_years )
    #significance_low = bootstrapping(np.reshape(np.array(var_anomalies), (200,365, len(lats), len(lons))), O3_lowest_index, var_min_zm, 200)
    
    # interpolate the values (to get values at the poles)
    
    if 90 and -90 not in lats:
        lats_90, lons_360, var_min_zm  = interpolate_pole(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_pole(0, lons, lats, var_max_zm)
      #  lats_90, lons_360, p_array = interpolate_pole(0, lons, lats, p_array)
      #  lats_90, lons_360, p_array_min = interpolate_pole(0, lons, lats, p_array_min)
        #lats_90, lons_360, significance_low = interpolate_pole(0, lons, lats, significance_low)
      #  lats_90, lons_360, p_array_max = interpolate_pole(0, lons, lats, p_array_max)
       # lats_90, lons_360, var_anomalies = interpolate_pole(0, lons, lats, var_anomalies)
        
    else:
        lats_90, lons_360, var_min_zm  = interpolate_lons(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_lons(0, lons, lats, var_max_zm)
      #  lats_90, lons_360, p_array = interpolate_lons(0, lons, lats, p_array)
      #  lats_90, lons_360, p_array_min = interpolate_lons(0, lons, lats, p_array_min)
       # lats_90, lons_360, significance_low = interpolate_lons(0, lons, lats, significance_low)
     #   lats_90, lons_360, p_array_max = interpolate_lons(0, lons, lats, p_array_max)
      #  lats_90, lons_360, var_anomalies = interpolate_lons(0, lons, lats, var_anomalies)
    
    m = Basemap(projection='ortho', lat_0=90,lon_0=0,resolution='l')
    lons_m,lats_m=np.meshgrid(lons_360,lats_90)
    xpt,ypt = m(lons_m,lats_m)
    
    
    #This function returns:
    #"var_min_zm": mean response in the 30 days after the ozone minimum
    #"var_max_zm": mean response in the 30 days after the ozone maximum
    # "xpt/ypt": coordinates for plotting
    # "var_anomalies": time series of anomalies of the variable of interest

    return var_min_zm, var_max_zm, xpt,ypt, var_anomalies

#__________________________________________________________________________________________________________


def plot_2D_anomalies_zonal_asym(nc_fid, var, lev, extreme_years, O3_highest_date, O3_lowest_date, O3_lowest_index, O3_years,model):
    
    # this function calculates zonal mean anomalies around the ozone extremes
    # either deviation from the zona mean climatology (var_zm_diff)
    # or daily zonal mean anomaly (daily deviation fromt the zonal mean value on the same day (var_zm_anomalies)
    # "var": name of variable to be plotted
    # "model": SOCOL, WACCM or MERRA
    # "lev": name of pressure level variable
   
    factor=1 # the variable will be mutiplied by this factor
    
    var=nc_fid[var]*factor
    lats=nc_fid['lat']
    lats=np.array(lats)
    
    lons=nc_fid['lon']
    lons=np.array(lons)

    time=nc_fid['time']
    time=np.array(time)


    if model=='MERRA':
        var_zm_anomalies = (var - var.mean(dim='lon'))/var.mean(dim='lon')
        zm_clim=var.groupby("time.dayofyear").mean("time").mean(dim='lon')
        var_zm_diff = var.groupby("time.dayofyear") - zm_clim
        var_anomalies=var.groupby("time.dayofyear")-var.groupby("time.dayofyear").mean("time")
    else:
        zm_clim=var.groupby("time.dayofyear").mean("time").mean(dim='lon')
        var_zm_diff = var.groupby("time.dayofyear") - zm_clim
        var_zm_anomalies = (var -var.mean(dim='lon'))/var.mean(dim='lon')
        var_anomalies=var.groupby("time.dayofyear")-var.groupby("time.dayofyear").mean("time")
     
    # FIND INDICES OF OZONE EXTREMES   
    
    ozone_low_date=np.empty((extreme_years))
    ozone_high_date=np.empty((extreme_years))
    
    for year in range(extreme_years):

        if np.array(O3_lowest_date[year]) in time:
            ozone_low_date[year]=int(np.reshape(np.array((np.where(time==np.array(O3_lowest_date[year])))), (1,)))
            ozone_high_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_highest_date[year]))), (1,)))

    # get data following the -5 to 5 days after an ozone extreme event
    
    var_min=np.empty((extreme_years,10,len(lats),len(lons)))
    var_min_zm_diff=np.empty((extreme_years,10,len(lats),len(lons)))
    var_max=np.empty((extreme_years,10,len(lats),len(lons)))
    var_min_anomalies=np.empty((extreme_years,10,len(lats),len(lons)))
     
    for i in range(extreme_years):   
        if np.array(O3_lowest_date[i]) in time:
                var_min[i,:,:,:]=var_zm_anomalies[int(ozone_low_date[i])-5:int(ozone_low_date[i]+5),:,:]
                var_max[i,:,:,:]=var_zm_anomalies[int(ozone_high_date[i])-5:int(ozone_high_date[i]+5),:,:]
                var_min_zm_diff[i,:,:,:]=var_zm_diff[int(ozone_low_date[i])-5:int(ozone_low_date[i]+5),:,:]
                var_min_anomalies[i,:,:,:]=var_anomalies[int(ozone_low_date[i])-5:int(ozone_low_date[i]+5),:,:]
                
                
    # average over the 10-day period

    var_min = np.nanmean(var_min[:,:,:,:], axis=1)
    var_max = np.nanmean(var_max[:,:,:,:], axis=1)
    var_min_zm_diff = np.nanmean(var_min_zm_diff[:,:,:,:], axis=1)
    var_min_anomalies = np.nanmean(var_min_anomalies[:,:,:,:], axis=1)
    
    var_min_zm = np.nanmean(var_min, axis=0)
    var_max_zm = np.nanmean(var_max, axis=0)
    var_min_zm_diff = np.nanmean(var_min_zm_diff, axis=0)
    var_min_anomalies = np.nanmean(var_min_anomalies, axis=0)
    
    
    # interpolate the values (to get values at the poles)
    
    if 90 and -90 not in lats:
        lats_90, lons_360, var_min_zm  = interpolate_pole(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_pole(0, lons, lats, var_max_zm)
        lats_90, lons_360, var_min_zm_diff  = interpolate_pole(0, lons, lats, var_min_zm_diff)
        lats_90, lons_360, var_min_anomalies  = interpolate_pole(0, lons, lats, var_min_anomalies)

        
    else:
        lats_90, lons_360, var_min_zm  = interpolate_lons(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_lons(0, lons, lats, var_max_zm)
        lats_90, lons_360, var_min_zm_diff  = interpolate_lons(0, lons, lats, var_min_zm_diff)
        lats_90, lons_360, var_min_anomalies  = interpolate_lons(0, lons, lats, var_min_anomalies)

    m = Basemap(projection='ortho', lat_0=90,lon_0=0,resolution='l')
    lons_m,lats_m=np.meshgrid(lons_360,lats_90)
    xpt,ypt = m(lons_m,lats_m)

    return var_min_zm, var_max_zm, var_min_zm_diff, var_min_anomalies, xpt,ypt, var_zm_anomalies
