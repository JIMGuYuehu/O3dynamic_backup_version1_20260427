import numpy as np
import xarray as xr
import random
from mpl_toolkits.basemap import Basemap

def find_ozone_extremes_baseonpressure(nc_fid, var, lev, years, extreme_years, model, p_top=30, p_bottom=100):
    """
    查找特定压力层范围内臭氧极值年份
    
    参数:
    nc_fid: 输入文件
    var: 臭氧变量名
    lev: 压力层变量名
    years: 输入文件中的年份数
    extreme_years: 需要的极值年份数
    model: 模式名称 (SOCOL, WACCM 或 MERRA)
    p_top: 顶层压力 (hPa)
    p_bottom: 底层压力 (hPa)
    """
    
    partial_column=True  # 计算部分臭氧柱
    
    O3=nc_fid[var]
    plev=nc_fid[lev]
    time=nc_fid['time']
    
    O3=O3.interp(lat=np.linspace(-90,90,73))     # 插值到2.5度
    
    if partial_column ==True:
        # 计算指定压力层范围的部分臭氧柱
        delta_p = np.zeros((len(O3.time), len(plev)))
        m_air = 28.964/(6.022E23)
        g = 980.6
            
        if plev[len(plev)-1] <= 1000 and model!='MERRA': 
            for level in range(1,len(plev)):
                delta_p[:,level].fill(plev[level] - plev[level-1])

            O3=O3.sel(lat=slice(60,90)) 
            weights = np.cos(np.deg2rad(O3.lat))
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')
            
            weights_p = xr.DataArray(delta_p*100, dims=['time','plev'], coords=[time,plev])
            O3 = O3 * weights_p * 10/ (g * m_air)
            
            # 使用输入的压力层范围
            O3=O3.sel(plev=slice(p_top,p_bottom)) 
            O3 = O3.sum(dim='plev')
            O3 = O3/2.687E16  # 转换为DU
            plev=plev.sel(plev=slice(p_top,p_bottom))
            
        if plev[len(plev)-1] > 1000:
            for level in range(1,len(plev)):
                delta_p[:,level].fill(plev[level] - plev[level-1])  
                
            O3=O3.sel(lat=slice(60,90))
            weights = np.cos(np.deg2rad(O3.lat))
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')
            
            weights_p = xr.DataArray(delta_p, dims=['time','plev'], coords=[time,plev])
            O3 = O3 * weights_p * 10/ (g * m_air)
            
            # 调整压力层范围（单位转换）
            O3=O3.sel(plev=slice(p_top*100, p_bottom*100))
            O3 = O3.sum(dim='plev')
            O3 = O3/2.687E16
            plev=plev.sel(plev=slice(p_top*100, p_bottom*100))
            plev=plev/100
            
        if model=='MERRA':
            for level in range(0,len(plev)-1):
                delta_p[:,level].fill(plev[level] - plev[level+1])   
            
            O3=O3.sel(lat=slice(60,90))*28.970/47.9982
            weights = np.cos(np.deg2rad(O3.lat))
            O3 = O3.weighted(weights)     
            O3=O3.mean(dim='lat')
            
            weights_p = xr.DataArray(delta_p*100, dims=['time','lev'], coords=[time,plev])
            O3 = O3 * weights_p * 10/ (g * m_air)
            
            # MERRA模式的压力层顺序是相反的
            O3=O3.sel(lev=slice(p_bottom,p_top)) 
            O3 = O3.sum(dim='lev')
            O3 = O3/2.687E16
            plev=plev.sel(lev=slice(p_bottom,p_top))
        
    if partial_column==False:

        # average over polar cap. Choose values at 100hPa
        if plev[len(plev)-1] > 1000:
            O3=O3.sel(lat=slice(60,90),plev=10000)
            
        if plev[len(plev)-1] <= 1000 and model!='MERRA':
            O3=O3.sel(lat=slice(60,90),plev=100) 
            
        if model=='MERRA':
            O3=O3.sel(lat=slice(60,90),lev=100)
     
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
  

def calc_parc_o3(var, p_top=30, p_bottom=100):
    """
    计算指定压力层范围的臭氧柱
    
    参数:
    var: xarray.DataArray, 臭氧数据
    p_top: float, 顶部压力层 (hPa)
    p_bottom: float, 底部压力层 (hPa)
    """
    m_air = 28.964/(6.022E23)
    g = 980.6

    if 'plev' in var.dims: 
        plev=var.plev
        level='plev'
    if 'lev' in var.dims: 
        plev=var.lev
        level='lev'
    if 'level' in var.dims: 
        plev=var.level
        level='level'
    
    time=var.time
    delta_p = np.zeros((len(time), len(plev)))
    
    if plev[0]<plev[len(plev)-1] and plev[len(plev)-1] <= 1000: 
        factor=100
        factor_2=1 # conversion Pa->hPa
    if plev[0]>plev[len(plev)-1] and plev[0] <=1000: 
        factor=100
        factor_2=1
    if plev[0]<plev[len(plev)-1] and plev[len(plev)-1] >1000: 
        factor=1
        factor_2=100
    if plev[0]>plev[len(plev)-1] and plev[0] >1000: 
        factor=1
        factor_2=100
    
    if plev[0]<plev[len(plev)-1]: # for pressure levels in hPa
        for levelx in range(1,len(plev)): 
            delta_p[:,levelx].fill(plev[levelx] - plev[levelx-1])

        weights_p = xr.DataArray(delta_p*factor, dims=['time',level], coords=[time,plev])
 
        O3 = var * weights_p * 10/ (g * m_air)
    
        if level=='plev': O3=O3.sel(plev=slice(p_top*factor_2, p_bottom*factor_2)) 
        if level=='lev': O3=O3.sel(lev=slice(p_top*factor_2, p_bottom*factor_2))
        if level=='level': O3=O3.sel(level=slice(p_top*factor_2, p_bottom*factor_2))

        O3 = O3.sum(dim=level)
        O3 = O3/2.687E16  #calculate DU
        
    if plev[0]>plev[len(plev)-1]: # for pressure levels in hPa
        for levelx in range(0,len(plev)-1): 
            delta_p[:,levelx].fill(plev[levelx] - plev[levelx+1])

        weights_p = xr.DataArray(delta_p*factor, dims=['time',level], coords=[time,plev])

        O3 = var * weights_p * 10/ (g * m_air)
        
        if level =='plev': O3=O3.sel(plev=slice(p_bottom*factor_2, p_top*factor_2)) 
        if level=='lev': O3=O3.sel(lev=slice(p_bottom*factor_2, p_top*factor_2)) 
        if level=='level': O3=O3.sel(level=slice(p_bottom*factor_2, p_top*factor_2)) 
            
        O3 = O3.sum(dim=level)
        O3 = O3/2.687E16  #calculate DU

    return O3.where(O3 != 0)

def analyse_ozone_extremes(nc_fid, var, extreme_years, surface_pressure,  O3_highest_date, O3_lowest_date, O3_lowest_index, O3_years,model):

    # This function calculates and plots the surface impact (pressure, temperature or precipitation) for high and low 
    # ozone years.
    
    # It requires input in form of daily surface fiels (lat - lon) of the variable "var".
    # "extreme_years": number of desired ozone extreme years.
    # if "surface_pressure=True": the input variable is surface pressure in Pa. The variable will be divided by 100 to get result in hPa
    # "O3_highest_date": date of the ozone maxima
    # "O3_highest_date": date of the ozone minima
    # "O3_lowest_index": index (day after March 1st) when ozone minima occur
    # "O3_years": years in which ozone minima occur
    # "model": SOCOL, WACCM or MERRA
    
    var=nc_fid[var]
    
    lons=nc_fid['lon']
    lons=np.array(lons)
    lats=nc_fid['lat']
    lats=np.array(lats)
    time=nc_fid['time']
    time=np.array(time)
    
    # Calculate anomalies
  
    if model=='MERRA':
        var_anomalies = var.groupby("time.dayofyear") - var.sel(time=nc_fid.time.dt.year.isin(np.linspace(1980,2019,40))).groupby("time.dayofyear").mean("time") 
        #for MERRA, exclude year 2020 from the anomaly calculation
    else:
        var_anomalies = var.groupby("time.dayofyear") - var.groupby("time.dayofyear").mean("time")
       # var_anomalies=(var.groupby("time.dayofyear") - var_mean.mean("lon"))/var.mean("lon") #zonal asymmetries
        

    # FIND INDICES OF OZONE EXTREMES   

    var_anomalies_xr=var_anomalies
    var_anomalies=np.array(var_anomalies)
    
    ozone_low_date=np.zeros((extreme_years))
    ozone_high_date=np.zeros((extreme_years))

    for year in range(extreme_years):
         ozone_low_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_lowest_date[year]))), (1,)))
         ozone_high_date[year]=int(np.reshape(np.array(np.where(time==np.array(O3_highest_date[year]))), (1,)))

        
    # get data following the 30 days after an ozone extreme event
    
    var_min=np.zeros((extreme_years,60,len(lats),len(lons)))
    var_max=np.zeros((extreme_years,60,len(lats),len(lons)))
     
    for i in range(extreme_years):   
            var_min[i,:,:,:]=var_anomalies[int(ozone_low_date[i])-0:int(ozone_low_date[i]+60),:,:]
            var_max[i,:,:,:]=var_anomalies[int(ozone_high_date[i])-0:int(ozone_high_date[i]+60),:,:]
         
    # average over the 30-day period

    var_min = np.nanmean(var_min[:,:,:,:], axis=1)
    var_max = np.nanmean(var_max[:,:,:,:], axis=1)
    
    # calculate significance of difference of high and low ozone years based on a t-test
    
  #  t_array, p_array = scipy.stats.ttest_ind(var_max,var_min, axis=0, equal_var=False)
  #  t_array_min, p_array_min = scipy.stats.ttest_1samp(var_min,0, axis=0)
 #   t_array_max, p_array_max = scipy.stats.ttest_1samp(var_min,0, axis=0)
  
    var_min_zm = np.nanmean(var_min, axis=0)
    var_max_zm = np.nanmean(var_max, axis=0)

   # calculate significance based on a bootstrapping test
    """
    if model=='MERRA' or model=='SOCOL':
        significance_low = bootstrapping_leap(nc_fid, var_anomalies_xr, O3_lowest_index, var_min_zm, O3_years )

    if model=='WACCM':
        significance_low = bootstrapping(np.reshape(np.array(var_anomalies), (39,365, len(lats), len(lons))), O3_lowest_index, var_min_zm, 39) #WACCM does not include leap years
    """
    significance_low=[]
    
    # interpolate the values (to get values at the poles)
    
    if 90 and -90 not in lats:
        lats_90, lons_360, var_min_zm  = interpolate_pole(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_pole(0, lons, lats, var_max_zm)
    #    lats_90, lons_360, p_array = interpolate_pole(0, lons, lats, p_array)
     #   lats_90, lons_360, p_array_min = interpolate_pole(0, lons, lats, p_array_min)
      #  lats_90, lons_360, significance_low = interpolate_pole(0, lons, lats, significance_low)
     #   lats_90, lons_360, p_array_max = interpolate_pole(0, lons, lats, p_array_max)
       # lats_90, lons_360, var_anomalies = interpolate_pole(0, lons, lats, var_anomalies)
        
    else:
        lats_90, lons_360, var_min_zm  = interpolate_lons(0, lons, lats, var_min_zm)
        lats_90, lons_360, var_max_zm  = interpolate_lons(0, lons, lats, var_max_zm)
     #   lats_90, lons_360, p_array = interpolate_lons(0, lons, lats, p_array)
     #   lats_90, lons_360, p_array_min = interpolate_lons(0, lons, lats, p_array_min)
      #  lats_90, lons_360, significance_low = interpolate_lons(0, lons, lats, significance_low)
      #  lats_90, lons_360, p_array_max = interpolate_lons(0, lons, lats, p_array_max)
      #  lats_90, lons_360, var_anomalies = interpolate_lons(0, lons, lats, var_anomalies)
    
    
    #new coordinates
    
    m = Basemap(projection='ortho', lat_0=90,lon_0=0,resolution='l')
    lons_m,lats_m=np.meshgrid(lons_360,lats_90)
    xpt,ypt = m(lons_m,lats_m)
    
    #function returns:
    
    #"var_min_zm": mean surface response in the 30 days following the ozone minima (lat-lon)
    # "xpt" & "ypt": coordinates used for plotting
    # "var_anomalies_xr": xarray containing the surface anomalies
    # "significance_low": significance of the surface response in low ozone years (lat-lon)
    
    #the output for high ozone years can be added here

    return var_min_zm, xpt,ypt, var_anomalies_xr, significance_low


#_________________________________________________________________________________________   

def interpolate_lons(time, lons, lats, array):
    
    # This function interpolates an array to 360° longitude
    
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
    
    

#__________________________________________________________________________________