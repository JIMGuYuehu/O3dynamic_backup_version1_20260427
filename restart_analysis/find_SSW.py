from netCDF4 import Dataset
import numpy as np
from get_winter_values import get_winter_values
import math
import xarray as xr



# THIS PROGRAM CALCULATES THE SSW FREQUENCY BASED ON THE DEFINITION IN CHARLTON AND POLVANI, 2007 
# (WIND REVERSAL AT 10 hPA AND 60°N)
# THE PROGRAM ALSO CALCULATES THE STANDARD ERROR OF THE SSW FREQUENCY BASED ON CHARLTON, 2006



def find_nearest(array, value):
     
    index =np.argmin(np.abs(array-value))
    
    return(index)



    
def find_SSW(nc_fid, years, startyear, var, lev):

    # This function calculated the SSW frequency using numpy. This function can only be used for data without leap years.
    # Input data should contain zonally averaged wind
    
    #nc_fid: Name of input file
    #years: number of years in the record
    #startyear: startyear of the data
    #var: name of wind variable
    #lev: name of pressure level coordinate

    

    # get latitudes and pressure levels
    lats = nc_fid.variables['lat'][:] 
    plev = nc_fid.variables[lev][:]
    

    # check orientation of lats (south - north or north - south)
    if lats[0]>0:
        lats_orientation = 'negative'
    else:
        lats_orientation = 'positive'
    
    
    # get pressure level in hPa if in Pa 
    if plev[len(plev)-1] > 1000:
        plev=plev/100

    
    # find pressure level of 10hPa
    SSW_level = find_nearest(plev,10)


    # get wind variable
    wind_u = nc_fid.variables[var]  # (time, plev, lats)
    
    
    # delete longitudinal dimension and get data at 10 hPa
    if len(np.shape(wind_u))==4:
        wind_u_SSW = wind_u[0:365*years, SSW_level,:,0] # choose wind at 10 hPa
        
        
    if len(np.shape(wind_u))==3:
        wind_u_SSW = wind_u[0:365*years, SSW_level,:] # choose wind at 10 hPa    
    


    # extract wind at 10hPa and 60°N --> LINEAR INTERPOLATION if 60 °N not in lats
    # find latitudes  60 °N

    SSW_lat = find_nearest(60,lats)
    
    if 60 not in lats:
        
        if lats_orientation == 'positive':
            if lats[SSW_lat] > 60:
                SSW_lat_upper = SSW_lat
                SSW_lat_lower = SSW_lat - 1
            else:     
                SSW_lat_upper = SSW_lat + 1
                SSW_lat_lower = SSW_lat 
        else:     
            if lats[SSW_lat] > 60:
                SSW_lat_upper = SSW_lat
                SSW_lat_lower = SSW_lat + 1
            else:
                SSW_lat_upper = SSW_lat - 1
                SSW_lat_lower = SSW_lat

        wind_u_SSW_upper = wind_u_SSW[:, SSW_lat_upper]
        wind_u_SSW_lower = wind_u_SSW[:, SSW_lat_lower]
        
        diff_lat  = lats[SSW_lat_upper] - lats[SSW_lat_lower]
        diff_lower = 60 - lats[SSW_lat_lower]
     
        wind_u_SSW = wind_u_SSW_lower + (wind_u_SSW_upper - wind_u_SSW_lower) * diff_lower/diff_lat
    
    else: 
    
        wind_u_SSW=wind_u_SSW[:, SSW_lat]
    
    
    wind_u_SSW=np.reshape(wind_u_SSW, (years*365))
    
    
    #_________________________________________________________________________________________________________________
    
    
    # find winter wind and winter dates, save in 2D array
    
    wind_u_winter, winter_indices, winter_dates=get_winter_values(wind_u_SSW, years, startyear)
    

    #___________________________________________________________________________________________________________________          
                
   
    
    # SSW defined as wind at 10hPa and 60°N turning easterly (see Charlton and Polvani, 2007)
    
    SSW_dates=[]   # saves dates, on which wind reversal happens          
    SSW_indices_structured=[] # saves year and day of the winter (starting from November) on which wind reversal happens in 2D array (needed for MC statistical test)
    SSW_indices=[] # saves indices of wind array, on which wind reversal happens      
    count_westerly_days=0                    
    year=0
    day=0                    
    #day=30 #mid-winter SSW
    SSW = False
    
    
    count_SSWs_per_winter=0 # counts number of SSWs per winter and appends it to list
    SSWs_per_winter=[]
    
    
    while True:
        
        if wind_u_winter[year,day] < 0:   

            # check, if this is a final warming (if winds turn westerly for at least 10 consecutive days before April 30th)
            
            count_westerly_days=0 
            
            for i in range(day+1, 181): 
                if wind_u_winter[year,i]>0:
                    count_westerly_days=count_westerly_days+1
                    
                    if count_westerly_days==10:
                        SSW_dates.append(winter_dates[year, day])   
                        SSW_indices.append(winter_indices[year,day])
                        SSW_indices_structured.append([year,day])
                        count_SSWs_per_winter=count_SSWs_per_winter+1
                        SSW=True
                        break
                    
                if wind_u_winter[year,i]<0:
                    count_westerly_days=0    
                    
           # block days until wind was westerly again for 20 consecutive days to not count one SSW twice
           
            count_westerly_days_2=0
            next_SSW=0

            if SSW==True:
                
                for j in range(day+1, 152): 
                #for j in range(day+1, 121): #mid-winter SSW
                    if wind_u_winter[year,j]>0:
                        count_westerly_days_2=count_westerly_days_2+1
                        if count_westerly_days_2==20:
                            next_SSW=j-day
                            break
                        
                    if wind_u_winter[year,j]<0:
                        count_westerly_days_2=0    
           
        
            
        if SSW==True: 
           
            if count_westerly_days_2<20: 
                
                year=year+1
                day=0
                #day=30 # mid-winter warming
                SSWs_per_winter.append(count_SSWs_per_winter)
                count_SSWs_per_winter=0
                
            if count_westerly_days_2==20:
                    
                day=day+next_SSW      
                
                
        if SSW==False: 
            
            if day+1 <= 150: # jump to next day until last day of March
         #   if day+1 <= 119: #mid-winter SSW 
                day=day+1
                
            else:    
                year=year+1
                day=0 
                #day=30 # mid-winter SSW
                SSWs_per_winter.append(count_SSWs_per_winter)
                count_SSWs_per_winter=0

                
        SSW=False    
        
        if year>=years-1: break
    
    
    
    years=years-1 # equals number of winters
    
    
    #_______________________________________________________________________________________________
    
    # calculate standard error of SSW frequency (based on Charlton, 2006)
    
    zero_SSWs = SSWs_per_winter.count(0)
    one_SSWs = SSWs_per_winter.count(1)
    two_SSWs = SSWs_per_winter.count(2)
    three_SSWs = SSWs_per_winter.count(3)
    four_SSWs = SSWs_per_winter.count(4)
    
    number_SSWs = len(SSW_dates)
    
    # sample mean frequency
    x = 1 * (one_SSWs/years) + 2 * (two_SSWs/years) + 3 * (three_SSWs/years) + 4 * (four_SSWs/years)
    
    # sample variance

    s_squared = math.pow((0-x),2) * (zero_SSWs/years) + math.pow((1-x),2) * (one_SSWs/years) + math.pow((2-x),2) * (two_SSWs/years) + math.pow((3-x),2) * (three_SSWs/years) + math.pow((4-x),2) * (four_SSWs/years)
    
    # standard error of SSW frequency
    
    e = math.sqrt(s_squared)/math.sqrt(number_SSWs) 
    

    
    #_____________________________________________________________________________________________-

    # calculate average amount of SSWs in each month

    SSW_Nov=[s for s in SSW_dates if '.11.' in s]
    SSW_Dec=[s for s in SSW_dates if '.12.' in s]
    SSW_Jan=[s for s in SSW_dates if '.1.' in s]
    SSW_Feb=[s for s in SSW_dates if '.2.' in s]
    SSW_Mar=[s for s in SSW_dates if '.3.' in s]
    
        
    years=years-1
    print("Total amount of SSWs per winter:" + "%.3f" % (x) + " ± " + str(e))
    print("SSWs per November:" + "%.3f" % (len(SSW_Nov)/years))
    print("SSWs per December:" + "%.3f" % (len(SSW_Dec)/years))
    print("SSWs per January:" +  "%.3f" % (len(SSW_Jan)/years))
    print("SSWs per February:" + "%.3f" % (len(SSW_Feb)/years))
    print("SSWs per March:" + "%.3f" % (len(SSW_Mar)/years))
    
    return SSW_dates, SSW_indices, SSW_indices_structured
    




   
def find_SSW_leap(nc_fid, years, var, lev,leap, model):


    # This function calculated the SSW frequency using xarray. Data can therefore contain leap years.
    # Input data should contain zonally averaged zonal wind.
    
    #nc_fid: Name of input file
    #years: number of years in the record
    #var: name of wind variable
    #lev: name of pressure level coordinate
    #leap: True or False: Does the data contain leap years?
    

    # get latitudes and pressure levels
    plev=nc_fid[lev]
    time=nc_fid['time']
    
    
    # get wind variable
    wind_u = nc_fid[var]  # (time, plev, lats)
    
    
    plev=np.array(plev)
    
    wind_u=wind_u.interp(lat=np.linspace(-90,90,73))


    #extract wind at 10hPa and 60°N 

    if model=='MERRA':
        wind_u=wind_u.sel(lat=60,lev=10, method='nearest') 
    else:
        if plev[len(plev)-1] > 1000: # if pressure levels are in Pa
            wind_u=wind_u.sel(lat=60,plev=1000, method='nearest')
        else: # if pressure levels are in hPa
            wind_u=wind_u.sel(lat=60,plev=10, method='nearest') 


   # choose winter months (Nov-March)
    wind_u_winter=wind_u.sel(time=nc_fid.time.dt.month.isin([1,2,3,4,11,12]))
    wind_u_winter=wind_u_winter.groupby("time.year")
    
    time=time.sel(time=nc_fid.time.dt.month.isin([1,2,3,4,11,12]))
    time=time.groupby("time.year")
    
    #rearrange array to get chronological wind from Nov-March instead of Jan,Febr,March,Nov,Dec
    winter_values=np.zeros((years-1, 182))
    winter_time=[]

   
    if leap==True: # This is for data that includes leap years (SOCOL and MERRA2)
        i=0
        for year, group in wind_u_winter:
            wind_u_year=group[:]
            wind_u_year=np.roll(wind_u_year, 61)
            winter_values[i,0:61]=np.reshape(wind_u_year[0:61], (61))
            i=i+1 
            if i==years-1:
                break
        
        i=0
        for year, group in wind_u_winter:
            wind_u_year=group[:]
            wind_u_year=np.array(wind_u_year)

        #    wind_u_year=np.roll(wind_u_year, 61)
            if i>=1:
                if year%4==0:
                    winter_values[i-1,61:182]=np.reshape(wind_u_year[0:121], (121))
                        
                else:
                    winter_values[i-1,61:181]=np.reshape(wind_u_year[0:120], (120))
                  #  winter_values[i-1,120]=np.mean(wind_u_year[120:121])
                  #  winter_values[i-1,121:182]=wind_u_year[120:181,0]
            i=i+1 
            if i==years:
                break
    
    if leap==False: # This is for data that does not have leap years (WACCM)
        i=0
        for year, group in wind_u_winter:
            wind_u_year=group[:]
            wind_u_year=np.roll(wind_u_year, 61)
            winter_values[i,0:61]=np.reshape(wind_u_year[0:61], (61))
            i=i+1 
            if i==years-1:
                break
        
        i=0
        for year, group in wind_u_winter:
            wind_u_year=group[:]
            wind_u_year=np.roll(wind_u_year, 61)
            if i>=1:
                  winter_values[i-1,61:181]=np.reshape(wind_u_year[0:120], (120))
            
            i=i+1 
            if i==years:
                break

# =============================================================================
#     i=0
#     for year, group in time:
#         winter_time_year=np.array(group[:])
#         winter_time.append(np.array(winter_time_year[121:182]))
#         i=i+1 
#         if i==years-1:
#             break
#     
#     i=0
#     for year, group in time:
#         winter_time_year=group[:]
#     
#         if i>=1:
#             if year%4==0:
#                 winter_time.append(np.array(winter_time_year[0:121]))
#             else:
#                 winter_time.append(np.array(winter_time_year[0:59]))
#                 winter_time.append(np.nan)
#                 winter_time.append(np.array(winter_time_year[59:120]))
#         i=i+1 
#         if i==years:
#             break
# 
#     print(np.array(winter_time))
#     winter_time=np.reshape(np.array(winter_time), (years-1, 182))
#     print(winter_time)
# =============================================================================
        
    SSW_days=[]   # saves dates, on which wind reversal happens          
    SSW_indices_structured=[] # saves year and day of the winter (starting from November) on which wind reversal happens in 2D array (needed for MC statistical test)
    SSW_indices=[] # saves indices of wind array, on which wind reversal happens      
    count_westerly_days=0        
    SSW_years=[]
            
    year=0
    day=0                    
   # day=30 #mid-winter SSW
    SSW = False
    
    
    
    
    count_SSWs_per_winter=0 # counts number of SSWs per winter and appends it to list
    SSWs_per_winter=[]
    
    wind_u_winter=winter_values
    
    while True:
        
        if wind_u_winter[year,day] < 0:   

            # check, if this is a final warming (if winds turn westerly for at least 10 consecutive days before April 30th)
            
            count_westerly_days=0 
            
            for i in range(day+1, 182): 
                if wind_u_winter[year,i]>0:
                    count_westerly_days=count_westerly_days+1
                    
                    if count_westerly_days==10:
                        SSW_days.append(day)   
                        SSW_indices.append([year,day])
                     #   if day<61: SSW_years.append(year) # get only February and March SSWs
                     #   else:  SSW_years.append(year+1)
                        SSW_years.append(year+1)
                       # SSW_indices_structured.append([year,day])
                        count_SSWs_per_winter=count_SSWs_per_winter+1
                        SSW=True
                        break
                    
                if wind_u_winter[year,i]<0:
                    count_westerly_days=0    
                    
           # block days until wind was westerly for 20 consecutive days 
           
            count_westerly_days_2=0
            next_SSW=0

            if SSW==True:
                
                for j in range(day+1, 152): #mid-winter SSW
               # for j in range(day+1, 121): 
                    if wind_u_winter[year,j]>0:
                        count_westerly_days_2=count_westerly_days_2+1
                        if count_westerly_days_2==20:
                            next_SSW=j-day
                            break
                        
                    if wind_u_winter[year,j]<0:
                        count_westerly_days_2=0    
           
        
            
        if SSW==True: 
           
            if count_westerly_days_2<20: 
                
                year=year+1
                day=0
                #day=30 # mid-winter warming
                SSWs_per_winter.append(count_SSWs_per_winter)
                count_SSWs_per_winter=0
                
            if count_westerly_days_2==20:
                    
                day=day+next_SSW      
                
                
        if SSW==False: 
            
            if day+1 <= 150: # jump to next day until last day of March
         #   if day+1 <= 119: #mid-winter SSW 
                day=day+1
                
            else:    
                year=year+1
                day=0 
                #day=30 # mid-winter SSW
                SSWs_per_winter.append(count_SSWs_per_winter)
                count_SSWs_per_winter=0

                
        SSW=False    
        
        if year>=years-1: break
    
    
    
    
    years=years-1 # equals number of winters
    

    
    #_______________________________________________________________________________________________
    
    # calculate standard error of SSW frequency (based on Charlton, 2006)
    
    zero_SSWs = SSWs_per_winter.count(0)
    one_SSWs = SSWs_per_winter.count(1)
    two_SSWs = SSWs_per_winter.count(2)
    three_SSWs = SSWs_per_winter.count(3)
    four_SSWs = SSWs_per_winter.count(4)
    
    number_SSWs = len(SSW_indices)
    
    # sample mean frequency
    x = 1 * (one_SSWs/years) + 2 * (two_SSWs/years) + 3 * (three_SSWs/years) + 4 * (four_SSWs/years)
    
    # sample variance

    s_squared = math.pow((0-x),2) * (zero_SSWs/years) + math.pow((1-x),2) * (one_SSWs/years) + math.pow((2-x),2) * (two_SSWs/years) + math.pow((3-x),2) * (three_SSWs/years) + math.pow((4-x),2) * (four_SSWs/years)
    
    # standard error of SSW frequency
    
    e = math.sqrt(s_squared)/math.sqrt(number_SSWs) 
    
    
    #_____________________________________________________________________________________________-

    # calculate average amount of SSWs in each month
    
    SSW_days=np.array(SSW_days)

    SSW_Nov = ( SSW_days <= 30).sum()
    SSW_Dec= (SSW_days > 30).sum() - (SSW_days > 61).sum()
    SSW_Jan=(SSW_days > 61).sum() - (SSW_days > 92).sum()
    SSW_Feb=(SSW_days > 92).sum() - (SSW_days > 121).sum()
    SSW_Mar=(SSW_days > 121).sum()
    
        
    years=years-1
    print("Total amount of SSWs per winter:" + "%.3f" % (x) + " ± " + str(e))
    print("SSWs per November:" + "%.3f" % (SSW_Nov/years))
    print("SSWs per December:" + "%.3f" % (SSW_Dec/years))
    print("SSWs per January:" +  "%.3f" % (SSW_Jan/years))
    print("SSWs per February:" + "%.3f" % (SSW_Feb/years))
    print("SSWs per March:" + "%.3f" % (SSW_Mar/years))
    
    
  #  np.savetxt('SSW_years_MERRA.txt', SSW_years)
    
    
    return SSW_indices, SSW_indices_structured, np.array(SSW_years)

    
# =============================================================================
#     # SSW defined as wind at 10hPa and 60°N turning easterly (see Charlton and Polvani, 2007)
#     
#     SSW_dates=[]   # saves dates, on which wind reversal happens          
#     SSW_indices_structured=[] # saves year and day of the winter (starting from November) on which wind reversal happens in 2D array (needed for MC statistical test)
#     SSW_indices=[] # saves indices of wind array, on which wind reversal happens      
#     SSW_years=[]
#     SSW_days=[]
#     count_westerly_days=0                    
#     year=0
#     day=0                    
#     #day=30 #mid-winter SSW
#     SSW = False
#     
#     
#     count_SSWs_per_winter=0 # counts number of SSWs per winter and appends it to list
#     SSWs_per_winter=[]
#         
#     
#     for year, group in wind_u_winter:
# 
#         day=0
#         wind_u_year=group[:]
#  
#     
#         
#         wind_u_year=np.roll(wind_u_year, 61)
# 
# 
# 
#         while day<=len(wind_u_year)-30:
#       
#             if wind_u_year[day] < 0 and day<len(wind_u_year)-30:   
#     
#           
#                 # check, if this is a final warming (if winds turn westerly for at least 10 consecutive days before April 30th)
#                 
#                 count_westerly_days=0 
#     
#                 for i in range(day+1, len(wind_u_year)): 
#                     if wind_u_year[i]>0:
#                         count_westerly_days=count_westerly_days+1
#                         
#                         if count_westerly_days==10:
#                          #   SSW_dates.append(winter_dates[year, day])   
#                             SSW_indices.append([year,day])
#                             SSW_years.append(year)
#                             SSW_days.append(day)    
#                         
#                          #   SSW_indices_structured.append([year,day])
#                             count_SSWs_per_winter=count_SSWs_per_winter+1
#                             SSW=True
#                 
#                             break
#                         
#                     if wind_u_year[i]<0:
#                         count_westerly_days=0    
#                         
#                # block days until wind was westerly for 20 consecutive days 
#                
#                 count_westerly_days_2=0
#                 next_SSW=0
#     
#                 if SSW==True:
#                     
#                     #for j in range(day+1, 152): mid-winter SSW
#                     for j in range(day+1, len(wind_u_year)-30): 
#                         if wind_u_year[j]>0:
#                             count_westerly_days_2=count_westerly_days_2+1
#                             if count_westerly_days_2==20:
#                                 next_SSW=j-day
#                                 break
#                             
#                         if wind_u_year[j]<0:
#                             count_westerly_days_2=0    
#                
#             
#                 
#             if SSW==True: 
#                
#                 if count_westerly_days_2<20: 
#                     
#                     #year=year+1
#                     day=len(wind_u_year)
#                    # day=30 # mid-winter warming
#                     SSWs_per_winter.append(count_SSWs_per_winter)
#                     count_SSWs_per_winter=0
#                     
#                 if count_westerly_days_2==20:
#                         
#                     day=day+next_SSW      
#                     
#                     
#             if SSW==False: 
#                 
#                 if day+1 <= len(wind_u_year)-30: # jump to next day until last day of March
#               #  if day+1 <= 119: #mid-winter SSW 
#                     day=day+1
#                     
#                 else:    
#     
#                     day=len(wind_u_year)
#                    # day=30 # mid-winter SSW
#                     SSWs_per_winter.append(count_SSWs_per_winter)
#                     count_SSWs_per_winter=0
#     
#                     
#             SSW=False    
#         
#       #  if year>=years-1: break
#     
#     
#  
#     years=years-1 # equals number of winters
#     
#     
#     
#     # get SSW dates 
#  
#     time=time.sel(time=nc_fid.time.dt.month.isin([11,12,1,2,3,4]))
#     time=time.groupby('time.year')
#     
#     
#     SSW_dates=[]
#     i=0
#     for year, group in time:
#         
#         time_SSW=group[:]
#         time_SSW=np.roll(time_SSW, 61)
#         
#      
#         
#         if year in SSW_years:
# 
#             if SSW_years.count(year)==1:
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i]]))
#                 i=i+1  
#             if SSW_years.count(year)==2:
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i]]))
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i+1]]))
#                 i=i+2
#             if SSW_years.count(year)==3:
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i]]))
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i+1]]))
#                 SSW_dates.append(np.array(time_SSW[SSW_days[i+2]]))
#                 i=i+3     
# 
#             
# 
#     
#     #_______________________________________________________________________________________________
#     
#     # calculate standard error of SSW frequency (based on Charlton, 2006)
#     
#     zero_SSWs = SSWs_per_winter.count(0)
#     one_SSWs = SSWs_per_winter.count(1)
#     two_SSWs = SSWs_per_winter.count(2)
#     three_SSWs = SSWs_per_winter.count(3)
#     four_SSWs = SSWs_per_winter.count(4)
#     
#     number_SSWs = len(SSW_years)
#     
#     # sample mean frequency
#     x = 1 * (one_SSWs/years) + 2 * (two_SSWs/years) + 3 * (three_SSWs/years) + 4 * (four_SSWs/years)
#     
#     # sample variance
# 
#     s_squared = math.pow((0-x),2) * (zero_SSWs/years) + math.pow((1-x),2) * (one_SSWs/years) + math.pow((2-x),2) * (two_SSWs/years) + math.pow((3-x),2) * (three_SSWs/years) + math.pow((4-x),2) * (four_SSWs/years)
#     
#     # standard error of SSW frequency
#     
#     e = math.sqrt(s_squared)/math.sqrt(number_SSWs) 
#     
#     SSW_dates=xr.DataArray(SSW_dates, coords=[SSW_dates], dims=['time'], name='time')
# 
# 
#     
#     time=SSW_dates
#     
#     #_____________________________________________________________________________________________-
# 
#     # calculate average amount of SSWs in each month
# 
# 
# 
#     SSW_Nov=np.array(SSW_dates.sel(time=time.dt.month.isin([11])))
#     SSW_Dec=np.array(SSW_dates.sel(time=time.dt.month.isin([12])))
#     SSW_Jan=np.array(SSW_dates.sel(time=time.dt.month.isin([1])))
#     SSW_Feb=np.array(SSW_dates.sel(time=time.dt.month.isin([2])))
#     SSW_Mar=np.array(SSW_dates.sel(time=time.dt.month.isin([3])))
#     
#         
#     
#     
#     years=years-1
#     print("Total amount of SSWs per winter:" + "%.3f" % (x) + " ± " + str(e))
#     print("SSWs per November:" + "%.3f" % (len(SSW_Nov)/years))
#     print("SSWs per December:" + "%.3f" % (len(SSW_Dec)/years))
#     print("SSWs per January:" +  "%.3f" % (len(SSW_Jan)/years))
#     print("SSWs per February:" + "%.3f" % (len(SSW_Feb)/years))
#     print("SSWs per March:" + "%.3f" % (len(SSW_Mar)/years))
#     
# =============================================================================
  #return SSW_dates, SSW_indices, SSW_indices_structured
     
    
    
#nc_fid=xr.open_dataset('/net/hydro/chemie/mfriedel/Data/ozone_extremes/WACCM/INT_O3_2000/U.101-300.zm.nc')
#find_SSW_leap(nc_fid, 200, 'U', 'plev', True, 'WACCM')  
# nc_fid=Dataset('/Users/mfriedel/Documents/Data/WACCM/CLIM_O3_2000/B2000WCN.NOCOUPL.e122.f19_g16.0-199_O3_TUZ.nc')
# find_SSW(nc_fid, 200, 2000, 'U', 'plev')
    


