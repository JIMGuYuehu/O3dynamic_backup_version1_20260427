from netCDF4 import Dataset, MFDataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import xarray as xr



# THIS PROGRAM CALCULATES THE FINAL WARMING DATE BASED ON DAILY VALUES
# (FINAL WIND REVERSAL AT 10 hPA AND 60°N)



def find_nearest(array, value):
     
    index =np.argmin(np.abs(array-value))
    
    return(index)


#___________________________________________________________________________________________________
    
    

def find_FW_vertical(nc_fid, years, var, lev, Data):

    #this function calculates the FSW date at different levels
    
    #"nc_fid": file containting wind variable
    #"years": number of years in the file
    #"var": name of wind variable
    #"lev": name of wind variable
    #"Data": name of model: MERRA, WACCM, or SOCOL
    
    plev = nc_fid[lev]
    
    if Data=='MERRA':
        plev=plev.sel(lev=slice(50,0.1))
    else:    
        if plev[len(plev)-1] > 1000:
            plev=plev.sel(plev=slice(10,5000))
            plev=plev/100
        else:
           plev=plev.sel(plev=slice(0.1,50))

        
    FW_dates_lev=np.zeros((len(plev),years))
    
    for i in range(len(plev)):
        FW_dates_lev[i,:] = np.reshape(np.array(find_FW_new_leap(nc_fid, years, var, lev, Data, plev[i])), (years)) 
        
  
    #np.savetxt('FW_dates_vertical_MERRA_thresh7.txt', FW_dates_lev) #save FSW dates as text file
    
    fig = plt.figure(figsize=(10,8))     
        
    plt.plot(np.nanmean(FW_dates_lev, axis=1),np.log(plev), color='k')
    plt.plot(np.nanmean(FW_dates_lev, axis=1)-np.std(FW_dates_lev, axis=1),np.log(plev),linestyle='--',color='k')
    plt.plot(np.nanmean(FW_dates_lev, axis=1)+np.std(FW_dates_lev, axis=1),np.log(plev),linestyle='--',color='k')
    plt.gca().invert_yaxis()
    plt.yticks([np.log(100),np.log(10),np.log(1)],('100','10','1'), fontsize=18)
    plt.xticks([0,30,59,90,120,151], ('January','February','March', 'Apr', 'May','June'), fontsize=18) 
    
   
    return FW_dates_lev
  
    
#___________________________________________________________________________________________________
    

    
def find_FW_new_leap(nc_fid, years, var, lev, Data, level):
    
    # This function uses the IMPROVED FW definition of Gerber & Butler, 2018:
    # The end of the vortex, or final warming (FW), occurs on the FIRST date when the winds reverse 
    # and do not return to westerly for more than 10 consecutive days. 
    
    #"nc_fid": file containing wind variable
    # "years": number of years in the file
    #"var": wind variable name
    #"lev": name of the pressure level variable
    #"level": level at which the FSW should be calculated
    
    plot=False

    plev = nc_fid[lev]

    # get wind variable
        
    wind_u = nc_fid[var] 

    wind_u=wind_u.interp(lat=np.linspace(-90,90,73))

    #extract wind at 10hPa and 60°N --> LINEAR INTERPOLATION if 60 °N not in lats
    
    if Data=='MERRA':
        wind_u=wind_u.sel(lat=60,lev=level, method='nearest') 
    else:    
        if plev[len(plev)-1] > 1000:
            wind_u=wind_u.sel(lat=60,plev=level*100, method='nearest')
        else:
            wind_u=wind_u.sel(lat=60,plev=level, method='nearest') 


    wind_u_FW=wind_u.sel(time=nc_fid.time.dt.month.isin([1,2,3,4,5,6])) #select months January - June

    wind_u_FW=wind_u_FW.groupby("time.year")    
    
    FW_dates=[]    
    
    day_start=0

    if level<=10: thresh = 0    #for pressure levels <=10 hPa: use wind threshold of 0 m/s
    else: thresh = 7            #for pressure levels >10 hPa: use adjusted wind threshold of 7 m/s to account for cold pole bias in models 

    for year, group in wind_u_FW: #loop through years
        
        day_max=len(np.array(group)) #how many days are in this year?
 
        for day in range(day_start,day_max): #check for every day of the year if the wind is still stronger than the threshold
            
            if group[day]<thresh:
                
                count_westerly_days=0 
                
                for i in range(day+1, day_max-10): #check if wind turns westerly for not more than 10 consecutive days until the end of June
                    
                    if group[i]>thresh:
                        
                        if all(j > thresh for j in group[i:i+10]):
                            count_westerly_days=11 
                            
                    else: continue 
               
                if count_westerly_days<=10: # if wind does not turn westerly for more than 10 consecutive days: FSW date is found
                    FW_dates.append(day) 
            
                    if Data == 'none':
                        day=day+1
                        if (year-1980)%4==0:
                
                            if day > 60 and day <= 91 :
                                FW_day = day - 31 - 29
                                FW_month = 3
                            if day > 91 and day <=121:
                                FW_day = day -31 - 29 -31
                                FW_month = 4
                            if day > 121:
                                FW_day = day - 31 -29 - 31 - 30
                                FW_month = 5
                                
                        if (year-1980)%4!=0:
                            if day > 60 and day <= 90 :
                                FW_day = day - 31 - 28
                                FW_month = 3
                            if day > 90 and day <=121:
                                FW_day = day -31 - 28 -31
                                FW_month = 4
                            if day > 121:
                                FW_day = day - 31 -28 - 31 - 30
                                FW_month = 5
                        print(str(FW_month) + '        ' + str(FW_day) + '      ' + str(year) + '      ' + str(count_westerly_days))
                    
                    break
        
            if day==day_max-1: FW_dates.append(np.nan)
   
    if plot==True:
        fig = plt.figure(figsize=(10,8)) 
        plt.hist(FW_dates, bins=20, color='grey',density=True)
        plt.ylabel('density', fontsize=18)
        plt.xticks([0,31,59,90,120,151,182], ('January','February','March', 'Apr', 'May','June'), fontsize=18) 
          
        x=np.linspace(np.nanmean(FW_dates)-50,np.nanmean(FW_dates)+50,100)
        plt.plot(x, norm.pdf(x, np.nanmean(FW_dates), np.nanstd(FW_dates)), color='k' , label='pdf')
        plt.title('Final warming date', fontsize=18)
        
        print('Mean final warming day (day of the year): ' + str(np.nanmean(FW_dates)) + " ± " + str(np.nanstd(FW_dates)))
    
    return FW_dates       
    

