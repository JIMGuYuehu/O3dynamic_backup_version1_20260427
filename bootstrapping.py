from netCDF4 import Dataset
import numpy as np
import random


def bootstrapping(array, indices, composite, years):
    
    #This function calculates the significance of 2D anomalies (lon-lat) in a 30-day window after the events compared to the climatology
    #based on a 500 sample bootstrapping 
    
    #note that this function only works for data that don't include leap years (i.e. WACCM). For data that include leap years, please use the function "bootstrapping_leap" (note that that function is less efficient)
    
    #"array": variable for which bootstrapping should be calculated
    #"indices": indices of the event of interest (SSW, ozone minima, ...)
    #"years": number of years in the dataset
    
    bootstrap = 500
    bootstrap_composite = np.zeros((bootstrap, len(array[0, 0, :, 0]), len(array[0, 0, 0, :])))
    
    print('starting bootstrapping...')
    
    i = 0
    while i < bootstrap:
        var_all = np.zeros((len(indices) * 30, len(array[0, 0, :, 0]), len(array[0, 0, 0, :])))
        for index in range(len(indices)):
            indice = int(indices[index]) + 59 #because index indicates day after March 1st, convert to day of the year
            random_year = random.randint(0, years - 1)
            var_all[index * 30:index * 30 + 30, :, :] = array[random_year, indice:indice + 30, :, :]

        bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
        i = i + 1

    print('bootstrapping completed.')
    
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    
    return significance


#_________________________________________________________________________________________________________________
    


def bootstrapping_leap(nc_fid, array, indices, composite, O3_years):
    
    
    #This function calculates the significance of 2D anomalies (lon-lat) in a 30-day window after the events compared to the climatology
    #based on a 500 sample bootstrapping for datasets with or without leap years 
    
    #"nc_fid": dataset containing variable of interest
    #"array": variable for which bootstrapping should be calculated
    #"indices": indices of the event of interest (SSW, ozone minima, ...)
    #"years": number of years in the dataset
    
    
    bootstrap = 500
    bootstrap_composite = np.zeros((bootstrap, len(array[0, :, 0]), len(array[0, 0, :])))
    
    print('starting bootstrapping...')
    i = 0
    while i < bootstrap:
        var_all = np.zeros((len(indices) * 30, len(array[0, :, 0]), len(array[0, 0, :])))
        for index in range(len(indices)):
            indice = int(indices[index]) + 59 #because index indicates day after March 1st, convert to day of the year
          #  print(indice)
            random_number = random.randint(0, len(O3_years) - 1)
            random_year = O3_years[random_number]
            array_random_year = array.sel(time=(nc_fid.time.dt.year.isin([random_year])))
            var_all[index * 30:index * 30 + 30, :, :] = array_random_year[indice:indice + 30, :, :]

        bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
        i = i + 1

    print('bootstrapping finished.')
    
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    return significance


#_________________________________________________________________________________________________________________
    



def bootstrapping_leap_alt_time(nc_fid, array, indices, composite, O3_years, model):
    
    
    #This function calculates the significance of 2D anomalies (height-time) in a -30/+60-day window around the events compared to the climatology
    #based on a 500 sample bootstrapping for datasets with or without leap years 
    
    #"nc_fid": dataset containing variable of interest
    #"array": variable for which bootstrapping should be calculated
    #"indices": indices of the event of interest (SSW, ozone minima, ...)
    #"composite": mean composite of anomalies of the variable of interest in the -30/+60-days window around the ozone minima
    #"O3_years": ozone minima years
    #"model": SOCOL, WACCM or MERRA
    
 
    bootstrap = 500
    if model == 'WACCM':
        leap = False
    else:
        leap = True
        
    if leap == True:
        bootstrap_composite = np.zeros((bootstrap, 90, len(array[0, :])))
        i = 0
        print('starting bootstrapping....')
        while i < bootstrap:
            var_all = np.zeros((len(indices), 90, len(array[0, :])))
            for index in range(len(indices)):
                indice = indices[index] + 59 #because index indicates day after March 1st, convert to day of the year
                random_number = random.randint(0, len(O3_years) - 2)
                random_year = O3_years[random_number]
                array_random_year = array.sel(time=(nc_fid.time.dt.year.isin([random_year])))
                if len(array_random_year.dims) == 3:
                    var_all[index, :, :] = array_random_year[indice - 29:indice + 61, :, 0]
                else:
                    var_all[index, :, :] = array_random_year[indice - 29:indice + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
            i = i + 1

    if leap == False:
        bootstrap_composite = np.zeros((bootstrap, 90, len(array[0, 0, :])))
        i = 0
        print('starting bootstrapping....')
        while i < bootstrap:
            var_all = np.zeros((len(indices), 90, len(array[0, 0, :])))
            for index in range(len(indices)):
                indice = indices[index] + 59 #because index indicates day after March 1st, convert to day of the year

                random_number = random.randint(0, len(O3_years) - 1)
                var_all[index, :, :] = array[random_number, indice - 29:indice + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
            i = i + 1

    print('bootstrapping succesfull.')
    
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    
    return significance


#_________________________________________________________________________________________________________________
  


def bootstrapping_space_2samp_1D(nc_fid1, nc_fid2, array1, array2, indices1, indices2, diff, O3_years1, O3_years2, model):
    
    #This function calculates whether the difference of two AO composites is statistically different.
    # First, a distribution of the AO difference in two different datasets by selecting the date of occurence of ozone minima
    #in random years in both datasets before taking the difference. The composite difference is then considered
    # statistically significant if it differs more than two standard deviations from the mean of the distribution.
    
    #"nc_fid1"/"nc_fid2": datasets containing variables of interest
    #"array1"/"array2": variables for which bootstrapping should be calculated
    #"indices1"/"indices2": indices of the event of interest in dataset 1 and 2
    #"diff": mean composite difference of anomalies 
    #"O3_years1"/"O3_years2": ozone minima years in dataset 1 and 2
    #"model": SOCOL, WACCM or MERRA
    
    
    bootstrap = 5000
    if model == 'WACCM':
        leap = False
    else:
        leap = True
    print('starting bootstrapping ....')
    
    if leap == True: # if the dataset contains leap years, set this to TRUE. Note that this is very inefficient.
        bootstrap_composite = np.zeros(bootstrap)
        i = 0
        
        while i < bootstrap:

            var_all1 = np.zeros(len(indices1) * 30)
            var_all2 = np.zeros(len(indices2) * 30)
            
            for index in range(len(indices1)):
                indice1 = indices1[index] + 59 #because index indicates day after March 1st, convert to day of the year
                indice2 = indices2[index] + 59 #because index indicates day after March 1st, convert to day of the year
                
                random_number1 = random.randint(0, len(O3_years1) - 1)
                random_year1 = O3_years1[random_number1]
                random_number2 = random.randint(0, len(O3_years2) - 1)
                random_year2 = O3_years2[random_number2]
                
                array_random_year1 = array1.sel(time=(nc_fid1.time.dt.year.isin([random_year1])))
                array_random_year2 = array2.sel(time=(nc_fid2.time.dt.year.isin([random_year2])))
                
                var_all1[index * 30:index * 30 + 30] = array_random_year1[indice1:indice1 + 30]
                var_all2[index * 30:index * 30 + 30] = array_random_year2[indice2:indice2 + 30]

            bootstrap_composite[i] = np.mean(var_all1) - np.mean(var_all2)
            i = i + 1

    if leap == False:
        bootstrap_composite = np.zeros(bootstrap)
        i = 0
        
        while i < bootstrap:
            var_all1 = np.zeros(len(indices1) * 30)
            var_all2 = np.zeros(len(indices2) * 30)
            
            for index in range(len(indices1)):
                
                indice1 = indices1[index] + 59 #because index indicates day after March 1st, convert to day of the year
                indice2 = indices2[index] + 59 #because index indicates day after March 1st, convert to day of the year
                
                random_number1 = random.randint(0, len(O3_years1) - 1)
                random_number2 = random.randint(0, len(O3_years2) - 1)
                
                var_all1[index * 30:index * 30 + 30] = array1[random_number1, indice1:indice1 + 30]
                var_all2[index * 30:index * 30 + 30] = array2[random_number2, indice2:indice2 + 30]

            bootstrap_composite[i] = np.mean(var_all1) - np.mean(var_all2)
            i = i + 1

    print('bootstrapping done.')
  
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - diff)
    ratio = diff_bootstrap / std_bootstrap
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    
    
    return (significance, ratio, bootstrap_composite)


#_________________________________________________________________________________________________________________
  

def bootstrapping_2samp(nc_fid_1, nc_fid_2, array1, array2, diff, indices1, indices2, O3_years1, O3_years2, model):
    
        #This function calculates whether the difference of two time-lat anomalies composites is statistically different.
    # First, a distribution of random composite differences between the two datasets is calculated by selecting the date of occurence of ozone minima
    #in random years in both datasets before taking the difference. The composite difference is then considered
    # statistically significant if it differs more than two standard deviations from the mean of the distribution.
    
    #"nc_fid1"/"nc_fid2": datasets containing variables of interest
    #"array1"/"array2": variables for which bootstrapping should be calculated
    #"indices1": indices of the events of interest in dataset 1
    #"indices2": indices of the events of interest in dataset 2
    #"O3_years1": ozone minima years in dataset 1
    #"O3_years2": ozone minima years in dataset 2
    #"model": SOCOL, WACCM or MERRA
    
    bootstrap = 500
    if model == 'WACCM':
        leap = False
    if model == 'SOCOL':
        leap = True
        
    print('starting bootstrapping....')
    
    if leap == True:
        bootstrap_composite = np.zeros((bootstrap, 90, len(array1[0, :])))
        i = 0
        while i < bootstrap:

            var_all_1 = np.zeros((len(indices1), 90, len(array1[0, :])))
            var_all_2 = np.zeros((len(indices1), 90, len(array2[0, :])))
            
            for index in range(len(indices1)):
                
                indice_1 = indices1[index] + 59 #because index indicates day after March 1st, convert to day of the year
                indice_2 = indices2[index] + 59 #because index indicates day after March 1st, convert to day of the year
                random_number_1 = random.randint(0, len(O3_years1) - 1)
                random_year_1 = O3_years1[random_number_1]
                array_random_year_1 = array1.sel(time=(nc_fid_1.time.dt.year.isin([random_year_1])))
                
                if len(array_random_year_1.dims) == 3:
                    var_all_1[index, :, :] = array_random_year_1[indice_1 - 29:indice_1 + 61, :, 0]
                else:
                    var_all_1[index, :, :] = array_random_year_1[indice_1 - 29:indice_1 + 61, :]
                    
                random_number_2 = random.randint(0, len(O3_years2) - 1)
                random_year_2 = O3_years2[random_number_2]
                array_random_year_2 = array2.sel(time=(nc_fid_2.time.dt.year.isin([random_year_2])))
                
                if len(array_random_year_2.dims) == 3:
                    var_all_2[index, :, :] = array_random_year_2[indice_2 - 29:indice_2 + 61, :, 0]
                else:
                    var_all_2[index, :, :] = array_random_year_2[indice_2 - 29:indice_2 + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all_1, axis=0) - np.mean(var_all_2, axis=0)
            i = i + 1

    if leap == False:
        
        bootstrap_composite = np.zeros((bootstrap, 90, len(array1[0, 0, :])))
        i = 0
        
        while i < bootstrap:
            var_all_1 = np.zeros((len(indices1), 90, len(array1[0, 0, :])))
            var_all_2 = np.zeros((len(indices1), 90, len(array2[0, 0, :])))
            
            for index in range(len(indices1)):
                indice1 = indices1[index] + 59 #because index indicates day after March 1st, convert to day of the year
                indice2 = indices2[index] + 59 #because index indicates day after March 1st, convert to day of the year
                
                random_number1 = random.randint(0, len(O3_years1) - 1)
                random_year1 = O3_years1[random_number1]
                random_number2 = random.randint(0, len(O3_years2) - 1)
                random_year2 = O3_years2[random_number2]
                
                var_all_1[index, :, :] = array1[random_number1, indice1 - 29:indice1 + 61, :]
                var_all_2[index, :, :] = array2[random_number2, indice2 - 29:indice2 + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all_1, axis=0) - np.mean(var_all_2, axis=0)
            i = i + 1

    print('boostrapping done...')
    
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - diff)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    
    return significance



#_________________________________________________________________________________________________________________


def bootstrapping_space_2samp(nc_fid1, nc_fid2, array1, array2, indices1, indices2, diff, O3_years1, O3_years2, model):
    
        
        #This function calculates whether the difference of two lat-lon anomalies composites is statistically different.
    # First, a distribution of random composite differences between the two datasets is calculated by selecting the date of occurence of ozone minima
    #in random years in both datasets before taking the difference. The composite difference is then considered
    # statistically significant if it differs more than two standard deviations from the mean of the distribution.
    
    #"nc_fid1"/"nc_fid2": datasets containing variables of interest
    #"array1"/"array2": variables for which bootstrapping should be calculated
    #"indices1": indices of the events of interest in dataset 1
    #"indices2": indices of the events of interest in dataset 2
    #"diff": difference of the two composites
    #"O3_years1": ozone minima years in dataset 1
    #"O3_years2": ozone minima years in dataset 2
    #"model": SOCOL, WACCM or MERRA
    
    bootstrap = 500
    if model == 'WACCM':
        leap = False
    else:
        leap = True
    print('starting bootstrapping ....')
    if leap == True:
        bootstrap_composite = np.zeros((bootstrap, len(array1[0, :, 0]), len(array1[0, 0, :])))
        i = 0
        while i < bootstrap:
            var_all1 = np.zeros((len(indices1) * 30, len(array1[0, :, 0]), len(array1[0, 0, :])))
            var_all2 = np.zeros((len(indices2) * 30, len(array2[0, :, 0]), len(array2[0, 0, :])))
            
            for index in range(len(indices1)):
                indice1 = indices1[index] + 59
                indice2 = indices2[index] + 59
                
                random_number1 = random.randint(0, len(O3_years1) - 4)
                random_year1 = O3_years1[random_number1]
                random_number2 = random.randint(0, len(O3_years2) - 4)
                random_year2 = O3_years2[random_number2]
                
                array_random_year1 = array1.sel(time=(nc_fid1.time.dt.year.isin([random_year1])))
                array_random_year2 = array2.sel(time=(nc_fid2.time.dt.year.isin([random_year2])))
                var_all1[index * 30:index * 30 + 30, :, :] = array_random_year1[indice1:indice1 + 30, :, :]
                var_all2[index * 30:index * 30 + 30, :, :] = array_random_year2[indice2:indice2 + 30, :, :]

            bootstrap_composite[i, :, :] = np.mean(var_all1, axis=0) - np.mean(var_all2, axis=0)
            i = i + 1

    if leap == False:
        bootstrap_composite = np.zeros((bootstrap, len(array1[0, 0, :, 0]), len(array1[0, 0, 0, :])))
        i = 0
        while i < bootstrap:
            var_all1 = np.zeros((len(indices1) * 30, len(array1[0, 0, :, 0]), len(array1[0, 0, 0, :])))
            var_all2 = np.zeros((len(indices2) * 30, len(array2[0, 0, :, 0]), len(array2[0, 0, 0, :])))
            
            for index in range(len(indices1)):
                indice1 = indices1[index] + 59
                indice2 = indices2[index] + 59
                
                random_number1 = random.randint(0, len(O3_years1) - 1)
                random_year1 = O3_years1[random_number1]
                random_number2 = random.randint(0, len(O3_years2) - 1)
                random_year2 = O3_years2[random_number2]
                
                var_all1[index * 30:index * 30 + 30, :, :] = array1[random_number1, indice1:indice1 + 30, :, :]
                var_all2[index * 30:index * 30 + 30, :, :] = array2[random_number2, indice2:indice2 + 30, :, :]

            bootstrap_composite[i, :, :] = np.mean(var_all1, axis=0) - np.mean(var_all2, axis=0)
            i = i + 1

    print('bootstrapping done.')
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - diff)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    return significance



#_________________________________________________________________________________________________________________


def bootstrapping_leap_alt_time(nc_fid, array, indices, composite, O3_years, model):
    
    bootstrap = 500
    if model == 'WACCM':
        leap = False
    else:
        leap = True
    if leap == True:
        
        bootstrap_composite = np.zeros((bootstrap, 90, len(array[0, :])))
        i = 0
        print('starting bootstrapping....')
        while i < bootstrap:
            var_all = np.zeros((len(indices), 90, len(array[0, :])))
            for index in range(len(indices)):
                indice = indices[index] + 59
                random_number = random.randint(0, len(O3_years) - 2)
                random_year = O3_years[random_number]
                array_random_year = array.sel(time=(nc_fid.time.dt.year.isin([random_year])))
                if len(array_random_year.dims) == 3:
                    var_all[index, :, :] = array_random_year[indice - 29:indice + 61, :, 0]
                else:
                    var_all[index, :, :] = array_random_year[indice - 29:indice + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
            i = i + 1

    if leap == False:
        bootstrap_composite = np.zeros((bootstrap, 90, len(array[0, 0, :])))
        i = 0
        print('starting bootstrapping....')
        while i < bootstrap:
            var_all = np.zeros((len(indices), 90, len(array[0, 0, :])))
            for index in range(len(indices)):
                indice = indices[index] + 59
                random_number = random.randint(0, len(O3_years) - 1)
                var_all[index, :, :] = array[random_number, indice - 29:indice + 61, :]

            bootstrap_composite[i, :, :] = np.mean(var_all, axis=0)
            i = i + 1

    print('bootstrapping succesfull.')
    
    mean_bootstrap = np.mean(bootstrap_composite, axis=0)
    std_bootstrap = np.std(bootstrap_composite, axis=0)
    diff_bootstrap = np.abs(mean_bootstrap - composite)
    significance = np.greater(np.abs(diff_bootstrap), 2 * std_bootstrap)
    
    return significance