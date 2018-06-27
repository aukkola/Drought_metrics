# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:34:05 2016

@author: annaukkola

"""


from netCDF4 import Dataset # to work with NetCDF files
import numpy as np
import glob
import sys 
import os



### Set paths ###

#root_path = "/Users/annaukkola/Documents/CABLE_GWSP3/"
root_path = "/srv/ccrc/data45/z3509830/CMIP5/"
lib_path  = root_path + '/scripts/Python/functions'


sys.path.append(os.path.abspath(lib_path))
from drought_metrics import *

   
   
### Set variable ###
var_name="pr"
var_path="Precip"

#######################
### Set experiments ###
#######################


experiment=['amip', 'historical']  #CMIP5 experiments
#experiment=['historical']  #CMIP5 experiments

#################
### Set years ###
#################

start_yr= 1989
end_yr  = 2005

# year indices to extract correct data years
yr_ind = np.array(range((start_yr-start_yr)*12, (end_yr-start_yr)*12 + 12))  #Need to add one to end year, silly python 



############################
#### Set drought limits ####
############################

#Set percentile for drought threshold
perc=10

bins=np.arange(1,12+1)  #drought durations to count

pet_lim=False       #Use PET to determine ET droughts? NOT implemented in this example
temp_lim=False      #Use a temperature threshold to determine droughts? NOT implemented in this example

#Use observations for calculating threshold?
obs_ref = False
obs_var = 'ET_mean'   #ET_mean, ET_min, ET_max


#Use threshold determined separately for each month?
#If set to false, determines one threshold from all data.
#Set to false if using annual data
monthly=True



##################
### Load files ###
##################


for k in range(len(experiment)):


    #List all model names
    models = os.listdir(root_path + '/Processed_data_' + start_yr + '_' + end_yr + '/' + 
                        experiment[k] + "/" + var_name + "/")


    #Reference data for calculating threshold
    for m in range(len(models)):

        
        #If using obs as reference (set to ET data, fix later if needed...)
        if obs_ref:
            obs_file = (root_path + '/LandFlux/Regridded/' + experiment[k] + "/" + 
                        models[m] + '/LandFluxEVAL.89-05.monthly.diagnostic_mm_month_' + 
                        obs_var + "_" + models[m] + "_regrid.nc")
        
    
        ### Find CMIP5 files ###
        files = glob.glob(root_path + '/Processed_data_' + start_yr + '_' + end_yr + '/'  + 
                          experiment[k] + "/" + var_name + "/" + models[m] + "/*regrid.nc")
        
        
        ### Load data ###
        
        #Model data
        fh = Dataset(files[0], mode='r')
        all_data = fh.variables[var_name][yr_ind]
        data     = all_data.data
        mask     = all_data.mask
        
        #Get lon and lat (name varies by CMIP5 model)
        try:
            lat = fh.variables['latitude'][:]
            lon = fh.variables['longitude'][:]
        except:
            lat = fh.variables['northing'][:]
            lon = fh.variables['easting'][:]
     
       
        #fh.close()
        
        #Mask
        miss_val = -99999.0
        data[mask==True] = miss_val
    
        #Read reference data used to calculate threshold
        if obs_ref:
            obs_fh = Dataset(obs_file, mode='r')
            control_ref = obs_fh.variables[obs_var][:].data
    
            print "Using OBS ref"
    
        else:
            control_ref = data        
            
            
    
        #Read PET data if using option
        if pet_lim == True:
        
            petfile = glob.glob(root_path + '/Priestley_Taylor_PET/' + experiment[k] + 
                                "/" + models[m] +'/*regrid_LFmasked.nc')
                
            pet_fh = Dataset(petfile[0], mode='r')
            pet_data = pet_fh.variables['PET'][yr_ind].data
            pet_fh.close()

        #Read temperature data if using option
        if temp_lim == True:
            tempfile = glob.glob(root_path + '/Processed_data/' + experiment[k] + 
                                 "/tas/" + models[m] +'/*regrid_LFmasked.nc')
                
            temp_fh = Dataset(tempfile[0], mode='r')
            temp_data = temp_fh.variables['tas'][yr_ind].data
            temp_fh.close()


        
        #Calculate AET/PET ratio if using pet_lim
        if pet_lim:
            pet_ratio_data = data/pet_data
        
            pet_ratio_data[pet_data==0.] = 0.    
        
    
        #Flip latitude as reads map upside down
        #Should implement this if using obs, pet or temp limits... Fix later
        #data = data[:,::-1,:]
    
    
        ### Initialise output arrays ###
        duration     = np.zeros((len(data), len(lat), len(lon))) + miss_val # * np.nan
        magnitude    = np.zeros((len(data), len(lat), len(lon))) + miss_val # * np.nan
        intensity    = np.zeros((len(data), len(lat), len(lon))) + miss_val # * np.nan
        timing       = np.zeros((len(data), len(lat), len(lon))) + miss_val # * np.nan    
        
        count_duration  = np.zeros((len(bins), len(lat), len(lon))) + miss_val # * np.nan
        count_magnitude = np.zeros((len(bins), len(lat), len(lon))) + miss_val # * np.nan
        count_intensity = np.zeros((len(bins), len(lat), len(lon))) + miss_val # * np.nan

        
        if monthly:
            threshold    = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan
        else:
            threshold    = np.zeros((len(lat), len(lon))) + miss_val # * np.nan
    
    
        ### Calculate metrics ###
    
        #Print progress
        print 'Experiment', (k+1), '/', len(experiment), ', model: ', m+1, "/", len(models)
        
        #Loop through grid cells
        for i in range(len(lat)):
                    
            for j in range(len(lon)):
            
                 #Calculate metrics if cell not missing
                 if any(~mask[:,i,j]):  
                 
                 
                     #No PET lim
                     if pet_lim == False:
                         metric = drought_metrics(mod_vec=data[:,i,j], lib_path=lib_path, perc=perc, 
                                                  obs_vec=control_ref[:,i,j], monthly=monthly, count=bins)
                                             
                    
                     ### Write metrics to variables ###
                     duration[range(np.size(metric['duration'])),i,j]   = metric['duration']  #total drought duration (months)
    
                     magnitude[range(np.size(metric['magnitude'])),i,j] = metric['magnitude'] #average magnitude
                    
                     intensity[range(np.size(metric['intensity'])),i,j] = metric['intensity'] #average intensity
            
                     timing[range(np.size(metric['timing'])),i,j]       = metric['timing']    #drought timing (month index)
 

                     count_duration[range(np.size(metric['count_duration'])),i,j]   = metric['count_duration']  #total drought duration (months)

                     count_magnitude[range(np.size(metric['count_magnitude'])),i,j] = metric['count_magnitude']  #total drought duration (months)

                     count_intensity[range(np.size(metric['count_intensity'])),i,j] = metric['count_intensity']  #total drought duration (months)

   
                     if monthly:
                         threshold[:,i,j] = metric['threshold'][0:12]    #drought timing (month index)
                     else:
                         threshold[i,j]   = metric['threshold']
    
    
        ##############################
        ### Write result to NetCDF ###
        ##############################
        
        
        #Creat output path
        out_path = (root_path + '/Drought_metrics/' + experiment[k] + "/" + 
                    var_path + '/Obs_' + str(obs_ref) )
        
        #Add obs reference if using
        if obs_ref:
            out_path = out_path + "_" + obs_var
            
        
        #Add percentile
        out_path = out_path + '/Perc_' + str(perc) + '/'
    
        
        #Create output directory if doesn't exist
        if not os.path.exists(out_path):    
            os.makedirs(out_path)
    
        
        #Create output file name
        out_file = (out_path + '/' + models[m] + '_drought_metrics_perc_' + str(perc)
        
        #If using PET limit
        if pet_lim:
            out_file = out_file + '_pet_ratio_' + str(pet_ratio)

        #If using temperature limit
        if temp_lim:
            out_file = out_file + '_' + str(temp_val) 
            
        #Finally add year and experiment info
        out_file = out_file + "_" + experiment[k] + '_' + str(start_yr) + '_' + str(end_yr) + '.nc'
        
        
        
        # Open a new netCDF file for writing
        ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC") 
        
        # Create the output data
        # Create the x, y and time dimensions
        ncfile.createDimension('lat', lat.shape[0])
        ncfile.createDimension('lon', lon.shape[0])
        ncfile.createDimension('time', len(data))
        ncfile.createDimension('count', len(bins))
            
        if monthly:
            ncfile.createDimension('month', 12)
    
    
    
        # Create dimension variables
    
        longitude = ncfile.createVariable("lon",  'f8', ('lon',))
        latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
        time      = ncfile.createVariable("time", 'i4', ('time',))
        count     = ncfile.createVariable("count",'i4', ('count',))
    
        if monthly:
            month = ncfile.createVariable("month", 'i4', ('month',))
    
        #Create data variables
        data_dur  = ncfile.createVariable('duration', 'i4',('time','lat','lon'), fill_value=miss_val)
        data_mag  = ncfile.createVariable('magnitude','f8',('time','lat','lon'), fill_value=miss_val)
        data_int  = ncfile.createVariable('intensity','f8',('time','lat','lon'), fill_value=miss_val)
        data_tim  = ncfile.createVariable('timing',   'i4',('time','lat','lon'), fill_value=miss_val)
        data_cdur = ncfile.createVariable('count_duration',  'i4',('count','lat','lon'), fill_value=miss_val)
        data_cmag = ncfile.createVariable('count_magnitude', 'f8',('count','lat','lon'), fill_value=miss_val)
        data_cint = ncfile.createVariable('count_intensity', 'f8',('count','lat','lon'), fill_value=miss_val)
    
        #Create data variable for threshold
        if monthly:
            data_thr = ncfile.createVariable('threshold', 'f8',('month','lat','lon'), fill_value=miss_val)
        else:
            data_thr = ncfile.createVariable('threshold', 'f8',('lat','lon'), fill_value=miss_val)
    
    
        #Set variable attributes
        longitude.units = 'degrees_east'
        latitude.units  = 'degrees_north'
        
        data_dur.long_name= 'drought event duration (no. months)'
        data_mag.long_name= 'drought event magnitude'
        data_int.long_name= 'drought event intensity'
        data_tim.long_name= 'drought event timing (month index)'
        data_thr.long_name= 'drought threshold (mm)'
        
        data_cdur.long_name= 'number of drought events with length specified by count'
        data_cmag.long_name= 'mean magnituge of drought events of length count (mm)'
        data_cint.long_name= 'mean intensity of drought events of length count (mm)'
    
       
        if monthly:
            month[:] = range(1,12+1)
     
        # Write data to dimension variables
        longitude[:]=lon
        latitude[:] =lat
        time[:]     = range(1,len(data)+1)
        count[:]    = bins
               
        if monthly:
            month[:] = range(1,12+1)
     
        #Write data to data variables
        data_dur[:,:,:] = duration    
        data_mag[:,:,:] = magnitude
        data_int[:,:,:] = intensity
        data_tim[:,:,:] = timing
        
        data_cdur[:,:,:] = count_duration
        data_cmag[:,:,:] = count_magnitude
        data_cint[:,:,:] = count_intensity
    
    
        if monthly:    
            data_thr[:,:,:] = threshold
        else:
            data_thr[:,:] = threshold
        
        
        # Close the file
        ncfile.close()










