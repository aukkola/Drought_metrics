# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:54:23 2016

@author: annaukkola
"""


#Calculates the frequency, duration and intensity of droughts
#Default threshold is 15th percentile
#Code always returns duration, optional outputs include timing (index of dry days), magnitude and intensity

#Options
# perc sets percentile for calculating threshold
# obs_vec: can use obs/control run to calculate threshold
# subset: can calculate threshold based on a subset of data
# monthly: calculates a separate threshold for each month
# pet_lim: additional constraint for identifying drought days (AET/PET ratio below pet_ratio)


def drought_metrics(mod_vec, lib_path, obs_vec=[float('nan')], perc=15, scale=3,                           
                    subset=float('nan'),
                    monthly=False, return_all_tsteps=False,
                    pet_lim=False,  pet_ratio=0.2, pet_vec=float('nan'),    #use AET/PET ratio?
                    temp_lim=False, temp_val=10, temp_vec=float('nan'),     #use temperature limit? temp above a threshold
                    mean_pet_lim=False,                                     #use standardised (AET + mean(PET)) / (PET + mean(PET)) ratio?
                    add_metrics=(['timing', 'magnitude', 'intensity', 'threshold', 
                    'count_duration', 'count_magnitude', 'count_intensity']),
                    count=[ 1,  2,  3,  4,  5,  6]):  
    
    
    import numpy as np
    
    ## source extra functions ##
    import sys 
    import os
    sys.path.append(os.path.abspath(lib_path))
    from drought_threshold import drought_threshold
    from find_consec import find_consec
    from find_end import find_end
    
    if mean_pet_lim == True:
        from mean_pet_ratio import mean_pet_ratio


   # from remove_short_events import *


    #################################
    ### Calculate running mean ts ###
    #################################
    
    #There is an option to turn the time series to a running sum series
    #using "scale". This is similar to the scale used for SPI
    
    if scale > 1:
        
        #Calculate rolling sum time series (the last values of this are bogus)
        mod_rollsum = np.convolve(mod_vec, np.ones((scale,)))[(scale-1):] 
        
        #Create of vector with first values as NA, and take the rest as rolling
        #sum (ignoring the last values as they are not correct)
        mod_vec = np.append(np.zeros(scale-1) * np.nan , 
                        mod_rollsum[0:(len(mod_rollmean) - scale+1)])
    
    
        #Similarly for obs_ref if using
        if all(np.isnan) == False:
            
            obs_rollsum = np.convolve(obs_vec, np.ones((scale,)))[(scale-1):] 
            
            obs_vec = np.append(np.zeros(scale-1) * np.nan , 
                            obs_rollsum[0:(len(obs_rollmean) - scale+1)])
    
    
    
    ###########################
    ### Calculate threshold ###
    ###########################
    
    #If obs_vec not supplied, set to use mod_vec for threshold calculation instead
    
    if all(np.isnan(obs_vec)):
        vec = mod_vec
    else:
        vec = obs_vec
    
    
    #Calculate threshold value
    threshold = drought_threshold(vec, perc, subset, monthly)

    
    #Calculate monthly standardised AET/PET ratio if using this option
    #Returns a list of ratio and threshold
    if mean_pet_lim == True:
        mean_ratio = mean_pet_ratio(aet_vec=vec, pet_vec=pet_vec, perc=perc, subset=subset, monthly=monthly)
    
    #Calculate monthly AET/PET ratio
    elif pet_lim == True:
        pet_ratio_data = vec/pet_vec
        
        #correct zero division
        pet_ratio_data[pet_vec==0.] = vec[pet_vec==0.]
     
     
    
    ##############################
    ### Calculate monthly mean ###
    ##############################

    # It is used later to determine drought magnitude and intensity
    #Takes scale into account when determining average of each month

    #Initialise
    sum_vec = np.zeros(12) * np.nan

    #loop months
    for k in range(12):

        #Create sequence of vector indices for each month
        ind = list(range(k, len(mod_vec), 12))


        #Scale larger than 1, calculate mean from current and previous months
        if scale > 1:

            #Remove first element(s) of ind when month no. less than scale (to avoid negative indices)
            if k < (scale-1):
                ind = ind[1:len(ind)]


            #Calculate average from current month and counting back scale no. of months
            temp_sum  = np.zeros(len(ind))

            for i in range(len(ind)):
                temp_sum[i]  = np.sum(mod_vec[(ind[i] - scale + 1):(ind[i]+1)])

            #Take average
            sum_vec[k]  = np.mean(temp_sum)


        #Scale equals 1
        else:
            #Calculate mean using values for each month
            sum_vec[k] = np.mean(mod_vec[ind])


    #Repeat threshold vector for calculating additional metrics
    sum_vec = np.tile(sum_vec, len(mod_vec)/12)    #repeat number of months



    #########################
    ### Find drought days ###    
    #########################
    
    #Month/day indices where mod_vec below threshold    
    

    
    #No monthly
    if monthly == False:

        #No PET lim
        if pet_lim == False:
            dry_days = np.where(mod_vec < threshold)[0]  #need to add zero so return a vector
            
        #PET lim
        elif pet_lim == True and temp_lim == False:
            dry_days = np.where((mod_vec < threshold) & (pet_ratio_data < pet_ratio))[0]

        

    #Monthly
    elif monthly == True:
            
        dry_days = []
    
        #loop months
        for k in range(len(threshold)):   
                            
            #Create sequence of indices for each month
            ind = list(range(k, len(vec), 12))                       
            
            #No PET lim
            if pet_lim == False:        
                dry_ind = np.where(mod_vec[ind] < threshold[k])[0]  #Extract relevant month from mod_vec and find correct threshold
            
            #PET lim
            elif pet_lim == True and temp_lim == False:
                dry_ind = np.where( (mod_vec[ind] < threshold[k]) & (pet_ratio_data[ind] < pet_ratio) )[0]
            
            #PET lim and temp lim
            elif pet_lim == True and temp_lim == True:
                dry_ind = np.where( (mod_vec[ind] < threshold[k]) & (pet_ratio_data[ind] < pet_ratio) & (temp_vec[ind] > temp_val))[0]
                
            #Mean PET lim (aet below threshold and standardised aet/pet ratio below threshold)
            elif mean_pet_lim == True:      
                dry_ind = np.where( (mod_vec[ind] < threshold[k]) & (mean_ratio['ratio'][ind] < mean_ratio['threshold'][k]) )[0]
                
            
            #Find month indices corresponding to dry indices (have to convert ind to array to extract elements)
            dry_days.extend(np.array(ind)[dry_ind])                   
            
            
    #sort ascending
    dry_days = np.sort(dry_days)
        
  
  


    #####################
    ### Find duration ###    
    #####################

    #Filter out consecutive days to find 1st drought day

    if len(dry_days) > 0:
        
        consec = find_consec(dry_days)
        
        start = dry_days[np.where(consec == 0)[0]]
        
        end = find_end(consec=consec, start_days=start, dry_days=dry_days)
        
        #Calculate duration
        duration = end - start +1
        
    else:
        duration = 0


    ### Timing ###
    
    timing = np.nan
    
    if 'timing' in add_metrics:
        
        timing = np.zeros(len(mod_vec))
        
        if len(dry_days) > 0:
            timing[dry_days]=1
    
    
    
    #Count number of drought events of certain length (as defined in "count")
 
    count_duration = np.nan
 
    if 'count_duration' in add_metrics:
        
        count_duration = np.zeros(len(count))
        
        for c in count:
            count_duration[c-1] = len(np.where(duration==c)[0])
      
      
    
    #Repeat threshold vector for calculating additional metrics
    if monthly == True:
        threshold = np.tile(threshold, len(mod_vec)/12)   #repeat number of yrs
    elif monthly == False:
        threshold = np.repeat(threshold, len(mod_vec))    #repeat number of months



    ######################
    ### Find magnitude ###
    ######################
    
    magnitude = np.nan
    count_magnitude = np.nan

    if 'magnitude' in add_metrics:
        
        #If found dry days
        if len(dry_days) > 0:
            
             #initialise
             magnitude = np.zeros(len(start))
             
             for k in range(len(start)):
                 
                 #More than one consec day
                 if end[k] - start[k] > 0:
                     magnitude[k] = sum( threshold[start[k]:(end[k]+1)] - mod_vec[start[k]:(end[k]+1)] )  #Need to add 1 to end day because of python indexing!!
                 #One consec day only
                 else:
                     magnitude[k] = threshold[start[k]] - mod_vec[start[k]]
                
        #If no dry days       
        else:  
            magnitude = 0.
                
   

        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'count_magnitude' in add_metrics:
            
            count_magnitude = np.zeros(len(count))

            #If found dry days
            if len(dry_days) > 0:
                
                for c in count:
            
                    #find indices for events of duration 'c'
                    ind = np.where(duration==c)[0]                    
                    
                    if len(ind) > 0:                    
                        count_magnitude[c-1] = np.mean(magnitude[ind])

            
            

    ######################
    ### Find intensity ###
    ######################

    intensity = np.nan
    rel_intensity = np.nan
    count_intensity = np.nan

    if 'intensity' in add_metrics:
        
        
        #If found dry days
        if len(dry_days) > 0:

            #initialise
            intensity = np.zeros(len(start))        
        
            #loop through drought periods
            for k in range(len(start)):
                
                #More than one consec day
                if end[k] - start[k] > 0:
                    
                    intensity[k] = max( threshold[start[k]:(end[k]+1)] - mod_vec[start[k]:(end[k]+1)] )  #maximum difference
                
                #One consec day only
                else:
                    
                    intensity[k] = threshold[start[k]] - mod_vec[start[k]]
                    
        #If no dry days
        else:
                intensity = 0.
                
        
        
        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'count_intensity' in add_metrics:
        
            count_intensity = np.zeros(len(count))
            
            #If found dry days
            if len(dry_days) > 0:

                for c in count:
            
                    #find indices for events of duration 'c'
                    ind = np.where(duration==c)[0]  
                    
                    if len(ind) > 0:
                        count_intensity[c-1] = np.mean(intensity[ind])



    ### List outputs ###


    outs = {'duration': duration, 'timing': timing, 'magnitude': magnitude, 'intensity': intensity, 
            'rel_intensity': rel_intensity, 'threshold': threshold, 'count_duration': count_duration, 
            'count_magnitude': count_magnitude, 'count_intensity': count_intensity}


    return outs;
    
    
    
    