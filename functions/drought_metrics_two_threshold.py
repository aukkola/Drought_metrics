# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:54:23 2016

@author: annaukkola
"""


#Calculates the frequency, duration and intensity of droughts
#Default threshold for onset is 10th percentile and termination is 50th percentile
#Drought starts when a month/year falls below onset percentile and ends when
#the time series next exceeds drought termination threshold.
#Code always returns duration, optional outputs include timing (index of dry days), 
#magnitude and intensity. Droughts can also be classified by their duration using
#count option. In this case, drought characteristics are returned separately for each
#drought duration

#Options
# perc_onset sets percentile for calculating threshold for drought onset
# perc_termination sets the percentile for calculating drought termination threshold
# obs_vec: can use obs/control run to calculate threshold
# subset: can calculate threshold based on a subset of data. Useful if want to use
# a specific baseline period for thresholds
# monthly: calculates a separate threshold for each month (i.e. accounts for seasonality)
# return_all_tsteps: returns a vector the same length as mod_vec rather than just
# drought time steps
# add_metrics: additional metrics to be outputted

def drought_metrics_two_threshold(mod_vec, lib_path, obs_vec=[float('nan')],  
                                  perc_onset=10, perc_termination=50, scale=3,                           
                                  subset=float('nan'),
                                  monthly=False, return_all_tsteps=False,
                                  add_metrics=(['timing', 'magnitude', 'intensity', 
                                  'threshold_onset', 'threshold_termination',
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
    from drought_end import drought_end

    #import pdb


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
                        mod_rollsum[0:(len(mod_rollsum) - scale+1)])
    
    
        #Similarly for obs_ref if using
        if ~all(np.isnan(obs_vec)):
            
            obs_rollsum = np.convolve(obs_vec, np.ones((scale,)))[(scale-1):] 
            
            obs_vec = np.append(np.zeros(scale-1) * np.nan, 
                                obs_rollsum[0:(len(obs_rollsum) - scale+1)])
    
    
    
    ###########################
    ### Calculate threshold ###
    ###########################
    
    #If obs_vec not supplied, set to use mod_vec for threshold calculation instead
    
    if all(np.isnan(obs_vec)):
        vec = mod_vec
    else:
        vec = obs_vec
    
    
    #Calculate threshold value for drought onset and termination
    threshold_onset = drought_threshold(vec, perc_onset, subset, monthly)
    
    threshold_termination = drought_threshold(vec, perc_termination, subset, monthly)


    #Repeat termination threshold vector
    threshold_termination = np.tile(threshold_termination, int(len(mod_vec)/12))    #repeat number of months

    
    ##############################
    ### Calculate monthly mean ###
    ##############################

    # It is used later to determine drought magnitude and intensity
    #Takes scale into account when determining average of each month

    if monthly:
        
        #Initialise
        sum_vec = np.zeros(12) * np.nan

        #loop months
        for k in range(12):

            #Create sequence of vector indices for each month
            ind = list(range(k, len(mod_vec), 12))

            #Calculate mean using values for each month
            sum_vec[k] = np.nanmean(mod_vec[ind])

        #Repeat mean vector for calculating additional metrics
        sum_vec = np.tile(sum_vec, int(len(mod_vec)/12))    #repeat number of months


    #Not monthly
    else:
        
        #Calculate mean of whole time series and repeat for the same length
        sum_vec = np.tile(np.nanmean(mod_vec), len(mod_vec))



    #########################
    ### Find drought days ###    
    #########################
    
    #First find months/days when mod_vec is below onset percentile
        
    #Not monthly
    if monthly == False:

        dry_days = np.where(mod_vec < threshold_onset)[0]  #need to add zero so return a vector
            
        
    #Monthly
    elif monthly == True:
            
        dry_days = []
    
        #loop months
        for k in range(len(threshold_onset)):   
                            
            #Create sequence of indices for each month
            ind = list(range(k, len(mod_vec), 12))                       
            
            #Extract relevant month from mod_vec and find correct threshold
            dry_ind = np.where(mod_vec[ind] < threshold_onset[k])[0]  
            
            #Find month indices corresponding to dry indices (have to convert ind to array to extract elements)
            dry_days.extend(np.array(ind)[dry_ind])                   
    
            
    #sort ascending
    dry_days = np.sort(dry_days)
        
  
    #Then loop through dry days to find drought ending days
    #(i.e. the first instance after drought onset when mod_vec goes above 
    #drought termination threshold)
    for k in range(len(dry_days)): 
        
        dry_days=np.append(dry_days, drought_end(vec=mod_vec, dry_day_ind=dry_days[k], 
                                                 threshold=threshold_termination))
        
    #sort ascending and remove duplicates
    dry_days = np.unique(np.sort(dry_days))
    


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
        duration = np.nan


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
      
      
    
    # #Repeat threshold vector for calculating additional metrics
    # if monthly == True:
    #     threshold_onset = np.tile(threshold, int(len(mod_vec)/12))   #repeat number of yrs
    # elif monthly == False:
    #     threshold = np.repeat(threshold, int(len(mod_vec)))    #repeat number of months



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
                     magnitude[k] = sum( sum_vec[start[k]:(end[k]+1)] - mod_vec[start[k]:(end[k]+1)] )  #Need to add 1 to end day because of python indexing!!
                 #One consec day only
                 else:
                     magnitude[k] = sum_vec[start[k]] - mod_vec[start[k]]
                
            
   

        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'count_magnitude' in add_metrics:
            
            count_magnitude = np.zeros(len(count)) * np.nan

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

    #Absolute (mm) and relative intensity (% departure from mean conditions)
    if 'intensity' or 'rel_intensity' in add_metrics:
                
        #If found dry days
        if len(dry_days) > 0:

            #initialise
            intensity = np.zeros(len(start))        
            rel_intensity = np.zeros(len(start))

            #Loop through drought periods
            for k in range(len(start)):


                #More than one consec day
                if end[k] - start[k] > 0:

                    #Initialise temporary intensity vector to match drought length
                    temp_int     = np.zeros(end[k]-start[k]+1)

                    #Create a vector of indices for current drought event
                    ind = np.arange(start[k], end[k]+1) #Need to add 1 to end day because of python indexing!!

                    #Calculate intensity for each month, taking into account scale
                    for d in range(len(ind)):                    
                    
                        #Absolute intensity (mean - drought month value)
                        temp_int[d]     = sum_vec[ind[d]] - mod_vec[ind[d]] 


                    #Average monthly magnitudes to get event intensity
                    intensity[k] = np.mean(temp_int)

                    #Relative intensity abs( (m - mean) / mean * 100)), where m is drought month value
                    #(using simplified version of this from Ned)
                    rel_intensity[k] = abs(( np.mean(mod_vec[ind]) / np.mean(sum_vec[ind]) -1)) * 100



                #One consec day only
                else:

                    #Absolute intensity (mean - drought month value)
                    intensity[k] = sum_vec[start[k]] - mod_vec[start[k]]

                    #Relative intensity abs( (m - mean) / mean * 100)), where m is drought month value
                    rel_intensity[k] = abs( (mod_vec[start[k]] - sum_vec[start[k]]) /
                                       sum_vec[start[k]] * 100 ) #Need to add 1 to end day because of python indexing!!
            
           
        
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


    #If returning metrics as time series, modify metrics here
    #(do not apply this to "count variables", perhaps should be changed later)
    if return_all_tsteps:
        
        if len(dry_days) > 0:
            
            #I'm sure there is a smarter way to do this?
            
            temp_dur     = np.zeros(len(mod_vec)) * np.nan
            temp_mag     = np.zeros(len(mod_vec)) * np.nan
            temp_int     = np.zeros(len(mod_vec)) * np.nan
            temp_rel_int = np.zeros(len(mod_vec)) * np.nan

            temp_dur[start]     = duration
            temp_mag[start]     = magnitude
            temp_int[start]     = intensity
            temp_rel_int[start] = rel_intensity

            duration=temp_dur
            magnitude=temp_mag
            intensity=temp_int
            rel_intensity=temp_rel_int
    


    #Compile outputs
    outs = {'duration': duration, 'timing': timing, 'magnitude': magnitude, 'intensity': intensity, 
            'rel_intensity': rel_intensity, 'threshold_onset': threshold_onset, 
            'threshold_termination': threshold_termination,
            'count_duration': count_duration, 
            'count_magnitude': count_magnitude, 'count_intensity': count_intensity,
            'tseries': mod_vec}


    return outs;
    
    
    
    