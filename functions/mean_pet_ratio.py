# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 19:08:55 2016

@author: annaukkola
"""


#Calculates percentile value of ratio (aet + mean_pet)/(pet + mean_pet)
#where mean_pet is annual mean PET, expressed here as monthly totals (mm/months)

def mean_pet_ratio(aet_vec, pet_vec, perc=15, subset=float('nan'), monthly=False):
    
    import numpy as np

    #Subset vector if wanted (allows threshold to be calculated for a subset of total time period)
    if ~np.isnan(subset):
        aet_vec = aet_vec[subset]
        pet_vec = pet_vec[subset]
        
        
        
    #Calculate mean annual PET (expressed in units mm/month)        
    mean_pet = sum(pet_vec)/len(pet_vec)        
    
    #Calculate ratio
    ratio = (aet_vec + mean_pet) / (pet_vec + mean_pet)            
       
        
    #Calculates threshold based on all data
    if monthly == False:

        #Determine percentile
        ratio_threshold = np.percentile(ratio, perc) 
        
        
    #Calculates a separate threshold for each month (assumes input is full 12 month periods)
    elif monthly == True:

        ratio_threshold = np.zeros(12) * np.nan
    
        for k in range(len(ratio_threshold)):
            
            ind = list(range(k, len(aet_vec), 12))  #Create sequence of indices for each month
             
            #Determine percentile
            ratio_threshold[k] = np.percentile(ratio[ind], perc)    #Calculate threshold using values for each month
    
    
    return {'ratio': ratio, 'threshold': ratio_threshold}