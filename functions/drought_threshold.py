# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:43:00 2016

@author: annaukkola
"""


#Calculates drought threshold


def drought_threshold(vec, perc=15, subset=float('nan'), monthly=False):
    
    import numpy as np

    #Subset vector if wanted (allows threshold to be calculated for a subset of total time period)
    if any(~np.isnan(subset)):
        vec=vec[subset]
            
    #Calculates threshold based on all data
    if monthly == False:

        threshold = np.nanpercentile(vec, perc) 
        
    #Calculates a separate threshold for each month (assumes input is full 12 month periods)
    elif monthly == True:

        threshold = np.zeros(12) * np.nan
    
        for k in range(len(threshold)):
            ind = list(range(k, len(vec), 12))  #Create sequence of indices for each month
            threshold[k] = np.nanpercentile(vec[ind], perc)    #Calculate threshold using values for each month
    
    
    return threshold;