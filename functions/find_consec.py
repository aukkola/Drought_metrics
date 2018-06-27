# -*- coding: utf-8 -*-
"""
Created on Tue May 10 17:45:13 2016

@author: annaukkola
"""

### Finds consecutive indices (drought days/months) ###

### Returns 1 for TRUE and 0 for FALSE !!


def find_consec(drought_days):
    
    import numpy as np
    
    
    consec=np.zeros(len(drought_days))
    
    for k in range(len(consec)):
        
        #Set first value of consec to TRUE if first or second day of time series, 
        #cannot search back in that case
        if k==0:
            #if drought_days[k] < 1:
            consec[k] = False
            #else:
             #   consec[k] = True
                
        else:
            if drought_days[k]-drought_days[k-1] == 1:
                consec[k] = True
            else:
                consec[k] = False
            
            
    return consec;
