#Finds drought end date using threshold_termination

def drought_end(vec, dry_day_ind, threshold):
    
    import numpy as np
    
    for k in range(dry_day_ind+1, len(vec)):
        
        #Check if time step exceeds threshold
        exceeds = vec[k] > threshold[k]
        
        if exceeds or k==len(vec): 
            return np.arange(dry_day_ind+1, k)
            
            