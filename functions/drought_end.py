#Finds drought end date using threshold_termination

def drought_end(vec, dry_day_ind, threshold):
    
    import numpy as np
    
    for k in range(dry_day_ind+1, len(vec)): #need len(vec) here because stupid python ignores last number of range
        
        #Check if time step exceeds threshold
        exceeds = vec[k] > threshold[k]
        
        if exceeds or k==(len(vec)-1): #need len(vec)-1 here 
            return np.arange(dry_day_ind+1, k)
            
            