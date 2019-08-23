# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:54:23 2016

@author: annaukkola
"""


#Calculates the frequency, duration and intensity of droughts
#Default threshold is 15th percentile
#Code always returns duration, optional outputs include timing (index of dry days), magnitude and intensity
# and the magnitude/intensity of droughts of lengths specified in 'count'


#Options
# mod_vec:      input data vector
# severity:     threshold for determing drought months (one of 'moderate', 'severe' or 'extreme' (see code below for more details))
# scale:        scale for calculating SPI/SPEI
# subset:       can calculate threshold based on a subset of data
# monthly:      calculates a separate threshold for each month
# add_metrics:  additional metrics to be returned by code (timing: index of drought months)
# count:        drought lengths for which to return frequency, magnitude and intensity


def drought_metrics_SPI(mod_vec, spi_vec, lib_path, severity, scale,
                        add_metrics=(['timing', 'magnitude', 'intensity', 'rel_intensity',
                        'count_duration', 'count_magnitude', 'count_intensity', 'count_rel_intensity',
                        'SPI_magnitude', 'SPI_count_magnitude',
                        'SPI_intensity', 'SPI_count_intensity']),
                        count=[ 1,  2,  3,  4,  5,  6], miss_val=float('nan'),
                        use_min_threshold=False):


    # source packages
    import sys
    import os
    import numpy as np

    #Set paths and import drought functions
    sys.path.append(os.path.abspath(lib_path))
    from find_consec import find_consec
    from find_end import find_end



    #########################
    ### Find drought days ###
    #########################


    # Set SPI threshold values according to selected severity
    # Criteria from http://www.wamis.org/agm/pubs/SPI/WMO_1090_EN.pdf

    if severity == 'moderate':
        threshold = -1.0
        min_threshold = -1.5
    elif severity == 'severe':
        threshold = -1.5
        min_threshold = -2.0
    elif severity == 'extreme':
        threshold = -2.0
        min_threshold = -100.0 #picked a random no.
    else:
        sys.exit("Check drought severity input, not defined correctly !!")



    #Suppress warning message for below command (which produces a warning when NaNs present in spi_vec)
    #A bit dangerous, should probably fix later
    np.seterr(invalid='ignore')


    ### Find dry months ###
    if use_min_threshold==False:
    	dry_days = np.where(spi_vec <= threshold)[0]
    else:
	dry_days = np.where((spi_vec <= threshold) & (spi_vec > min_threshold))[0]

    #sort ascending
    dry_days = np.sort(dry_days)



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


    #####################
    ### Find duration ###
    #####################

    #Filter out consecutive days to find 1st drought day

    if len(dry_days) > 0:

        consec = find_consec(dry_days)

        start = dry_days[np.where(consec == 0)[0]]

        end = find_end(consec=consec, start_days=start, dry_days=dry_days)

        #Calculate duration (difference of end and start days PLUS scale-1)
        duration = (end - start +1) + (scale -1)

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




    ######################
    ### Find magnitude ###
    ######################

    #Calculates magnitude as monthly mean minus monthly value
    #Calculate magnitude using current and previous months according to scale parameter


    magnitude = miss_val
    count_magnitude = miss_val

    if 'magnitude' in add_metrics:

        #If found dry days
        if len(dry_days) > 0:

             #initialise
             magnitude = np.zeros(len(start))


             for k in range(len(start)):


                 #More than one consec day
                 if end[k] - start[k] > 0:

                     #Initialise temporary magnitude vector to match drought length
                     temp_mag = np.zeros(end[k]-start[k]+1)

                     #Create a vector of indices for current drought event
                     ind = np.arange(start[k], end[k]+1)


                     #Calculate magnitude for each month, taking into account scale
                     for d in range(len(ind)):
                         temp_mag[d] = sum_vec[ind[d]] - sum(mod_vec[(ind[d]-scale+1) : (ind[d]+1)]) #Need to add 1 to end day because of python indexing!!


                     #Sum monthly magnitudes to get total drought magnitude for current event
                     magnitude[k] = sum(temp_mag)



                 #One consec day only
                 else:
                     magnitude[k] = sum_vec[start[k]] - sum(mod_vec[(start[k]-scale+1) : (start[k]+1)])

        #If no dry days
        else:
            magnitude = 0.



        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'count_magnitude' in add_metrics:

            count_magnitude = np.zeros(len(count)) + miss_val

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

    #Calculates intensity as max(monthly mean minus monthly value)

    intensity = miss_val
    count_intensity = miss_val

    if 'intensity' or 'rel_intensity' in add_metrics:


       #If found dry days
        if len(dry_days) > 0:


             #initialise
             intensity = np.zeros(len(start))
             rel_intensity = np.zeros(len(start))

             for k in range(len(start)):


                 #More than one consec day
                 if end[k] - start[k] > 0:

                     #Initialise temporary intensity vector to match drought length
                     temp_int     = np.zeros(end[k]-start[k]+1)
                     temp_rel_int = np.zeros(end[k]-start[k]+1)

                     #Create a vector of indices for current drought event
                     ind = np.arange(start[k], end[k]+1)

                     #Calculate magnitude for each month, taking into account scale
                     for d in range(len(ind)):

                         #Absolute intensity
                         temp_int[d]     = sum_vec[ind[d]] - sum(mod_vec[(ind[d]-scale+1) : (ind[d]+1)]) #Need to add 1 to end day because of python indexing!!


                     #Average monthly magnitudes to get event intensity
                     intensity[k] = np.mean(temp_int)

                     #Relative intensity abs( (m - mean) / mean * 100)), where m is drought month value
                     rel_intensity[k] = abs( (sum(mod_vec[(ind[0]-scale+1) : (ind[-1]+1)]) - sum(sum_vec[ind])) /
                                        sum(sum_vec[ind] )* 100 )
                     
                     
                 #One consec day only
                 else:

                     #Absolute intensity
                     intensity[k] = sum_vec[start[k]] - sum(mod_vec[(start[k]-scale+1) : (start[k]+1)])

                     #Relative intensity abs( (m - mean) / mean * 100)), where m is drought month value
                     rel_intensity[k] = abs( (sum(mod_vec[(start[k]-scale+1) : (start[k]+1)]) - sum_vec[start[k]]) /
                                        sum_vec[start[k]] * 100 ) #Need to add 1 to end day because of python indexing!!



        #If no dry days
        else:
            intensity = 0.
            rel_intensity = 0.


        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'count_intensity' or 'count_rel_intensity' in add_metrics:

            count_intensity     = np.zeros(len(count)) + miss_val
            count_rel_intensity = np.zeros(len(count)) + miss_val

            #If found dry days
            if len(dry_days) > 0:

                for c in count:

                    #find indices for events of duration 'c'
                    ind = np.where(duration==c)[0]

                    if len(ind) > 0:
                        count_intensity[c-1]     = np.mean(intensity[ind])
                        count_rel_intensity[c-1] = np.mean(rel_intensity[ind])




    ##########################
    ### Find SPI magnitude ###
    ##########################

    #Calculates magnitude of drought-month SPI values
    #Don't need to use scale, already reflected in monthly SPI values

    #Initialise here if not wanted outputs, to avoid errors creating output dict
    spi_magnitude = miss_val
    spi_count_magnitude = miss_val

    if 'SPI_magnitude' in add_metrics:

        #If found dry days
        if len(dry_days) > 0:

             #initialise
             spi_magnitude = np.zeros(len(start)) + miss_val

             for k in range(len(start)):

                 #More than one consec day
                 if end[k] - start[k] > 0:
                     spi_magnitude[k] = abs(sum(spi_vec[start[k]:(end[k]+1)] ))  #Need to add 1 to end day because of python indexing!!
                 #One consec day only
                 else:
                     spi_magnitude[k] =  abs(spi_vec[start[k]])

        #If no dry days
        else:
            spi_magnitude = 0.



        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'SPI_count_magnitude' in add_metrics:

            spi_count_magnitude = np.zeros(len(count)) + miss_val

            #If found dry days
            if len(dry_days) > 0:

                for c in count:

                    #find indices for events of duration 'c'
                    ind = np.where(duration==c)[0]

                    if len(ind) > 0:
                        spi_count_magnitude[c-1] = np.mean(spi_magnitude[ind])




    ##########################
    ### Find SPI intensity ###
    ##########################

    #Calculates intensity of drought-month SPI values
    #Don't need to use scale, already reflected in monthly SPI values

    #Initialise here if not wanted outputs, to avoid errors creating output dict
    spi_intensity = miss_val
    spi_count_intensity = miss_val

    if 'SPI_intensity' in add_metrics:

        #If found dry days
        if len(dry_days) > 0:

             #initialise
             spi_intensity = np.zeros(len(start)) + miss_val

             for k in range(len(start)):

                 #More than one consec day
                 if end[k] - start[k] > 0:
                     spi_intensity[k] = abs(np.mean(spi_vec[start[k]:(end[k]+1)] ))  #Need to add 1 to end day because of python indexing!!
                 #One consec day only
                 else:
                     spi_intensity[k] =  abs(spi_vec[start[k]])

        #If no dry days
        else:
            spi_intensity = 0.



        #Calculate mean magnitude of drought events of certain length (as defined in "count")
        if 'SPI_count_intensity' in add_metrics:

            spi_count_intensity = np.zeros(len(count)) + miss_val

            #If found dry days
            if len(dry_days) > 0:

                for c in count:

                    #find indices for events of duration 'c'
                    ind = np.where(duration==c)[0]

                    if len(ind) > 0:
                        spi_count_intensity[c-1] = np.mean(spi_intensity[ind])

    ### List outputs ###






    outs = {'duration': duration, 'timing': timing, 'magnitude': magnitude, 'intensity': intensity,
            'rel_intensity': rel_intensity, 'count_duration': count_duration, 'count_magnitude': count_magnitude,
            'count_intensity': count_intensity, 'count_rel_intensity': count_rel_intensity,
            'SPI_magnitude': spi_magnitude, 'SPI_count_magnitude': spi_count_magnitude,
            'SPI_intensity': spi_intensity, 'SPI_count_intensity': spi_count_intensity}


    return outs;
