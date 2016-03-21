#! /usr/bin/env python

""" Determines from timestamps whether a pointing started with grism or 
a direct image.

Authors:
    
    C.M. Gosmeyer, 29 Feb. 2016


Outputs:

"""

import numpy as np
import os
from astropy.time import Time
from astropy.io import ascii
from clear_tools import set_paths

def find_pointing_start(asn_table_name):
    """

    Parameters:
        asn_table_name : string
            For example, 'gs2-01-189-g102'. Targname-visit-PA-filter

    Returns:
        0 : if visit starts with direct image
        1 : if visit starts with grism

    Outputs:
    
    """
    # Parse out the asn tables's target and visit

    path_to_outputs = set_paths.paths['path_to_outputs']

    # Read in files.info.
    data = ascii.read(os.path.join(path_to_outputs, 'files.info'))

    # Makes lists out of FILE, TARGNAME, TIME-OBS
    filenames = np.array(data['FILE'])
    targnames = np.array(data['TARGNAME'])
    dateobses = np.array(data['DATE-OBS'])
    timeobses = np.array(data['TIME-OBS'])
    filters = np.array(data['FILTER'])

    # Sort out only the rows with same target and visit as asn table.
    targname_sort = np.where(targnames == asn_table_name.split('-')[0].upper())

    visit_sort = []
    for filename, i in zip(filenames[targname_sort], np.arange(len(filenames[targname_sort]))):
    	if filename[4:6] == asn_table_name.split('-')[1]:
    		visit_sort.append(i)

    visit_sort = np.array(visit_sort)
    dateobses_sorted = (dateobses[targname_sort])[visit_sort]
    timeobses_sorted = (timeobses[targname_sort])[visit_sort]
    filters_sorted = (filters[targname_sort])[visit_sort]

    # Now loop through time-obs list and find which entry came very first. 
    # Set the first_image initially to be first entry. (and in most cases 
    # this should remain true after loop)
    first_time = Time(dateobses_sorted[0] + ' ' + timeobses_sorted[0], format='iso', out_subfmt='date_hms')
    first_image = filters_sorted[0]
    print "Initial start image: {} and start time: {}".format(first_image, first_time)
    for dateobs, timeobs, filt in zip(dateobses_sorted[1:], timeobses_sorted[1:], filters_sorted[1:]):
    	# Turn into a date-time object.
        next_time = Time(dateobs + ' ' + timeobs, format='iso', out_subfmt='date_hms')
        min_time = min(first_time, next_time)
        if min_time == next_time:
        	first_time = next_time
        	first_image = filt
        	print "New start image: {} and start time: {}".format(first_image, first_time)

    # If a Direct image (filter=f105w) taken first, return 0.
    if first_image.lower() == 'f105w':
        return 0
    # If a Grism (filter=g102) taken first, return 1.
    elif first_image.lower() == 'g102':
	    return 1
