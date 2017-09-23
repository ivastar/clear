
"""
Utilities for CLEAR grizli pipeline.

Author:

    C.M. Gosmeyer, Sept 2017

"""

import glob
import numpy as np
import os
import sys
import time


def make_timestamp_dir(dest):
    """Creates time-stamped directory. YYYY.MM.DD.hh.mm.ss<_num>

    Parameters
    ----------
    dest : string
        Path to where the time-stamp directory should be created.

    Returns
    -------
    path_to_time_dir : string
        Path to and including the time-stamped directory.

    Outputs
    -------
    Directory at 'dest' with a time-stamped name.

    """
    time_tuple = time.localtime()

    time_array = np.array(time_tuple, dtype='str')
    time_array = np.array(['0'+t if len(t)==1 else t for t in time_array])     
        
    time_dir = '{}.{}.{}.{}.{}.{}'.format(time_array[0],time_array[1],
        time_array[2],time_array[3],time_array[4],time_array[5])
    path_to_time_dir = os.path.join(dest, time_dir)
    
    # If one does not exist for today, create the time-stamp dir.
    if not os.path.isdir(path_to_time_dir):
        os.mkdir(path_to_time_dir)
        return path_to_time_dir
    
    # If one already exists, create it with underscore index 1-100.
    else:
        for num in range(1,100):
            path_to_time_dir_num = path_to_time_dir + '_' + str(num)
            if not os.path.isdir(path_to_time_dir_num):
                os.mkdir(path_to_time_dir_num)
                return path_to_time_dir_num    


#-------------------------------------------------------------------------------

def store_outputs(path_outputs, store_type):
    """
    This does not deal with sorting extractions. Only the output from 
    prep and interlacing steps. 
    It should be easy to go into this directory and recreate the needed
    extractions, functioning just as outputs always does, only this is
    a frozen version of outputs that cannot be overwritten.


    Parameters
    ----------
    store_type : string
        Either 'prep' or 'inter', to be appended to directory name.

    Returns
    -------
    path_outputs_timestamp : string
        Path to the outputs

    """
    path_outputs = make_timestamp_dir(path_outputs)
    path_outputs_split = path_outputs.split('/')
    #path_outputs_split[-1] = path_outputs_split[-1] + '_' + store_type
    path_outputs_timestamp = '/'.join(path_outputs_split)

    return path_outputs_timestamp


#-------------------------------------------------------------------------------

def retrieve_latest_outputs(path_outputs, store_type):
    """
    Returns path+name to the most recent output sub-directory for 'store_type'.

    Parameters
    ----------
    store_type : string
        Either 'prep' or 'inter', appended to directory name.

    Returns
    -------

    """
    # glob over the directories, return the one with greatest (most recent) 
    # creation time.
    directories = glob.glob(os.path.join(path_outputs, '*{}'.format(store_type)))
    max_birthday = 0
    for directory in directories:
        birthday = os.path.getctime(directory)
        if birthday > max_birthday:
            max_birthday = birthday

    return max_birthday





