#! /usr/bin/env python

"""

Author:

    C.M. Gosmeyer, Aug. 2017



"""
from __future__ import print_function

import argparse
import grizli
import drizzlepac
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import time

from astropy.io import fits
from set_paths import paths
#from logging import logging

overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11']}


#------------------------------------------------------------------------------- 

def find_files():
    """
    """
    files = glob.glob(os.path.join(paths['path_to_RAW'], '*flt.fits')) \
        + glob.glob(os.path.join(paths['path_to_RAW_3DHST'], '*flt.fits'))
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)

    for i in range(len(visits)):
        print(dict(visits[i]))

    return visits, filters

#------------------------------------------------------------------------------- 

def clear_grizli_pipeline():
    """ Main wrapper on pre-processing, interlacing and extracting steps.
    """

    visits, filters = find_files()


#-------------------------------------------------------------------------------  

def parse_args():
    """Parses command line arguments.
        
    Returns
    -------
    args : object
        Containing the image and destination arguments.
            
    """

    fields_help = "List the fields over which to run pipeline. Default is all. "
    ref_filter_help = "The reference image filter. Choose either F105W or F125W. Default is F105W. "
    mag_lim_help = "The magnitude limit for extraction. Default is 25."
    do_steps_help = "List the processing steps to run. Default is all five. Steps are NOT independent. " 
    do_steps_help += "If choose to run a step alone, be sure products exist from the previous step. "
    do_steps_help += "1 - Interlace visits. "
    do_steps_help += "2 - Create contamination model. "
    do_steps_help += "3 - Extract traces. "
    do_steps_help += "4 - Stack traces. "
    do_steps_help += "5 - Fit redshifts and emission lines of traces. "
    cats_help = "List of catalogs over which to run pipeline. Use in combination with mag_lim. "
    cats_help += "Default is 'full', which is generally used when extracting by mag_lim. "
    cats_help += "Catalog options are 'full', and its subsets, 'emitters', 'quiescent', 'sijie', and 'zn'. "
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--fields', dest = 'fields',
                        action = 'store', type = str, required = False,
                        help = fields_help, nargs='+', default=['GS1', 'GS2', 'GS3', 'GS4', 'GS5', 'ERSPRIME', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])       
    parser.add_argument('--steps', dest = 'do_steps',
                        action = 'store', type = int, required = False,
                        help = do_steps_help,  nargs='+', default=[1,2,3,4,5])    
    parser.add_argument('--mlim', dest = 'mag_lim',
                        action = 'store', required = False,
                        help = mag_lim_help, default=25)
    parser.add_argument('--ref', dest = 'ref_filter',
                        action = 'store', type = str, required = False,
                        help = ref_filter_help,  default='F105W')
    parser.add_argument('--cats', dest = 'cats',
                        action = 'store', type = str, required = False,
                        help = cats_help, nargs = '+', default='full')

    args = parser.parse_args()
     
    return args

#------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------- 

if __name__=='__main__':
    clear_grizli_pipeline()


