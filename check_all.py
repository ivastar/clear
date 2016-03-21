#! /usr/bin/env python

"""Module to check all data or subset of . 

Use:
    
    >>> python check_all.py --v (optional)

    --v : 

Authors:
    
    C.M. Gosmeyer, Feb. 2016

Example:

    >>> python check_all.py --v 44 05

Outputs:


"""

from astropy.io import ascii
import numpy as np
import argparse
import glob
import os
import shutil
import mywfc3
from mywfc3 import reprocess_wfc3
#import unicorn


#-------------------------------------------------------------------------------  

def read_filesinfo():
    """ Assumes in RAW directory.
    """
    if os.path.isfile('files.info'):
        
        files, filters = np.loadtxt('files.info', unpack=True, \
                                              usecols=(0,5), dtype='S24') 

    else: 
        print "'files.info' does NOT exist."

    return files, filters

#-------------------------------------------------------------------------------  

def check_all(visits=[]):
    """
    """
    print visits
    # First read in the filelist. 
    filenames, filters = read_filesinfo()

    for filename in filenames:
    	imaname = (filename.split('_flt'))[0] + '_ima.fits'
        print "IMA name: {}".format(imaname)
        visit = int(imaname[4:6])
        print "Visit: {}".format(visit)

        if ((visits != []) and (visit in visits)) or (visits == []):
            reprocess_wfc3.show_MultiAccum_reads(raw = imaname)
             

    pngs = glob.glob("*png")
    for png in pngs:
        print "Moving {}".format(png)
        shutil.move(png, '../checks/')

#-------------------------------------------------------------------------------  

def parse_args():
    """Parses command line arguments.
    
    Parameters:
        nothing
        
    Returns:
        args : object
            Containing the image and destination arguments.
            
    Outputs:
        nothing.
    """

    visits_help = "List of visits to loop over. Default is to loop over all."
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--v', dest = 'visits',
                        action = 'store', type = int, required = False,
                        help = visits_help, nargs='+', default=[])       
    args = parser.parse_args()
     
    return args

#-------------------------------------------------------------------------------  

if __name__=="__main__":

    args = parse_args()

    visits = args.visits

    check_all(visits=visits)