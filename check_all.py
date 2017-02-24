#! /usr/bin/env python

"""Module to create diagnostic PNGs for checking ramps in a visit or 
subset of visits. 

*Step 2 of Prep.*

Checks:

    Flip through the PNGs. If bad reads found, you will need to mark them 
    to be popped in :mod:`reprocess_all.py`.
    For earth-flat persistence, or any anamoly that affects all the reads,
    will need to create masks.


HOW TO CREATE MASKS
++++++++++++++++++++
- Be in RAW.
- In ipython,
  > import threedhst.dq
  > dq = threedhst.dq.checkDQ(asn_direct_file='icxt07020_asn.fits', asn_grism_file='icxt07020_asn.fits', path_to_flt='./')
- Once DS9 opens, select the FLT then,
  "flag grism"
  draw polygons
  "flag grism" again
  A *reg file for that FLT will be saved.

HOW TO APPLY MASKS
+++++++++++++++++++
- Masks are applied in the prep step
   threedhst.prep_flt_astrodrizzle.prep_direct_grism_pair
   which I wrapped in :mod:`align_all.py`.

Use:
    
    Be in RAW directory.

    >>> python check_all.py --v (optional)

    --v : the visit number.

Authors:
    
    C.M. Gosmeyer, Feb. 2016

Example:

    >>> python check_all.py --v 44 05 A4

"""

import argparse
import glob
import numpy as np
import os
import shutil

from astropy.io import ascii

import mywfc3
from mywfc3 import reprocess_wfc3


#-------------------------------------------------------------------------------  

def read_filesinfo():
    """ Reads in the "files.info" file. Assumes you are in the RAW directory.

    Returns
    -------
    files : numpy array
        Names of the FLT files.
    filters : numpy array
        Names the files' filters. 

    """
    if os.path.isfile('files.info'):
        
        files, filters = np.loadtxt('files.info', unpack=True, \
                                              usecols=(0,5), dtype='S24') 

    else: 
        print "'files.info' does NOT exist."

    return files, filters


#-------------------------------------------------------------------------------  

def check_all(visits=[]):
    """ Creates diagonstic PNGs of the ramps of each file in the given visits.
    To be used for searching for satellite trails, earth-limb, etc., so can 
    either pop the read (if affects single read) or create a mask (if affects
    all the reads). 

    Parameters
    ----------
    visits : list of ints/strings
        The two-digit numbers of the visits whose observations' ramps will 
        be printed to PNGs. e.g., [01, 11, 'a4']    
    """
    print visits
    # First read in the filelist. 
    filenames, filters = read_filesinfo()

    for filename in filenames:
    	imaname = (filename.split('_flt'))[0] + '_ima.fits'
        print "IMA name: {}".format(imaname)
        visit = imaname[4:6]
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
        
    Returns
    -------
    args : object
        Containing the image and destination arguments.
    """

    visits_help = "List of visits to loop over. Default is to loop over all."
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--v', dest = 'visits',
                        action = 'store', type = str, required = False,
                        help = visits_help, nargs='+', default=[])       
    args = parser.parse_args()
     
    return args


#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------  

if __name__=="__main__":

    args = parse_args()

    visits = args.visits

    check_all(visits=visits)