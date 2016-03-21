#! /usr/bin/env python

"""Module to align all data or subset of visits. 

Use:
    
    Be in the outputs directory.

    >>> python align_all.py --v (optional)

    --v : the visit number.

Authors:
    
    C.M. Gosmeyer, Feb. 2016

Example:

    >>> python check_all.py --v 44 05

Outputs:

Notes:

    Based on `clear/clear_grism.prep_clear`.
"""

import argparse
import os
import glob
from astropy.io import fits as pyfits
import numpy as np

import threedhst
import unicorn

import clear_tools
from clear_tools.set_paths import paths

import threedhst.prep_flt_astrodrizzle as init
import unicorn.interlace_test as test


def align_all(visits=[], make_asn = False):
    """
    """
    print paths
    print visits

    if make_asn:
        unicorn.candels.make_asn_files()
        # Lower case the filenames
        #asn_files = glob.glob('*asn.fits')
        #for asn in asn_files:
        #   os.rename(asn, asn.lower())

    direct_files = glob.glob('*F105W_asn.fits')
    grism_files = glob.glob('*G102_asn.fits')

    print direct_files
    
  

    for direct, grism in zip(direct_files, grism_files):
        print direct, grism

        if 'GN' in direct:
            cat_file = 'goodsn_radec.cat'
            print "Using radec cat for NORTH, {}".format(cat_file)
        elif 'GS' in direct:
            cat_file = 'goodss_3dhst.v4.1.radec.cat'
            print "Using radec cat for SOUTH, {}".format(cat_file)


        visit = int(grism[4:6])
        print "Visit: {}".format(visit)
        # Only reprocess the visits in the list, or reprocess
        # everything if the list is empty.

        if ((visits != []) and (visit in visits)) or (visits == []):

            init.prep_direct_grism_pair(direct_asn=direct, grism_asn=grism, 
                radec=paths['path_to_ref_files']+'REF/'+cat_file,
                raw_path = paths['path_to_RAW'], mask_grow=8, 
                scattered_light=False, final_scale = None,
                skip_direct=False, ACS=False, align_threshold=6)




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



if __name__=="__main__":

    args = parse_args()

    visits = args.visits

    align_all(visits=visits)