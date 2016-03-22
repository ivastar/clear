#! /usr/bin/env python

"""Module to align all data or subset of visits. 
It runs TWEAKREG and aligns the direct images to a catalog, flags cosmic
rays, and subtracts background.
For the grisms, it does a more complicated background subtraction (see 
Gabe's ISR). If persistence files available, it will mask the affected 
pixels.
When it is complete, should look at all the FLTs. They all should be 
flatfielded, have background substracted (so average is about 0). No 
cartwheel should be visible (but Deathstar still there). No streaks. 
No gradients.


*Step 4 of Prep.*

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
import shutil
import numpy as np

import threedhst
import unicorn

from clear.set_paths import paths

import threedhst.prep_flt_astrodrizzle as init
import unicorn.interlace_test as test


def cleanup_outputs():
    """ Removes all extraneous files from last run of this script.
    """

    path_to_outputs = paths['path_to_outputs']

    for f in ['sex_stderr', 'threedhst_auto.param', 'threedhst_auto.sex',
              'threedhst_auto-2.param', 'threedhst_auto-2.sex', 
              'tweakreg.log', 'astrodrizzle.log', 'default.nnw', 
              'gauss_4.0_7x7.conv', 'grism.conv']
        f = path_to_outputs + f
        if os.path.isfile(f):
            os.remove(f)


def align_all(visits=[], make_asn = False):
    """
    """
    print paths
    print visits

    shutil.copy('../RAW/files.info', 'files.info')
    cleanup_outputs()

    if make_asn:
        unicorn.candels.make_asn_files()
        # Lower case the filenames
        #asn_files = glob.glob('*asn.fits')
        #for asn in asn_files:
        #   os.rename(asn, asn.lower())

    direct_files = glob.glob('*F105W_asn.fits')

    print direct_files

    for direct in direct_files:
        grism = direct.replace('F105W', 'G102')
        print " "
        print direct, grism

        visit = int(grism[4:6])
        print "Visit: {}".format(visit)
        # Only reprocess the visits in the list, or reprocess
        # everything if the list is empty.

        if ((visits != []) and (visit in visits)) or (visits == []):

            if 'GN' in direct:
                cat_file = 'goodsn_radec.cat'
                print "Using radec cat for NORTH, {}".format(cat_file)
            elif 'GS' in direct:
                cat_file = 'goodss_3dhst.v4.1.radec.cat'
                print "Using radec cat for SOUTH, {}".format(cat_file)

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

    align_all(visits=visits, make_asn = True)