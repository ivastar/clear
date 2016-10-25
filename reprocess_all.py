#! /usr/bin/env python

"""Module to reprocess all data or subset of visits. 
Pops reads designated reads. 
Check the new FLTs against old. They should cleaner, especially if
specified a read to be removed.

*Step 3 of Prep.*

Use:

    Be in RAW directory.
    
    >>> python reproces_all.py --v (optional)

    --v : the visit number.

Authors:
    
    C.M. Gosmeyer, Feb. 2016

Example:

    >>> python reprocess_all.py --v 44 05

Outputs:


Notes:

    Based on `unicorn.prepare_fixed_flt.py`.
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
    """ Assumes in RAW directory.
    """
    if os.path.isfile('files.info'):

        files, filters = np.loadtxt('files.info', unpack=True, \
                                              usecols=(0,5), dtype='S24') 

    else: 
        print "'files.info' does NOT exist."

    return files, filters


#-------------------------------------------------------------------------------  

def reprocess_all(visits=[]):
    """
    Big list of all of the FLTs reprocessed from the IMA files

    Satellites:
    icxt04e4q: read 6
    icxt05h7q: read 1
    icxt31r3q: read 2
    icxt63lzq: read 10

    Earth shine:
    44,45
    (for now don't pop reads)

    """
    print visits
    # First read in the filelist. 
    filenames, filters = read_filesinfo()

    # Only retain those wtih grisms.
    grism_indices = np.where(filters == 'G102')
    filenames = filenames[grism_indices]

    for filename in filenames:
        print "Filename: {}".format(filename)
        pop_reads=[]
        visit = int(filename[4:6])
        print "Visit: {}".format(visit)
        # Only reprocess the visits in the list, or reprocess
        # everything if the list is empty.
        if ((visits != []) and (visit in visits)) or (visits == []):
            rawname = filename.split('_flt.fits')[0] + '_raw.fits'
            print "Reprocessing raw file {}".format(rawname)
            if 'icxt04e4q' in rawname:
                # Satelite 
                pop_reads=[6]
            elif 'icxt05h7q' in rawname:
                # Satelite
                pop_reads=[1]
            elif 'icxt31r3q' in rawname:
                # Satelite
                pop_reads=[2]
            elif 'icxt63lzq' in rawname:
                # Satelite
                pop_reads=[10]
            elif 'icxt34ekq' in rawname:
                # Satelite?
                pop_reads=[1]
            elif 'icxt16jzq' in rawname:
                # CR explosion?
                pop_reads=[6]
            elif 'icxt01ciq' in rawname:
                # Satelite
                pop_reads=[9]
            elif 'icxt01clq' in rawname:
                # streak-thing? 
                pop_reads=[2]
            elif 'icat07bxq' in rawname:
                # 3DHST, satelite.
                pop_reads=[3]
            elif 'icat21cfq' in rawname:
                # 3DHST, CR barf
                pop_reads=[6]
            elif 'icat11c2q' in rawname:
                # 3DHST, CR thing
                pop_reads=[6]
            elif 'icxt27hhq' in rawname:
                # Dispersed satellite trail. Possibly persisting into next read.
                pop_reads=[3]
            elif 'icxt64qvq' in rawname:
                # Dispersed satellite trail.
                pop_reads=[7]
            elif 'icxt50i1q' in rawname:
                # Satellite trail.
                pop_reads=[12]
            elif 'icxt50hwq' in rawname:
                # Satellite trail.
                pop_reads=[3]
            elif 'icxt61gzq' in rawname:
                # Satellite trail?
                pop_reads=[4]
            elif 'icxt61gwq' in rawname:
                pop_reads=[5,6,7,8,9,10,11]

            reprocess_wfc3.make_IMA_FLT(raw=rawname, pop_reads=pop_reads)


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
#-------------------------------------------------------------------------------  

if __name__=="__main__":

    args = parse_args()

    visits = args.visits

    reprocess_all(visits=visits)