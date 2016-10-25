#! /usr/bin/env python

"""Module to retrieve new FITS and persistence products from /wfc3 
central store and copy them into /astro directories.
Gzips all the original FLTs and moves them to 'flt_zipped' directory.

*Step 1 of Prep.*

Authors:
    
    C.M. Gosmeyer, Mar. 2016

Use:
    
    Be in the RAW directory.

    >>> python retrieve_all.py --v (optional)

    --v : the visit number.

Example:

    >>> python retrieve_all.py --v 44 05

"""

import argparse
import glob
import os
import shutil
import subprocess

from set_paths import paths


#-------------------------------------------------------------------------------  

def make_filesinfo():
    """ Creates files.info file in RAW.
    """
    # Retain the old files.info. User can delete if desired.
    path_to_RAW = paths['path_to_RAW']
    path_to_software = paths['path_to_software']
    if os.path.isfile(path_to_RAW + 'files.info'):
        os.rename(path_to_RAW + 'files.info', path_to_RAW + 'files.info.retain')
    cwd = os.getcwd()

    os.chdir('{}threedhst/threedhst/bin/'.format(path_to_software))
    subprocess.call(['./flt_info.sh'])
    os.chdir(cwd)
    shutil.copy(path_to_software+'threedhst/threedhst/bin/files.info', path_to_RAW + 'files.info')


#-------------------------------------------------------------------------------  

def retrieve_all(visits=[], overwrite=True):
    """

    Parameters:
        overwrite : {True, False}
            Set to True if want to overwrite existing files.
            False by default.
    """

    wfc3_centralstore = '/grp/hst/wfc3a/GO_Links/14227/'
    path_to_RAW = paths['path_to_RAW']
    path_to_PERSIST = paths['path_to_PERSIST']

    # Glob over all the Visit directories in /wfc3 central store.
    visit_dirs = glob.glob(wfc3_centralstore+'Visit*')

    for visit_dir in visit_dirs:
       
        visit = int(((visit_dir.split('/'))[len(visit_dir.split('/'))-1])[5:7])
        print "Visit: {}".format(visit)

        # Only retrieve the visits in the list, or reprocess
        # everything if the list is empty.
        if ((visits != []) and (visit in visits)) or (visits == []):
            fits_list = glob.glob(visit_dir+'/*fits')
            for fits_path in fits_list:
                fits = fits_path.split('/')[-1]
                #print "Copy from: {}".format(fits_path)
                #print "Copy to: {}".format(path_to_RAW + fits)
  
                # Check whether file already exists. 
                if os.path.isfile(path_to_RAW+fits) and overwrite:
                    # Copy into RAW
                    print "Overwriting {}.".format(fits)
                    shutil.copy(fits_path, path_to_RAW + fits)
                    if 'flt.fits' in fits:
                        os.system('gzip {}'.format(path_to_RAW + fits))
                        shutil.copy(fits_path, path_to_RAW + fits)
                elif os.path.isfile(path_to_RAW+fits) and not overwrite:
                    print "Not copying {} because will overwrite.".format(fits)
                else:
                    print "Copying {}.".format(fits)
                    shutil.copy(fits_path, path_to_RAW + fits)
                    # gzip all the FLTs
                    if 'flt.fits' in fits:
                        os.system('gzip {}'.format(path_to_RAW + fits))
                        shutil.copy(fits_path, path_to_RAW + fits)

            # Copy persistence files in PERSIST
            persist_list = glob.glob(visit_dir+'/Persist/*persist.fits')
            print persist_list
            for persist_path in persist_list:
                persist = persist_path.split('/')[-1]
                shutil.copy(persist_path, path_to_PERSIST + persist)
                print "Copying {}.".format(persist)
        
    # Move all *flt.fits.gz to flt_zipped dir.
    #(turning off for now because align_all.py wants them in main dir now...)
    #gzip_list = glob.glob(path_to_RAW + '*flt.fits.gz')
    #for gzip in gzip_list:
    #    shutil.move(gzip, path_to_RAW + 'flt_zipped/')

    make_filesinfo()


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

    visits_help = "List of visits to retrieve from wfc3 central store. Default is all."
        
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

    retrieve_all(visits=visits)
    make_filesinfo()
