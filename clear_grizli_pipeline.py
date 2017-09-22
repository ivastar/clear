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
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

from grizli.prep import process_direct_grism_visit
from grizli import config
#from grizli.config import path_raw, path_outputs, path_persistence
from analysis_tools.tables import bypass_table 
from astropy.io import fits
from astropy.table import Table
from set_paths import paths
from utils import store_outputs, retrieve_latest_outputs
from record_imports.record_imports import meta_record_imports
from record_imports.log import setup_logging, log_info, log_metadata

path_raw = paths['path_to_RAW']
path_outputs = paths['path_to_outputs']
path_ref = os.path.join(paths['path_to_ref_files'], 'REF')


# Is grizli smart enough to know that these program 13420 fields overlap with 
# the given CLEAR fields?
overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11']}


#-------------------------------------------------------------------------------


# put this all into a huge hideous class, clear_grizli_pipeline? 

def find_files(fields=['GN2']):
    """

    Parameters
    ----------
    fields : list of strings
        The CLEAR fields; retain ability to specify the individual pointings
        so that can easily re-run single ones if find an issue.

    Returns
    -------
    visits : OrderedDict
        Keys of 'files' and 'products'; values of list of FLT files and product name.
    filters : OrderedDict
        Keys of filters; values of OrderedDicts with keys of orient and values of
        lists of FLT files.

    """
    files = glob.glob(os.path.join(paths['path_to_RAW'], '*flt.fits')) \
        + glob.glob(os.path.join(paths['path_to_RAW_3DHST'], '*flt.fits'))
    info = grizli.utils.get_flt_info(files)
    # 'info' is an astropy table.

    # Creating a new table and inserting only the rows I want is quite annoying
    # So to conserve our sanity, just convert 'info' into an ordered dictionary
    info_dict = bypass_table.decompose_table(info, return_type=dict, include_meta=True)

    new_info_list = []

    # Convert table to ordered dictionary, put it in a list, convert to numpy array
    # to convert it back into a table. awesome
    for row in range(len(info_dict['TARGNAME'])):
        if info_dict['TARGNAME'][row] in fields:
            new_info_list.append([info_dict[key][row] for key in info_dict.keys() if key != 'meta'])
        # Now, if one of the fields is GOODS-N, need be sure all the Barro programs
        # are also retained. 
        for field in [field for field in fields if 'N' in field]:
            if info_dict['TARGNAME'][row] in overlapping_fields[field]:
                new_info_list.append([info_dict[key][row] for key in info_dict.keys() if key != 'meta'])
                # Break so the Barro programs are added only once.
                break

    # Convert 'info' back into a table
    # I couldn't simply do dtype=info.dtype, so instead of hitting my head on the
    # keyboard for hours, just hard code it
    new_info_tab = Table(np.array(new_info_list), names=info.colnames, meta=info.meta, 
        dtype=['S18', 'S5', 'S8', 'S10', 'S8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'])

    # huzzah it works
    visits, filters = grizli.utils.parse_flt_files(info=new_info_tab, uniquename=True)

    return visits, filters

#-------------------------------------------------------------------------------

@log_metadata
def prep(visits, prime_filt='F105W', prime_grism='G102'):
    """
    This is akin to `align_all.py`. It copies the FLTs from 'RAW/' and 
    performs background subtraction and flat-fielding, extracts visit-level
    catalogs and seg maps from direct images, and produces aligned, background-
    subtracted FLTs and drizzled mosaics of both direct and grism images.

    Parameters
    ----------

    Outputs 
    -------
    The individual visit outputs are

    * icat21dlq_crclean.fits
    [are the FLTs modified?]

    In the directory specified by path_to_outputs the combined for direct 
    images are

    * gn2-cxt-51-345.0-f105w_asn.fits
    * gn2-cxt-51-345.0-f105w_bkg.fits
    * gn2-cxt-51-345.0-f105w_drz_sci.fits - drizzled direct mosaic
    * gn2-cxt-51-345.0-f105w_drz_wht.fits 
    * gn2-cxt-51-345.0-f105w_seg.fits
    * gn2-cxt-51-345.0-f105w_wcs.fits 
    * gn2-cxt-55-022.0-f105w.cat   
    * gn2-cxt-55-022.0-f105w.cat.radec
    * gn2-cxt-55-022.0-f105w.cat.reg 
    * gn2-cxt-53-309.0-f105w_shifts.log 

    For grism images

    * gn2-cxt-51-345.0-g102_asn.fits 
    * gn2-cxt-51-345.0-g102_drz_sci.fits
    * gn2-cxt-51-345.0-g102_drz_wht.fits 
    * gn2-cxt-53-309.0-g102_1_sky_background.info 

    Astrometry:
    * _wcs.png - diagnostic of [...]
    * _wcs.log - 
    * _wcs.fits - 

    Grism sky background subtraction:
    * _column.png - diagnostic of [...]
    * _column.dat - produces PNG plot




    Aligned, background-subtracted FLTs 

    Drizzled mosaics


    [which are analogs of the old pipeline's files?]

    """
    

    # Match the direct and the grism visits.
    # Going by order as in the example won't work. 
    # Might be able to match 'product'.split('.')[0] values from 'visit' dictionary
    print(visits)

    # this ain't right
    for visit1 in visits:
        product1 = visit1['product']
        filt1 = product1.split('.')[1]
        basename1 = product1.split('.')[0]
        if prime_filt.lower() in filt1.lower():
            # If the filter indicates stare images, search for the grisms
            for visit2 in visits:
                # this really ain't right 
                # should be able to pull this out of the dictionary as keys
                product2 = visit2['product']
                filt2 = product2.split('.')[1]
                basename2 = product2.split('.')[0]
                field = product2.split('-')[0]
                if basename1 == basename2 and prime_grism.lower() in filt2.lower():
                    # Point to correct, field-dependent radec catalog
                    print("Matched direct to grism products: {} {}".format(product1, product2))
                    if 'n' in field.lower():
                        radec_catalog = 'goodsn_radec.cat'
                        print("Using radec cat for NORTH: {}".format(radec_catalog))
                    elif 's' in field.lower():
                        radec_catalog = 'goodss_3dhst.v4.1.radec.cat'
                        print("Using radec cat for SOUTH: {}".format(radec_catalog))

                    # Do the prep steps.
                    status = process_direct_grism_visit(
                        direct=visit1,
                        grism=visit2,
                        radec=os.path.join(path_ref, radec_catalog),
                        path_raw=path_raw,
                        align_mag_limits=[14,23])


#-------------------------------------------------------------------------------

@log_info
@log_metadata
def clear_grizli_pipeline(path_outputs_timestamp):
    """ Main wrapper on pre-processing, interlacing and extracting steps.
    """
    # cd into outputs directory, as running grizli requires
    os.chdir(path_outputs_timestamp)
    logging.info("cd into {}".format(path_outputs_timestamp))

    # Find the files in RAW
    visits, filters = find_files()

    # In outputs run the prepsteps
    prep(visits=visits, prime_filt='F105W', prime_grism='G102')

    # Do the interlacing; need have option which outputs subdir to use? (nominally, all will be same)

    # Do the extractions; need have option which outputs subdir to use? 

    # Do the fitting; need have option which outputs subdir to use? 


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
                        help = fields_help, nargs='+', 
                        default=['GS1', 'GS2', 'GS3', 'GS4', 'GS5', 'ERSPRIME', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])       
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
    # Create the output directory
    path_outputs_timestamp = store_outputs(path_outputs=path_outputs, 
        store_type='')
    # Setup logging to save in the output directory
    setup_logging(__file__, path_logs=path_outputs_timestamp)
    #meta_record_imports(__file__, print_or_log='print') #but how direct this to a log file?
    # Call the main pipeline function.
    clear_grizli_pipeline(path_outputs_timestamp=path_outputs_timestamp)


