#! /usr/bin/env python

""" Grizli pipeline for the CLEAR collaboration.  

Author:

    C.M. Gosmeyer, Aug. 2017

Notes:

    * This loops over CLEAR 'fields' (GN1, GS1, ERSPRIME, etc).
      It's not actually necssary, since grizli is smart enough to
      match overlapping exposures without grouping them yourself.
      But I decided to allow the user to select the specific fields
      to make easier re-running troublesome fields and easier for 
      me to process fields on my laptop (doing the full chunks of 
      Goods-S or Goods-N would set it on fire). The user is still
      able to do all fields by default, however. See the arg parse
      options.

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
imoprt pickle
import shutil
import time

from astropy.io import fits
from astropy.table import Table
from set_paths import paths
from utils import store_outputs, retrieve_latest_outputs, tobool
from clear_inspection_tools import flt_residuals

# My grizli fork and other personal packages
from analysis_tools.tables import bypass_table 
from grizli.prep import process_direct_grism_visit
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
from record_imports.record_imports import meta_record_imports
from record_imports.log import setup_logging, log_info, log_metadata

# Set global paths.
PATH_RAW = paths['path_to_RAW']
PATH_OUTPUTS = paths['path_to_outputs']
PATH_REF = os.path.join(paths['path_to_ref_files'], 'REF')

global PATH_OUTPUTS_TIMESTAMP

# should overlapping_fields and pointing be put in their own module?
overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11']}

class Pointing():
    """ Generalization of GN1, GS1, ERSPRIME, etc
    """
    def __init__(self, field, ref_filter):
        if 'N' in field.upper():
            self.pad = 500 # really only necessary for GDN
            self.radec_catalog = 'goodsn_radec.cat'
            self.seg_map = 'Goods_N_plus_seg.fits'
            if '125' in ref_filter:
                self.catalog = 'GoodsN_plus_merged.cat'
                self.ref_image = 'goodsn_3dhst.v4.0.F125W_orig_sci.fits'
            elif '105' in ref_filter:
                self.catalog = 'goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat'
                self.ref_image = 'goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'
                
        elif 'S' in field.upper():
            self.pad = 200 # grizli default
            self.radec_catalog = 'goodss_3dhst.v4.1.radec.cat'
            self.seg_map = 'Goods_S_plus_seg.fits'
            if '125' in ref_filter:
                self.catalog = 'GoodsS_plus_merged.cat'
                self.ref_image = 'goodss_3dhst.v4.0.F125W_orig_sci.fits'  
            elif '105' in ref_filter:
                self.catalog = 'goodss-F105W-astrodrizzle-v4.3_drz_sub_plus.cat'
                self.ref_image = 'goodss-F105W-astrodrizzle-v4.3_drz_sci.fits'
            

#-------------------------------------------------------------------------------


def find_files(fields):
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
    # Just convert 'info' into an ordered dictionary
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
    # I couldn't simply do dtype=info.dtype, so just hard code it
    new_info_tab = Table(np.array(new_info_list), names=info.colnames, meta=info.meta, 
        dtype=['S18', 'S5', 'S8', 'S10', 'S8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'])

    visits, filters = grizli.utils.parse_flt_files(info=new_info_tab, uniquename=True)

    return visits, filters


#-------------------------------------------------------------------------------

@log_metadata
def prep(visits, ref_filter='F105W', ref_grism='G102'):
    """
    This is akin to `align_all.py`. It copies the FLTs from 'RAW/' and 
    performs background subtraction and flat-fielding, extracts visit-level
    catalogs and seg maps from direct images, and produces aligned, background-
    subtracted FLTs and drizzled mosaics of both direct and grism images.

    Parameters
    ----------
    visits :

    ref_filter : string
        The reference image's filter.
    ref_grism : string
        The grism.

    Outputs 
    -------
    The individual visit outputs are

    * icat21dlq_crclean.fits
    [are the FLTs modified?]

    In the directory specified by path_to_outputs the combined for direct 
    images are

    * flt.fits - modified, only outputs needed for model step

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

    for visit1 in visits:
        product1 = visit1['product']
        filt1 = product1.split('.')[1]
        basename1 = product1.split('.')[0]
        if ref_filter.lower() in filt1.lower():
            # If the filter indicates stare images, search for the grisms
            for visit2 in visits:
                # would be better to pull this out of the dictionary as keys
                product2 = visit2['product']
                filt2 = product2.split('.')[1]
                basename2 = product2.split('.')[0]
                field = product2.split('-')[0]
                if basename1 == basename2 and ref_grism.lower() in filt2.lower():
                    # Point to correct, field-dependent radec catalog
                    print("Matched direct to grism products: {} {}".format(product1, product2))
                    p = Pointing(field=field, ref_filter=ref_filter)
                    radec_catalog = p.radec_catalog
                    print("Using radec catalog: {}".format(radec_catalog))

                    # Do the prep steps.
                    logging.info(" ")
                    logging.info("process_direct_grism_visit on direct {} and grism {}"\
                        .format(visit1, visit2))
                    status = process_direct_grism_visit(
                        direct=visit1,
                        grism=visit2,
                        radec=os.path.join(PATH_REF, radec_catalog),
                        path_raw=PATH_RAW,
                        align_mag_limits=[14,23])

    # need make some inspection tools


#-------------------------------------------------------------------------------

@log_metadata
def model(visits, field='', ref_filter='', use_prep_path='.', load_only=False):
    """ Model the contamination.

    Parameters
    ----------
    visits :

    fields : list of strings
        Generally will be the same as in the prep step, unless for some reason
        you wish to run over a subset.
    ref_filter : string
        The reference image's filter.
    use_prep_path : string
        Path to the pre-processed files. By default the working directory.

    Outputs
    -------
    * icat08hiq_flt.01.wcs.fits
    * icat08hiq.01.GrismFLT.fits
    * icat08hiq.01.GrismFLT.pkl - contains what? 


    """

    grism_files = []
    all_flt_files = []
    for i in range(len(visits)):
        # e.g., visits[i]['product'] = 'gn2-cxt-51-345.0-g102'  
        print(visits[i]['product'])
        field_in_contest = visits[i]['product'].split('-')[0].upper()
        print(field_in_contest)   
        if '-g1' in visits[i]['product']:
            # Only add to list the grism files IF they are among the specified
            # fields
            if field_in_contest in overlapping_fields[field] or \
                field_in_contest == field:
                grism_files.extend(visits[i]['files'])
                all_flt_files.extend(visits[i]['files'])
        elif '-f1' in visits[i]['product']:
            if field_in_contest in overlapping_fields[field] or \
                field_in_contest == field:
                all_flt_files.extend(visits[i]['files'])

    p = Pointing(field=field, ref_filter=ref_filter)
    
    # If modeling on a later run than when generating prep step files,
    # copy the pre-processed FLTs to the current time-stamp directory.
    # Brammer says only *FLTs produced by prep are needed, since for CLEAR
    # we are providing the ref_file, seg_file, and catalog.

    if use_prep_path != '.':
        all_flt_files = [os.path.join(use_prep_path, flt) for flt in all_flt_files]
        # Copy only the FLT files for the given field.
        for flt in all_flt_files:
            # assuming have cd'd into the outputs location
            logging.info("Copying {} from {} to {}"\
                .format(flt, use_prep_path, PATH_OUTPUTS_TIMESTAMP))
            shutil.copy(flt, '.')

        # Create a readme listing origin of the prep files
        with open('README.txt', 'w') as readme:
            readme.write('{}\n'.format(time.ctime()))
            readme.write('FLT files copied from {}\n'.format(use_prep_path))


    # Do modeling
    logging.info(" ")
    logging.info("GroupFLT on field {}".format(field))
    grp = GroupFLT(grism_files=grism_files, direct_files=[], 
              ref_file=os.path.join(PATH_REF, p.ref_image),
              seg_file=os.path.join(PATH_REF, p.seg_map),
              catalog=os.path.join(PATH_REF, p.catalog),
              pad=p.pad,
              cpu_count=8)

    if not load_only:
        grp.compute_full_model(mag_limit=26)  # mag limit for contam model
        grp.refine_list(poly_order=2, mag_limits=[16, 24])

        # Save flt - model as PNGs.
        flt_residuals(grp)

        # Save grp.
        grp.save_full_data()

    return grp
        

#-------------------------------------------------------------------------------

@log_metadata
def fit(field='', mag_lim=35, use_model_path='.'):
    """
    Extract, stack, and fit (z and em lines).

    Parameters
    ----------
    mag_lim : int
        The magnitude limit of sources to extract and fit.
        By default '35', unrealistically high so to include everything.

    [will this create 1D and 2D traces for individual extractions, as well
    as the stack?]

    [do I really want the stacking and fitting to be in same function?]

    """
    #if use_model_path != '.':
    #   copy stuff?

    # don't need a new timestamp dir? Because always will be copying Extractions and
    # fits to the Extractions directory? (and finals to RELEASE?)

    # Load templates with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)

    # Load individual line templates for fitting the line fluxes
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)

    for id, mag in zip(np.array(grp.catalog['NUMBER']), np.array(grp.catalog['MAG_AUTO'])):
        if mag <= mag_lim:
            print(id, mag)
            # Extract the 2D traces
            beams = grp.get_beams(id, size=80) #size??
            if beams != []:
                mb = MultiBeam(beams, fcontam=1, group=field)

                # Save a FITS file with the 2D cutouts (beams) from the individual exposures
                mb.write_master_fits()

                # Fit polynomial model for initial continuum subtraction
                wave = np.linspace(2000,2.5e4,100)
                poly_templates = grizli.utils.polynomial_templates(wave, order=7)
                pfit = mb.template_at_z(
                    z=0, 
                    templates=poly_templates, 
                    fit_background=True, 
                    fitter='lstsq', 
                    get_uncertainties=2)

                # Drizzle grisms / PAs
                hdu, fig = mb.drizzle_grisms_and_PAs(
                    fcontam=0.2, 
                    flambda=False, 
                    kernel='point', 
                    size=32, 
                    zfit=pfit)

                # Save drizzled ("stacked") 2D trace as PNG and FITS
                fig.savefig('{0}_{1:05d}.stack.png'.format(field, id))
                hdu.writeto('{0}_{1:05d}.stack.fits'.format(field, id), clobber=True)

                # Fit the emission lines and redshifts
                #[this produces a file?]
                out = grizli.fitting.run_all(
                    id, 
                    t0=templ0, 
                    t1=templ1, 
                    fwhm=1200, 
                    zr=[0.5, 2.3], 
                    dz=[0.004, 0.0005], 
                    fitter='nnls', 
                    group_name=field,
                    prior=None, 
                    fcontam=0.,
                    pline=pline, 
                    mask_sn_limit=7, 
                    fit_beams=True, 
                    fit_stacks=False,  
                    root=field+'_',
                    fit_trace_shift=False, 
                    verbose=True, 
                    phot=None, 
                    scale_photometry=False, 
                    show_beams=True)

                mb, st, fit, tfit, line_hdu = out

                # make plots and save

#-------------------------------------------------------------------------------

@log_info
@log_metadata
def clear_grizli_pipeline(fields, ref_filter='F105W', mag_lim=25,
    do_steps=['prep', 'model', 'fit'], use_prep_path='.', use_model_path='.'):
    """ Main wrapper on pre-processing, modeling and extracting/fitting steps.

    Parameters
    ----------
    fields : list of strings
        The CLEAR fields; retain ability to specify the individual pointings
        so that can easily re-run single ones if find an issue.
    ref_filter : string
        The reference image's filter.
    mag_lim : int
        The magnitude limit of sources to extract and fit.
    do_steps : list of strings
        The pipeline steps to run.     

    """
    if use_prep_path != '.':
        use_prep_path = os.path.join(PATH_OUTPUTS, use_prep_path)
    if use_model_path != '.':
        use_model_path = os.path.join(PATH_OUTPUTS, use_model_path)

    # cd into outputs directory, as running grizli requires
    os.chdir(PATH_OUTPUTS_TIMESTAMP)
    #os.chdir(os.path.join(PATH_OUTPUTS, '2017.09.25.14.26.28'))
    logging.info("cd into {}".format(PATH_OUTPUTS_TIMESTAMP))

    # Find the files in RAW
    visits, filters = find_files(fields=fields)

    # In outputs run the prepsteps
    if 'prep' in do_steps:
        logging.info(" ")
        logging.info("PERFORMING PRE-PROCESSING STEP")
        logging.info("...")
        prep(visits=visits, ref_filt='F105W', ref_grism='G102')

    for field in fields:
        # Do the modeling; need have option which outputs subdir to use? (nominally, all will be same)
        if 'model' in do_steps:
            logging.info(" ")
            logging.info("PERFORMING MODELING STEP FOR FIELD {}"\
                .format(field.upper()))
            logging.info("...")
            grp = model(visits=visits, field=field, ref_filter=ref_filter, 
                mag_lim=mag_lim, use_prep_path=use_prep_path)
            # make this return a directory path?
            #print(grp)

        # Do the fitting; need have option which outputs subdir to use? 
        if 'fit' in do_steps:
            logging.info(" ")
            logging.info("PERFORMING EXTRACTING/STACKING/FITTING STEP FOR FIELD {}"\
                .format(field.upper()))
            logging.info("...")
            if 'model' not in do_steps:
                # fetch the grp
                # search in use_model_path for grp file?
                #grp = model(visits=visits, field=field, ref_filter=ref_filter, 
                #    mag_lim=mag_lim, use_prep_path='.', load_only=True)
            #fit()
        

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
    do_steps_help = "List the processing steps to run. Default is all three. Steps are NOT independent. " 
    do_steps_help += "If choose to run a step alone, be sure products exist from the previous step. "
    do_steps_help += "  prep "
    do_steps_help += "  model "
    do_steps_help += "  fit "
    rerun_help = "If the last run just crashed because you're testing things, "
    rerun_help += "re-run in the last outputs directory."
    rerun_help += "Do NOT do this for actual runs, to keep your results clean."
    prepdir_help = "Timestamp directory containing pre-processed files. '.' by default."
    modeldir_help = "Timestamp directory containing model files. '.' by default."
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--fields', dest = 'fields',
                        action = 'store', type = str, required = False,
                        help = fields_help, nargs='+', 
                        default=['GS1', 'GS2', 'GS3', 'GS4', 'GS5', 
                        'ERSPRIME', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])       
    parser.add_argument('--steps', dest = 'do_steps',
                        action = 'store', type = int, required = False,
                        help = do_steps_help,  nargs='+', 
                        default=['prep', 'model', 'fit'])    
    parser.add_argument('--mlim', dest = 'mag_lim',
                        action = 'store', required = False,
                        help = mag_lim_help, default=25)
    parser.add_argument('--ref', dest = 'ref_filter',
                        action = 'store', type = str, required = False,
                        help = ref_filter_help,  default='F105W')
    parser.add_argument('--rerun', dest = 'rerun',
                        action = 'store', type = str, required = False,
                        help = rerun_help,  default='False')
    parser.add_argument('--prepdir', dest = 'prepdir',
                        action = 'store', type = str, required = False,
                        help = prepdir_help,  default='.')
    parser.add_argument('--modeldir', dest = 'modeldir',
                        action = 'store', type = str, required = False,
                        help = modeldir_help,  default='.')

    args = parser.parse_args()
     
    return args

#------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------- 

if __name__=='__main__':
    # Parse user inputs
    # (not yet implemented)
    args = parse_args()
    fields = args.fields
    do_steps = args.do_steps
    mag_lim = args.mag_lim
    ref_filter = args.ref_filter
    rerun = tobool(args.rerun)
    prepdir = args.prepdir
    modeldir = args.modeldir

    if rerun:
        PATH_OUTPUTS_TIMESTAMP = retrieve_latest_outputs(path_outputs=PATH_OUTPUTS)
    else:
        # Create the output directory
        PATH_OUTPUTS_TIMESTAMP = store_outputs(path_outputs=PATH_OUTPUTS, 
           store_type='')

    # Setup logging to save in the output directory
    setup_logging(__file__, path_logs=PATH_OUTPUTS_TIMESTAMP, stdout=True)

    #meta_record_imports(__file__, print_or_log='print') #but how direct this to a log file?

    # Call the main pipeline function.
    clear_grizli_pipeline(fields=['GN2'], ref_filter='F105W', mag_lim=25, 
        do_steps=['model'], prepdir='2017.09.22.15.30.29', modeldir=modeldir)

