""" Module containing all wrapper functions to run CLEAR interlacing and
extractions.

Add catalogs or overlapping fields to the global lists and dictionaries 
up top.  

Note what I call a 'field' here means a 'pointing' of either GOODS-N or
GOODS-S, for example, GS1 or ERSPRIME.

Use:
    
    Be in the outputs directory.  In the Main, change the field, ref_image,
    and step numbers.  Then,

    >>> python clear_pipeline.py --fields (required) --steps (optional) --mlim (optional) --ref (optional)

    --fields : List the fields over which to run pipeline. Default is all. For example,

    --steps : List the processing steps to run. Default is all five. Steps are NOT independent.
              If choose to run a step alone, be sure products exist from the previous step.
              1 - Interlace visits
              2 - Create contamination models
              3 - Extract traces
              4 - Stack traces
              5 - Fit redshifts and emission lines of traces
    
    --mlim : The magnitude limit for extraction. Default is 25.

    --ref : The reference image filter. Choose either F105W or F125W. Default is F105W.

Example:

    >>> python clear_pipeline.py --fields GS1 ERSPRIME --steps 1 2 

Authors:
    
    C.M. Gosmeyer, 2016

Future improvements:

    1. Add command-line option to extract either based on magnitude OR on catalog
    2. Logging would awesome, but probably not worth effort at this point.


Here are the 14227 visits contained in each field. 

ERSPRIME: 19, 20, 21, 22, 23
GS1: 07, 08, 09, 10, 11, 12
GS2: 01, 02, 03, 04, 05, 06
GS3: 30, 31, 32, 33, 34, 35
GS4: 24, 25, 26, 27, 28, 29
GS5: 13, 14, 15, 16, 17, 18
GN1: 46, 47, 48, 49, 50
GN2: 51, 52, 53, 54 (A4), 55
GN3: 56, 57, 58, 59, 60
GN4: 61, 62, 63, 64, 65
GN5: 41, 42, 43, 44, 45
GN7: 36, 37, 38, 39, 40

GN6 is a myth
"""

from __future__ import print_function

import argparse
import astropy.io.fits as pyfits
import glob
import numpy as np
import os
import shutil
import time

import threedhst
import unicorn

from unicorn import interlace_test

from astropy.table import Table
from astropy.io import ascii

from find_pointing_start import find_pointing_start
from set_paths import paths

# Define catalogs for S and N.
# These are only needed to extract subsets of source IDs.
quiescent_cats = {'N' : ['UVJ_quiescent_goodsn.dat'], 
                  'S' : ['UVJ_quiescent_goodss.dat'], 
                  'name' : ['quiescent']}
#emitters_cats = {'N' : ['Steves_source_goodsn_w_ids.dat'], 
#                 'S' : ['Steves_source_goodss_w_ids.dat'], 
#                 'name' : ['emitters']}
emitters_cats = {'N' : ['GOODSN_plus_translation_all.dat'], 
                 'S' : ['GOODSS_plus_translation_all.dat'], 
                 'name' : ['emitters']}
ivas_cat = {'N' : ['Ivas_goodsn.dat'], 
            'name' : ['ivas']}
zn_cats = {'N' : ['added_sources_N_key_z{}.dat'.format(str(s)) for s in [3,4,5,6,7,8]], 
             'S' : ['added_sources_S_key_z{}.dat'.format(str(s)) for s in [3,4,5,6,7,8]],
             'name' : ['z{}'.format(str(s)) for s in [3,4,5,6,7,8]]}

sijie_cats = {'N':['noExtractionListNorth.txt'], 'name':['sijiegoodsn'],
              'S':['NoExtracationGoodsList.txt'], 'name':['sijiegoodss']}

# Full catalogs, to be used if extracting based on magnitude. 
full_cats = {'N' : ['GoodsN_plus.cat'], 
                 'S' : ['GoodsS_plus.cat'], 
                 'name' : ['full']}

all_cats = {'emitters':emitters_cats, 'full':full_cats, 
            'quiescent':quiescent_cats, 
            'sijie':sijie_cats, 'zn':zn_cats, }

## Associate CLEAR Goods-N pointings with overlapping Barro pointings.
# To be extra confusing, each Barro 'visit' is also a single 'pointing'.
# But the dictionary below maps them to a single CLEAR pointing, so that they
# can be stacked appropriately.
overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11']}


#-------------------------------------------------------------------------------  

def interlace_clear(field, ref_filter):
    """ Interlaces given field with the given reference image filter.

    Parameters
    ----------
    field : string
        The GOODS field to process. 
    ref_filter : string
        Filter of the reference image.

    ** Step 1 of Interlace steps. **
    ** Takes a field and loops over all pointings (visits) **

    We interlace instead of drizzle.
    Drizzling with just a few images produces correlated noise and loss of 
    resoluation as the uneven output grid gets smoothed over.
    Correlated noise can mimic emission or absorption features.
    Interlacing improves the sampling of the PSF by a factor of two without 
    having to interpolate. 
    Can interlace these images because (a) the relative pointing errors 
    between images of a given set are small, 0.1 pixels and (b) the dithers 
    between images are small, under 10 pixels and the relative distortion at 
    that scale is small.
    The primary drawback is that if one of the four images is missing there 
    will be a hole (empty pixels) in the final interlaced image.
    A second drawback is that only observations taken at same rotation angle 
    can be combined.

    Produces
    --------
    - *inter.fits
    - *inter.reg
    - *inter_seg.fits
    - *ref_inter.fits
    - *radec.dat

    Checks
    ------
    - Make sure *G102_ref_inter.fits and *inter_seg.fits are all aligned!
    - Load each *inter.reg into each *G102_ref_inter.fits to check 
      objects overlap catalog.
    - These checks must be done by Visit, not field, since Visits are rotated.
      So do for example,

    ds9 GN7-38-315-F105W_inter.fits GN7-38-315-G102_ref_inter.fits GN7-38-315-G102_inter_seg.fits &


    """

    print("Processing field {}".format(field))

    from unicorn.reduce import adriz_blot_from_reference as adriz_blot

    path_to_REF = paths['path_to_ref_files'] + 'REF/'
      
    NGROWX = 200
    NGROWY = 30
    if 'GDN' in field:
        pad = 500
    else:
        pad = 60

    if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
        REF_EXT = 0
        SEG_IMAGE = os.path.join(path_to_REF, 'Goods_S_plus_seg.fits') #'goodss_3dhst.v4.0.F160W_seg.fits'
        if ref_filter == 'F125W':
            print("Using ref F125W!")
            CATALOG = os.path.join(path_to_REF, 'GoodsS_plus_merged.cat') #'goodss_3dhst.v4.0.F125W_conv_fix.cat'
            REF_IMAGE = os.path.join(path_to_REF, 'goodss_3dhst.v4.0.F125W_orig_sci.fits')
            REF_FILTER = 'F125W'
        elif ref_filter == 'F105W':
            print("Using ref F105W!")
            CATALOG = os.path.join(path_to_REF, 'goodss-F105W-astrodrizzle-v4.3_drz_sub_plus.cat')
            # The new ref image. 
            REF_IMAGE = os.path.join(path_to_REF, 'goodss-F105W-astrodrizzle-v4.3_drz_sci.fits')
            REF_FILTER = 'F105W'

    elif 'GN' in field or 'GDN' in field:
        REF_EXT = 0
        SEG_IMAGE = os.path.join(path_to_REF, 'Goods_N_plus_seg.fits') #'goodsn_3dhst.v4.0.F160W_seg.fits'
        if ref_filter == 'F125W':
            print("Using ref F125W!")
            CATALOG = os.path.join(path_to_REF, 'GoodsN_plus_merged.cat') #'goodsn_3dhst.v4.0.F125W_conv.cat'
            REF_IMAGE = os.path.join(path_to_REF, 'goodsn_3dhst.v4.0.F125W_orig_sci.fits')
            REF_FILTER = 'F125W'
        elif ref_filter == 'F105W':
            print("Using ref F105W!")
            CATALOG = os.path.join(path_to_REF, 'goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat')
            # The new ref image.
            REF_IMAGE = os.path.join(path_to_REF, 'goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits')
            REF_FILTER = 'F105W'


    grism = glob.glob(field+'*G102_asn.fits')
    print("grism: {}".format(grism))

    for i in range(len(grism)):
        pointing=grism[i].split('_asn')[0]
        print(pointing)
        
        # Find whether pointing begins with a direct image (0) or grism (1).
        ref_exp = 0 #find_pointing_start(pointing)
        print("ref_exp: {}, pointing: {}".format(ref_exp, pointing))

        adriz_blot(
            pointing=pointing, 
            pad=pad, 
            NGROWX=NGROWX, 
            NGROWY=NGROWY, 
            growx=2, 
            growy=2, 
            auto_offsets=True, 
            ref_exp=ref_exp, 
            ref_image=REF_IMAGE, 
            ref_ext=REF_EXT, 
            ref_filter=REF_FILTER, 
            seg_image=SEG_IMAGE, 
            cat_file=CATALOG, 
            grism='G102')     
       
        if 'GN5-42-346' in pointing:
             # These stare images had a bad dither. Set to 1x1 binning.
            print("Binning 1x1!")
            growx=1
            growy=1
        else:
            growx=2
            growy=2 

        # Interlace the direct images. Set ref_exp to always be zero.                                                                             
        unicorn.reduce.interlace_combine(
            pointing.replace('G102','F105W'), 
            view=False, 
            use_error=True, 
            make_undistorted=False,
            pad=pad, 
            NGROWX=NGROWX, 
            NGROWY=NGROWY,
            ddx=0, 
            ddy=0, 
            growx=growx, 
            growy=growy, 
            auto_offsets=True, 
            ref_exp=0)
        # Interlace the grism images.
        unicorn.reduce.interlace_combine(
            pointing, 
            view=False, 
            use_error=True, 
            make_undistorted=False, 
            pad=pad, 
            NGROWX=NGROWX, 
            NGROWY=NGROWY, 
            ddx=0, 
            ddy=0, 
            growx=2, 
            growy=2, 
            auto_offsets=True, 
            ref_exp=ref_exp)


    print("*** interlace_clear step complete ***")


#-------------------------------------------------------------------------------  

def model_clear(field, mag_lim=None):
    """ Creates model contam images. 

    ** Step 2. of Interlace steps. **

    Parameters
    ----------
    field : string
        The GOODS field to process. 
    mag_lim : int 
        The magnitude down from which to extract.  If 'None' ignores
        magnitude filter.

    Produces
    --------
    - *inter_model.pkl
    - *inter_model.fits
    - *inter_0th.reg
    - *maskbg.dat
    - *maskbg.png

    Checks
    ------
    - Contam data is saved in *pkl and *inter_model.fits files.
    the latter are meant to look like the observational *inter.fits.
    - Display the *inter_model.fits and *inter.fits.
    Set scale the same.
    The spectra should match; there shouldn't be extra objects;
    the start-end of spectra should be the same; images should
    be of same brightness.

    For example,

    ds9 GDN7-07-335-G102_inter_model.fits GDN7-07-335-G102_inter.fits &


    ** Note that the *mask* files do not overwrite **
    ** Delete these files first if want to do a rerun **

    """
    grism_asn = glob.glob(field+'*G102_asn.fits')
    
    if mag_lim == None or mag_lim < 24:
        contam_mag_lim = 24
    else:       
        contam_mag_lim = mag_lim

    for i in range(len(grism_asn)):
        root = grism_asn[i].split('-G102')[0]
        direct = 'F105W'
        grism = 'G102'
        m0 = unicorn.reduce.GrismModel(
            root=root,
            direct=direct,
            grism=grism)
        model_list = m0.get_eazy_templates(
            dr_min=0.5, 
            MAG_LIMIT=25)
        model = unicorn.reduce.process_GrismModel(
            root=root, 
            model_list=model_list,
            grow_factor=2, 
            growx=2, 
            growy=2, 
            MAG_LIMIT=contam_mag_lim, 
            REFINE_MAG_LIMIT=21, 
            make_zeroth_model=False, 
            use_segm=False, 
            model_slope=0, 
            direct=direct, 
            grism=grism, 
            BEAMS=['A', 'B', 'C', 'D','E'], 
            align_reference=False)
        if not os.path.exists(os.path.basename(model.root) + '-G102_maskbg.dat'):
            model.refine_mask_background(
                threshold=0.002, 
                grow_mask=14, 
                update=True, 
                resid_threshold=4, 
                clip_left=640, 
                save_figure=True, 
                interlace=True)

    print("*** model_clear step complete ***")


#-------------------------------------------------------------------------------  

def extract_clear(field, tab, mag_lim=None):
    """ Extracts all sources given in tab from the given field.

    ** Step 3. of Interlace steps. **

    Parameters
    ----------
    field : string
        The GOODS field to process. 
    tab : dictionary
        Values for each source read from catalog.
    mag_lim : int 
        The magnitude down from which to extract.  If 'None' ignores
        magnitude filter.

    Produces
    --------
    - <field>-<visit>-<orient>-G102_<id>.1D.fits (ie, GN3-56-175-G102_28015.1D.fits)
    - <field>-<visit>-<orient>-G102_<id>.1D.png
    - <field>-<visit>-<orient>-G102_<id>.2D.fits 
    - <field>-<visit>-<orient>-G102_<id>.2D.png

    Checks
    ------
        *2D.png should show some extracted spectra (most may be empty)

    """  

    grism = glob.glob(field+'*G102_asn.fits')[:1]

    # Keep magnitude limit for contam models from being too low. 
    if mag_lim == None or mag_lim < 24:
        contam_mag_lim = 24
    else:
        contam_mag_lim = mag_lim

    for i in range(len(grism)):
        root = grism[i].split('-G102')[0]
        # Takes a subset of the full catalog, retaining
        # only objects visible in the field (model.cat.id)
        model, ids = return_model_and_ids(
            root=root, contam_mag_lim=contam_mag_lim, tab=tab)

        #print(model.cat.mag)
        #print(model.cat.id)

        for id, mag in zip(model.cat.id, model.cat.mag):
            # If extracting by magnitude limit.
            print("magnitude limit is {}".format(mag_lim))
            if mag_lim != None:
                print("extracting down to magnitude limit {}".format(mag_lim))
                if (id in ids and mag <= mag_lim):   
                    print("id, mag: ", id, mag)
                    try:
                        # In spite of name, also creates 1D FITS.
                        model.twod_spectrum(
                            id=id, 
                            grow=1, 
                            miny=-36, 
                            maxy=None, 
                            CONTAMINATING_MAGLIMIT=contam_mag_lim, 
                            refine=False, 
                            verbose=False, 
                            force_refine_nearby=False, 
                            USE_REFERENCE_THUMB=True,
                            USE_FLUX_RADIUS_SCALE=3, 
                            BIG_THUMB=False, 
                            extract_1d=True)
                        model.show_2d(savePNG=True, verbose=True)
                        print("Extracted {}".format(id))
                    except:
                        continue

            # If extracting by catalog.
            elif mag_lim == None:
                print("magnitude limit is undefined. extracting all sources")
                if (id in ids):
                    print("id, mag: ", id, mag)
                    try:
                        # In spite of name, also creates 1D FITS.
                        model.twod_spectrum(
                            id=id, 
                            grow=1, 
                            miny=-36, 
                            maxy=None, 
                            CONTAMINATING_MAGLIMIT=contam_mag_lim, 
                            refine=False, 
                            verbose=False, 
                            force_refine_nearby=False, 
                            USE_REFERENCE_THUMB=True,
                            USE_FLUX_RADIUS_SCALE=3, 
                            BIG_THUMB=False, 
                            extract_1d=True)
                        model.show_2d(savePNG=True, verbose=True)
                        print("Extracted {}".format(id))
                    except:
                        continue

    print("*** extract_clear step complete ***")        


#-------------------------------------------------------------------------------  

def stack_clear(field, tab, catname, ref_filter, mag_lim=None):
    """ Stacks the extractions for all the sources in the given field.

    Parameters
    ----------
    field : string
        The GOODS field to process. Needs to know to reference GN or GS 
        catalog.
    tab : dictionary
        Values for each source from the catalog.
    catname : string
        Name of the catalog, for naming output directory.
    ref_filter : string
        Filter of the reference image.
    mag_lim : int 
        The magnitude down from which to extract.  If 'None' ignores
        magnitude filter.

    Produces
    --------
    - <field>-G102_<id>.1D.fits (ie, GN3-G102_28121.1D.fits)
    - <field>-G102_<id>.1D.png
    - <field>-G102_<id>.2D.fits
    - <field>-G102_<id>.2D.png

    Checks
    ------
        In *stack.png check that the contam models reasonably match actual
        image and that the contam-flux is reasonable.
        python> !open *stack.png

    Notes
    -----
        This should stack ALL the *2D* files present in the directory 
        that have the same id, REGARDLESS of the field name.
    """

    from unicorn.hudf import Stack2D

    grism = glob.glob(field+'*G102_asn.fits')

    # Keep magnitude limit for contam models from being too low. 
    if mag_lim == None or mag_lim < 24:
        contam_mag_lim = 24
    else:
        contam_mag_lim = mag_lim

    for i in range(len(grism)):
        root = grism[i].split('-G102')[0]
        model, ids = return_model_and_ids(
            root=root, contam_mag_lim=contam_mag_lim, tab=tab)

        for id, mag in zip(model.cat.id, model.cat.mag):
            # If extracting by magnitude limit.
            if mag_lim != None:
                if (id in ids and mag <= mag_lim):   
                    print("id, mag: ", id, mag)
                    try:
                        search='*-*-*-G102'
                        print('searching %s*%05d.2D.fits'%(search, id))
                        spec = Stack2D(
                            id=np.int(id), 
                            inverse=False, 
                            scale=[1,99], 
                            fcontam=2.,
                            ref_wave = 1.05e4,
                            root='{}-G102'.format(field), 
                            search='*-*-*-G102', 
                            files=None, 
                            go=True, 
                            new_contam=False)
                    except:
                        continue

            # If extracting by catalog.
            elif mag_lim == None:
                if (id in ids):
                    print("id: ", id)
                    try:
                        search='*-*-*-G102'
                        print('searching %s*%05d.2D.fits'%(search, id))
                        spec = Stack2D(
                            id=np.int(id), 
                            inverse=False, 
                            scale=[1,99], 
                            fcontam=2.,
                            ref_wave = 1.05e4,
                            root='{}-G102'.format(field), 
                            search='*-*-*-G102', 
                            files=None, 
                            go=True, 
                            new_contam=False)
                    except:
                        continue


    #cleanup_extractions(field=field, cat=cat, catname=catname, ref_filter=ref_filter)

    print("*** stack_clear step complete ***")
    

#------------------------------------------------------------------------------- 

def fit_redshifts_and_emissionlines(field, tab, mag_lim=None):
    """ Fits redshifts and emission lines. 

    Parameters
    ----------
    field : string
        The GOODS field to process. Needs to know to reference GN or GS 
        catalog.
    tab : dictionary
        Values for each source from the catalog.
    mag_lim : int 
        The magnitude down from which to extract.  If 'None' ignores
        magnitude filter.


    Notes
    -----
        Based on lines 182-207 of unicorn/aws.py

    """
    # add way to copy all 1D and 2D files from latest extraction for field? 

    # Keep magnitude limit for contam models from being too low. 
    if mag_lim == None or mag_lim < 24:
        contam_mag_lim = 24
    else:
        contam_mag_lim = mag_lim

    # Find the unique root for the pointing, for its stacked 2D.fits.
    twods = glob.glob('{}-G102_*.2D.fits'.format(field))
    print("twods: {}".format(twods))

    if len(twods) >= 1:
        pointing_root = np.unique(np.array([twod.split('_')[0] for twod in twods]))[0]

        # Get IDs of all the sources stacked.
        stacked_ids = np.array([twod.split('_')[1].split('.2D.fits')[0] for twod in twods], dtype=int)
        print("stacked_ids: ")
        print(stacked_ids)

        for id in stacked_ids:
            obj_root='{}_{:05d}'.format(pointing_root, id)
            print("obj_root: ", obj_root)
            try:
                # Redshift fit
                gris = interlace_test.SimultaneousFit(
                    obj_root,
                    lowz_thresh=0.01, 
                    FIGURE_FORMAT='png')     
            except (ValueError) as err: 
                print(err)
                print("Error in {}; skipping...".format(obj_root))
                continue

            #if gris.status is False:
            #    continue

            if not os.path.exists(obj_root + '.new_zfit.pz.fits'):
                print("Fitting z...")
                try:
                    gris.new_fit_constrained()
                    gris.new_save_results()
                    gris.make_2d_model()
                except:
                    continue
            if not os.path.exists(obj_root+'.linefit.fits'):
                print("Fitting em lines...")
                #try:
                    # Emission line fit
                gris.new_fit_free_emlines(ztry=None, NSTEP=600)
                #except:
                #    continue
    else:
        print("Skipping {}-G102_*.2D.fits because none found".format(field))


    print("*** fit redshifts and emission lines step complete ***")


#-------------------------------------------------------------------------------

def return_model_and_ids(root, contam_mag_lim, tab):
    """ Returns the model object and list of ids of sources present in field.

    Parameters
    ----------
    root : string
        Root string <field>-<visit>-<orient>.
    contam_mag_lim : int
        The mag limit to bring contamination model.
    tab : dictionary
        Values for each source from the catalog.

    Returns
    -------
    model : GrismModel
        The GrismModel object for interlaced field.
    ids : list
        List of ids of sources present in field.

    """
    model = unicorn.reduce.process_GrismModel(
        root=root, 
        grow_factor=2, 
        growx=2, 
        growy=2, 
        MAG_LIMIT=contam_mag_lim,
        REFINE_MAG_LIMIT=21, 
        make_zeroth_model=False, 
        use_segm=False,
        model_slope=0, 
        direct='F105W', 
        grism='G102', 
        BEAMS=['A', 'B', 'C', 'D','E'],
        align_reference=False)

    try:
        ids = tab['ID']
    except:
        try:
            # Needed for the 'plux_z*' catalogs.
            ids_raw = tab['id']
            ids = []
            for id in ids_raw:
                ids.append(id.split('_')[-1])
        except:
            # For 'Goods*plus' catalogs.
            ids = tab['NUMBER']

    return model, ids


#-------------------------------------------------------------------------------  

def sort_outputs(field, overlapping_field, catname, ref_filter, mag_lim=None):
    """ Sorts outputs into following tree of directories in Extractions.
    The primary branch of each (the pointing/field) will be tarred. 

    POINTING_REFFILTER_DATE
        ORIENT1
            ...
        ORIENT2
            1D
                FITS
                PNG
            2D
                FITS
                PNG
        ...
        COMBINED
            1D
                FITS
            2D
                FITS
                PNG
            LINEFIT
                DAT
                FITS
                PNG
                CHAIN_PNG
            ZFIT
                DAT
                FITS
                PNG
                PZ_FITS
                2D_PNG
                TILT_DAT
                TILT_PNG

    Parameters
    ----------
    field : string
        The GOODS field to process. Needs to know to reference GN or GS 
        catalog.
    overlapping_field : string
        Should be None, unless the field is from 3DHST, so that the 1D and 2D
        files can be appropriately sorted. 
    catname : string
        Name of the catalog, for naming output directory.
    ref_filter : string
        Filter of the reference image.
    mag_lim : int 
        The magnitude down from which to extract.  If 'None' ignores
        magnitude filter.

    """ 

    path_to_Extractions = paths['path_to_Extractions']

    # Create the top-level name for the extractions.
    # Name by catalog and magnitude limit. 
    basename = 'cat{}_maglim{}'.format(catname.upper(), str(mag_lim)) 

    topdir = os.path.join(path_to_Extractions, field, '{}_ref{}_{}'.format(basename, ref_filter,
        time.strftime('%Y-%m-%d')))

    if not os.path.isdir(topdir):
        os.mkdir(topdir)
  

    if overlapping_field == None:
        # Normal case.
        # Glob all the ORIENT-specific 1D FITS files.
        orient_1D = glob.glob('{}-[0-9]*1D.fits'.format(field.upper())) + \
            glob.glob('{}-[0-9]*1D.png'.format(field.upper()))
        # Glob all the ORIENT-specific 2D FITS and PNG files.
        orient_2D = glob.glob('{}-[0-9]*2D.fits'.format(field.upper())) + \
            glob.glob('{}-[0-9]*2D.png'.format(field.upper()))

    else:
        # Barro field
        # Glob all the ORIENT-specific 1D FITS files.
        orient_1D = glob.glob('{}-[0-9]*1D.fits'.format(overlapping_field.upper())) + \
            glob.glob('{}-[0-9]*1D.png'.format(overlapping_field.upper()))
        # Glob all the ORIENT-specific 2D FITS and PNG files.
        orient_2D = glob.glob('{}-[0-9]*2D.fits'.format(overlapping_field.upper())) + \
            glob.glob('{}-[0-9]*2D.png'.format(overlapping_field.upper()))

    # If field is GN2, need make a special case so can catch visit A4.
    if field.upper() == 'GN2':
        orient_1D += glob.glob('{}-A4*1D.fits'.format(field.upper())) + \
            glob.glob('{}-A4*1D.png'.format(field.upper()))
        orient_2D += glob.glob('{}-A4*2D.fits'.format(field.upper())) + \
            glob.glob('{}-A4*2D.png'.format(field.upper()))


    # These files need be sorted and moved by visit, 'bottomdir', and extension.
    check_and_create_dirs(file_list=orient_1D, topdir=topdir, 
        bottomname='1D', orient=True, overwrite_ext=None)
    check_and_create_dirs(file_list=orient_2D, topdir=topdir, 
        bottomname='2D', orient=True, overwrite_ext=None)

    # This finishes sorting the ORIENTS.
    # Next on to the COMBINED, which are more numerous and varied. 
    comb_1D = glob.glob('{}-G102*1D.fits'.format(field.upper())) + \
        glob.glob('{}-G102*1D.png'.format(field.upper()))
    comb_2D = glob.glob('{}-G102*2D.fits'.format(field.upper())) + \
        glob.glob('{}-G102*stack.png'.format(field.upper()))
    comb_linefit = glob.glob('{}-G102*linefit.dat'.format(field.upper())) + \
        glob.glob('{}-G102*linefit.fits'.format(field.upper())) + \
        glob.glob('{}-G102*linefit.png'.format(field.upper()))
    comb_linefit_chain = glob.glob('{}-G102*linefit.chain.png'.format(field.upper()))
    comb_zfit = glob.glob('{}-G102*zfit.dat'.format(field.upper())) + \
        glob.glob('{}-G102*zfit.fits'.format(field.upper())) + \
        glob.glob('{}-G102*zfit.png'.format(field.upper()))
    comb_zfit_pz = glob.glob('{}-G102*zfit.pz.fits'.format(field.upper())) 
    comb_zfit_2D = glob.glob('{}-G102*zfit.2D.png'.format(field.upper())) 
    comb_zfit_tilt_dat = glob.glob('{}-G102*zfit_tilt.dat'.format(field.upper())) 
    comb_zfit_tilt_png = glob.glob('{}-G102*zfit_tilt.png'.format(field.upper())) 

    # Sort and move files. 
    check_and_create_dirs(file_list=comb_1D, topdir=topdir, 
        bottomname='1D', orient=False, overwrite_ext=None)
    check_and_create_dirs(file_list=comb_2D, topdir=topdir, 
        bottomname='2D', orient=False, overwrite_ext=None)
    check_and_create_dirs(file_list=comb_linefit, topdir=topdir, 
        bottomname='LINEFIT', orient=False, overwrite_ext=None)
    check_and_create_dirs(file_list=comb_linefit_chain, topdir=topdir, 
        bottomname='LINEFIT', orient=False, overwrite_ext='CHAIN_PNG')
    check_and_create_dirs(file_list=comb_zfit, topdir=topdir, 
        bottomname='ZFIT', orient=False, overwrite_ext=None)
    check_and_create_dirs(file_list=comb_zfit_pz, topdir=topdir, 
        bottomname='ZFIT', orient=False, overwrite_ext='PZ_FITS')
    check_and_create_dirs(file_list=comb_zfit_2D, topdir=topdir, 
        bottomname='ZFIT', orient=False, overwrite_ext='2D_PNG')
    check_and_create_dirs(file_list=comb_zfit_tilt_dat, topdir=topdir, 
        bottomname='ZFIT', orient=False, overwrite_ext='TILT_DAT')
    check_and_create_dirs(file_list=comb_zfit_tilt_png, topdir=topdir, 
        bottomname='ZFIT', orient=False, overwrite_ext='TILT_PNG')

    # Now tar the top-level directory.
    shutil.make_archive(os.path.join(path_to_Extractions, field, 
        '{}_{}_ref{}'.format(field, basename, ref_filter)), 
        'gztar', topdir)


#-------------------------------------------------------------------------------  

def check_and_create_dirs(file_list, topdir, bottomname, orient=True, 
    overwrite_ext=None):
    """ Creates levels of directories, sorts, and moves files.

    Parameters
    ----------
    file_list : list of strings
        List of files to sort and move.
    topdir : string
        The top-most path+directory. 
    bottomdir : string
        Name of the bottom directory. E.g., '1D', 'ZFIT', etc.
    orient : {True, False}
        Keep True if the directory between top and bottom is to be
        named for a visit. Turn False if that directory is to be
        'COMBINED'.
    overwrite_ext : string
        Keep 'None' if you want to directory under 'bottomdir' to be named
        after the file's extension. E.g., 'FITS', 'PNG', 'DAT'.
        Otherwise, enter the string you would like this directory to be
        named. 
    """
    for f in file_list:    
        visit = f.split('-')[1]
        program = f.split('-')[0]

        print("program: ", program)

        if orient:
            # Create directory for that visit, if it doesn't exist.
            # Will need special case for images from G. Barro's program 13420.
            if 'GDN' in program:
                visitdir = os.path.join(topdir, 'barro-{}'.format(visit))
            else:
                visitdir = os.path.join(topdir, 'clear-{}'.format(visit))
        else:
            # Assume the directory with be COMBINED
            visitdir = os.path.join(topdir, 'COMBINED')

        if not os.path.isdir(visitdir):
            os.mkdir(visitdir)

        # Now check whether the bottomdir directory exists.
        bottomdir = os.path.join(visitdir, bottomname)
        if not os.path.isdir(bottomdir):
            os.mkdir(bottomdir)

        # Finally! Sort by the extension.
        if overwrite_ext == None:
            ext = (f.split('.')[-1]).upper()
        else:
            ext = overwrite_ext

        extdir = os.path.join(bottomdir, ext)
        if not os.path.isdir(extdir):
            os.mkdir(extdir)

        print("Moving {} to {}".format(f, extdir))
        shutil.move(f, os.path.join(extdir, f))

    

#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------  

def clear_pipeline_main(fields, do_steps, cats, mag_lim, ref_filter):
    """ Main for the interlacing and extracting step. 

    Parameters
    ----------
    fields : list of strings
        The GOODS fields to process. Needs to know to reference GN or GS 
        catalog.
    do_steps : list of ints 
        The step numbers to perform.
        1 - Interlace visits
        2 - Create contam models
        3 - Extract traces
        4 - Stack traces
        5 - Fit redshifts and emission lines of traces
    cats : dictionary
        Dictionaries of the field (N/S) and catalog name.
        To extract all mags to this limit from full catalog, select 
        'full_cats'.
    mag_lim : int 
        The magnitude down from which to extract. If 'None', Then
        defaults to 24. 
    ref_filter : string
        Filter of the reference image.    
    """
    path_to_REF = paths['path_to_ref_files'] + 'REF/'

    for field in fields:
        print("***Beginning field {}***".format(field))
        print("")

        # Choose the field for the catalogs, where the catalog options include emitters and quiescent.
        if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
            cats_field = cats['S']
        elif 'GN' in field or 'GDN' in field:
            cats_field = cats['N']

        for cat, catname in zip(cats_field, cats['name']):
            print("***Beginning catalog {}***".format(cat))
            print("")

            if 'GN' in field:
                # add primary CLEAR pointing to fields
                overlapping_fields_all = overlapping_fields[field]
                overlapping_fields_all.append(field)
                for overlapping_field in overlapping_fields_all:
                    # Remember that a Barro 'field' really is just a visit of
                    # the CLEAR field! Buuuut steps 1, 2, & 3 treat them as full-fledged
                    # fields. 
                    print("***Beginning overlapping 13420 field {}***".format(overlapping_field))
                    print("")
                    if 1 in do_steps:
                        interlace_clear(field=overlapping_field, ref_filter=ref_filter)
                    if 2 in do_steps:
                        model_clear(field=overlapping_field, mag_lim=mag_lim)
                    if 3 in do_steps:
                        tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                        extract_clear(field=overlapping_field, tab=tab, mag_lim=mag_lim)
                if 4 in do_steps:
                    tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                    stack_clear(field=field, tab=tab, catname=catname, ref_filter=ref_filter, mag_lim=mag_lim) 
                if 5 in do_steps:
                    print("Starting z and em line fitting!")
                    tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                    fit_redshifts_and_emissionlines(field=field, tab=tab)

                    for overlapping_field in overlapping_fields_all:  
                        ## remember to remove this after done making 1D.PNG files!
                        #onedfiles = glob.glob('{}*1D.fits'.format(overlapping_field.upper()))
                        #for onedfile in onedfiles:
                        #    unicorn.reduce.Interlace1D(file=onedfile, PNG=True) 
                        ##
                        sort_outputs(field=field, overlapping_field=overlapping_field, 
                            catname=catname, ref_filter=ref_filter, mag_lim=mag_lim)

            else:
                # GS fields
                if 1 in do_steps:
                    interlace_clear(field=field, ref_filter=ref_filter)
                if 2 in do_steps:
                    model_clear(field=field, mag_lim=mag_lim)

                if 3 in do_steps:
                    tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                    extract_clear(field=field, tab=tab, mag_lim=mag_lim)              
                if 4 in do_steps:
                    tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                    stack_clear(field=field, tab=tab, catname=catname, ref_filter=ref_filter, mag_lim=mag_lim) 
                if 5 in do_steps:
                    tab = Table.read(os.path.join(path_to_REF, cat), format='ascii')
                    fit_redshifts_and_emissionlines(field=field, tab=tab)
                    sort_outputs(field=field, overlapping_field=None, catname=catname, 
                        ref_filter=ref_filter, mag_lim=mag_lim)

            if 4 in do_steps and 5 not in do_steps:
                print("")
                sort_outputs(field=field, overlapping_field=None, catname=catname, 
                    ref_filter=ref_filter, mag_lim=mag_lim)           
 

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


if __name__=='__main__':

    args = parse_args()
    fields = args.fields
    do_steps = args.do_steps
    mag_lim = args.mag_lim
    ref_filter = args.ref_filter
    cat_names = args.cats

    print("CLEAR pipeline running on fields {},\nover steps {},\nat mag-limit {},\nfor reference filter {}\nover catalogs {}.\n"\
        .format(fields, do_steps, mag_lim, ref_filter, cat_names))

    ## switching between catalog and mag-limit extraction should be made a command-line option
    # all_cats = [quiescent_cats, emitters_cats, ivas_cat, zn_cats]
    # mag_lim = None

    # match to the catalog dictionary
    cats_list = []
    for cat in cat_names:
        cats_list.append(all_cats[cat])
    print(cats_list)

    # From argparse, None is a string.
    if mag_lim.lower() == 'none':
        mag_lim = None

    for cat in cats_list: #[full_cats]:
        clear_pipeline_main(
            fields=fields, 
            do_steps=do_steps, 
            cats=cat, 
            mag_lim=mag_lim, 
            ref_filter=ref_filter)

