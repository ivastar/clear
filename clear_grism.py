import os
import glob
import numpy as np
import astropy.io.fits as pyfits
import shutil
import time
import zipfile

import threedhst
import unicorn

from astropy.table import Table
from astropy.io import ascii

from find_pointing_start import find_pointing_start
from set_paths import paths

# Define catalogs for S and N.
quiescent_cats = {'N' : ['UVJ_quiescent_goodsn.dat'], 
                  'S' : ['UVJ_quiescent_goodss.dat'], 
                  'name' : ['quiescent']}
emitters_cats = {'N' : ['Steves_source_goodsn_w_ids.dat'], 
                 'S' : ['Steves_source_goodss_w_ids.dat'], 
                 'name' : ['emitters']}
ivas_cat = {'N' : ['Ivas_goodsn.dat'], 
            'name' : ['ivas']}
plus_cats = {'N' : ['added_sources_N_key_z{}.dat'.format(str(s)) for s in [3,4,5,6,7,8]], 
             'S' : ['added_sources_S_key_z{}.dat'.format(str(s)) for s in [3,4,5,6,7,8]],
             'name' : ['plus_z{}'.format(str(s)) for s in [3,4,5,6,7,8]]}

all_cats = [quiescent_cats, emitters_cats, ivas_cat, plus_cats]

## Associate CLEAR Goods-N pointings with overlapping 3DHST pointings.
overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11']}

#
#-------------------------------------------------------------------------------  

def prep_clear(field, make_asn = True, check_background = False, 
    run_prep = True):

    import threedhst.prep_flt_astrodrizzle as init
    import unicorn.interlace_test as test

    if make_asn:

        unicorn.candels.make_asn_files()

    if check_background:

        #os.chdir('/Users/imomcheva/Work/CLEAR/RAW')
        os.chdir(paths['path_to_RAW'])
        asn_files = glob.glob('*asn.fits')
        for asn in asn_files:
            root = asn.split('_asn')[0]
            png = glob.glob(root+'*orbit.*')
            if (len(png) > 0):
                continue
            try:
                mywfc3.bg.show_orbit_limbangle(asn = [root])
            except:
                fp = open(root+'_orbit.failed','w')
                fp.write('---')
                fp.close()
        os.chdir(paths['path_to_PREPARE'])
        
    if remove_satelites:
        
        #os.chdir('/Users/imomcheva/Work/CLEAR/RAW/')
        os.chdir(paths['path_to_RAW'])
        unicorn.prepare.make_IMA_FLT(raw='icxt31r3q_raw.fits', pop_reads=[2])    
        os.chdir(paths['path_to_PREPARE'])
        
    if run_prep:

        os.chdir(paths['path_to_PREPARE'])

        grism_files = glob.glob(field + '*G102_asn.fits')
        direct_files = glob.glob(field + '*F105W_asn.fits')

        for direct, grism in zip(direct_files, grism_files):
            init.prep_direct_grism_pair(direct_asn=direct, grism_asn=grism, 
                radec=paths['path_to_ref_files']+'REF/goodss_3dhst.v4.1.radec.cat',
                raw_path = paths['path_to_RAW'], mask_grow=8, 
                scattered_light=False, final_scale = None,
                skip_direct=False, ACS=False, align_threshold=6)


#-------------------------------------------------------------------------------  

def interlace_clear(field):
    """
    ** Step 1. of Interlace steps. **

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

    Produces:
    - *inter.fits
    - *inter.reg
    - *inter_seg.fits
    - *ref_inter.fits
    - *radec.dat


    Checks:
   - Make sure *G102_ref_inter.fits and *inter_seg.fits are all aligned!
   - Load each *inter.reg into each *G102_ref_inter.fits to check 
     objects overlap catalog.
   - These checks must be done by Visit, not field, since Visits are rotated.
     So do for example,

   ds9 GN7-38-315-F105W_inter.fits GN7-38-315-G102_ref_inter.fits GN7-38-315-G102_inter_seg.fits &


    Parameters:
        field : string
            The GOODS field to process. 
    """

    print "Processing field {}".format(field)

    from unicorn.reduce import adriz_blot_from_reference as adriz_blot

    path_to_REF = paths['path_to_ref_files'] + 'REF/'
      
    NGROWX = 200
    NGROWY = 30
    if 'GDN' in field:
        pad = 500
    else:
        pad = 60

    if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
        CATALOG = path_to_REF + 'goodss_3dhst.v4.0.F125W_conv_fix.cat'
        REF_IMAGE = path_to_REF + 'goodss_3dhst.v4.0.F125W_orig_sci.fits'
        REF_EXT = 0
        SEG_IMAGE = path_to_REF + 'goodss_3dhst.v4.0.F160W_seg.fits'

    elif 'GN' in field or 'GDN' in field:
        CATALOG = path_to_REF + 'goodsn_3dhst.v4.0.F125W_conv.cat'
        REF_IMAGE = path_to_REF + 'goodsn_3dhst.v4.0.F125W_orig_sci.fits'
        REF_EXT = 0
        SEG_IMAGE = path_to_REF + 'goodsn_3dhst.v4.0.F160W_seg.fits'

    REF_FILTER = 'F125W'

    grism = glob.glob(field+'*G102_asn.fits')
    print "grism: {}".format(grism)

    for i in range(len(grism)):
        pointing=grism[i].split('_asn')[0]
        print pointing
        
        # Find whether pointing begins with a direct image (0) or grism (1).
        ref_exp = 0 #find_pointing_start(pointing)
        print "ref_exp: {}, pointing: {}".format(ref_exp, pointing)

        adriz_blot(pointing=pointing, pad=pad, NGROWX=NGROWX, 
            NGROWY=NGROWY, growx=2, growy=2, auto_offsets=True, 
            ref_exp=ref_exp, ref_image=REF_IMAGE, ref_ext=REF_EXT, 
            ref_filter=REF_FILTER, seg_image=SEG_IMAGE, 
            cat_file=CATALOG, grism='G102')     
       
        if 'GN5-42-346' in pointing:
             # These stare images had a bad dither. Set to 1x1 binning.
            print "Binning 1x1!"
            growx=1
            growy=1
        else:
            growx=2
            growy=2 

        # Interlace the direct images. Set ref_exp to always be zero.                                                                             
        unicorn.reduce.interlace_combine(pointing.replace('G102','F105W'), 
            view=False, use_error=True, make_undistorted=False, pad=pad, 
            NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
            growx=growx, growy=growy, auto_offsets=True, ref_exp=0)
        # Interlace the grism images.
        unicorn.reduce.interlace_combine(pointing, view=False, 
            use_error=True, make_undistorted=False, pad=pad, 
            NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
            growx=2, growy=2, auto_offsets=True, ref_exp=ref_exp)



    print "*** interlace_clear step complete ***"


 #-------------------------------------------------------------------------------  

def model_clear(field):
    """
    ** Step 2. of Interlace steps. **

    Creates model contam images. 

    Produces:
    - *inter_model.pkl
    - *inter_model.fits
    - *inter_0th.reg
    - *maskbg.dat
    - *maskbg.png

    Checks:
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


    Parameters:
        field : string
            The GOODS field to process. 
    """
    grism = glob.glob(field+'*G102_asn.fits')
    
    for i in range(len(grism)):
        root = grism[i].split('-G102')[0]
        #m0 = unicorn.reduce.GrismModel(root='goodss-01')
        #model_list = m0.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
        model = unicorn.reduce.process_GrismModel(root=root, 
            grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, 
            REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, 
            model_slope=0, direct='F105W', grism='G102', 
            BEAMS=['A', 'B', 'C', 'D','E'], align_reference=False)
        if not os.path.exists(os.path.basename(model.root) + '-G102_maskbg.dat'):
            model.refine_mask_background(threshold=0.002, grow_mask=14, 
                update=True, resid_threshold=4, clip_left=640, 
                save_figure=True, interlace=True)

    print "*** model_clear step complete ***"


#-------------------------------------------------------------------------------  

def extract_clear(field, tab):
    """

    ** Step 3. of Interlace steps. **


    Parameters:
        field : string
            The GOODS field to process. 
        tab : dictionary
            Values for each source.

    Checks:
        *2D.png should show some extracted spectra (most may be empty)

    """  

    grism = glob.glob(field+'*G102_asn.fits')
    #if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
    #    tab = Table.read(path_to_REF + cat['S'], format='ascii')
    #elif 'GN' in field or 'GDN' in field:
    #    tab = Table.read(path_to_REF + cat['N'], format='ascii')

    for i in range(len(grism)):
        root = grism[i].split('-G102')[0]
        model = unicorn.reduce.process_GrismModel(root=root, 
            grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, 
            REFINE_MAG_LIMIT=21, make_zeroth_model=False, 
            use_segm=False, model_slope=0, direct='F105W', grism='G102', 
            BEAMS=['A', 'B', 'C', 'D','E'], align_reference=False)
        for id in [id for id in tab['ID'] if id in model.cat.id]:
            try:
                model.twod_spectrum(id=id, grow=1, miny=-36, maxy=None, 
                                    CONTAMINATING_MAGLIMIT=23, 
                                    refine=False, verbose=False, 
                                    force_refine_nearby=False, 
                                    USE_REFERENCE_THUMB=True,
                                    USE_FLUX_RADIUS_SCALE=3, 
                                    BIG_THUMB=False, extract_1d=True)
                model.show_2d(savePNG=True, verbose=True)
                print 'Extracted {}'.format(id)
            except:
                continue

    print "*** extract_clear step complete ***"              


#-------------------------------------------------------------------------------  

def stack_clear(field, tab, cat, catname):
    """

    Parameters:
        field : string
            The GOODS field to process. Needs to know to reference GN or GS 
            catalog.
        tab : dictionary
            Values for each source.
        cat : dictionary
            Keys are 'N' or 'S' and values are the string names of 
            the catalog files.
        catname : string
            Name of the catalog, for naming output directory.

    Checks:
        In *stack.png check that the contam models reasonably match actual
        image and that the contam-flux is reasonable.
        python> !open *stack.png

    Notes:
        This should stack ALL the *2D* files present in the directory 
        that have the same id, REGARDLESS of the field name.
    """
    
    # Need a function that knows to search over all correlating 3dHST pointings
    # when given a specific CLEAR pointing.

    from unicorn.hudf import Stack2D

    grism = glob.glob(field+'*G102_asn.fits')
    #if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
    #    tab = Table.read(path_to_REF + cat['S'], format='ascii')
    #elif 'GN' in field or 'GDN' in field:
    #    tab = Table.read(path_to_REF + cat['N'], format='ascii')
    for i in range(len(grism)):
        root = grism[i].split('-G102')[0]
        model = unicorn.reduce.process_GrismModel(root=root, 
                                                  grow_factor=2, growx=2, 
                                                  growy=2, MAG_LIMIT=24,
                                                  REFINE_MAG_LIMIT=21, 
                                                  make_zeroth_model=False, 
                                                  use_segm=False,
                                                  model_slope=0, 
                                                  direct='F105W', 
                                                  grism='G102', 
                                                  BEAMS=['A', 'B', 'C', 'D','E'],
                                                  align_reference=False)
        for id in tab['ID']:
            if (id in model.cat.id):
                #and (not os.path.exists(field+'-G102_{}.2D.fits'.format(id))):
                try:
                    search='*-*-*-G102'
                    print 'searching %s*%05d.2D.fits'%(search, id)
                    spec = Stack2D(id=np.int(id), inverse=False, 
                                   scale=[1,99], fcontam=2.,
                                   ref_wave = 1.05e4,
                                   root='{}-G102'.format(field), 
                                   search='*-*-*-G102', files=None, 
                                   go=True, new_contam=False)
                except:
                    continue

    cleanup_extractions(field=field, cat=cat, catname=catname)

    print "*** stack_clear step complete ***"
    

#------------------------------------------------------------------------------- 

def cleanup_extractions(field, cat, catname):
    """ Moves all *1D*, *2D*, and *stack* files to a directory in Extractions
    named /<field>/<catalog>_yyyy-mm-dd/.
    Then tars the directory in Extractions.

    Parameters:
        field : string
            The GOODS field to process. Needs to know to reference GN or GS 
            catalog.
        cat : dictionary
            Keys are 'N' or 'S' and values are the string names of 
            the catalog files.
        catname : string
            Name of the catalog, for naming output directory.

    """
    path_to_Extractions = paths['path_to_Extractions']

    files = glob.glob('*1D*') + glob.glob('*2D*') + glob.glob('*stack*')
    print files

    dirname = os.path.join(path_to_Extractions, field, '{}_{}'.format(catname, 
        time.strftime('%Y-%m-%d')))

    #dirname = 'extractions_{}_{}'.format(catname, time.strftime('%d%B%Y'))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    print "Moving extractions from catalog {} to {}".format(catname, dirname)

    for f in files:
        shutil.move(f, os.path.join(dirname,f))

    # Now tar the directory.
    shutil.make_archive(os.path.join(path_to_Extractions, field, '{}_extractions_{}'.format(field, catname)), 
        'gztar', dirname)


#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------  

def clear_pipeline_main(fields, do_steps, cats):
    """
    """
    path_to_REF = paths['path_to_ref_files'] + 'REF/'

    for field in fields:
        print "***Beginning field {}***".format(field)
        print ""

        # Choose the field for the catalogs, where the catalog options include emitters and quiescent.
        if 'GS' in field or 'GDS' in field or 'ERSPRIME' in field:
            cats_field = cats['S']
        elif 'GN' in field or 'GDN' in field:
            cats_field = cats['N']

        for cat, catname in zip(cats_field, cats['name']):
            print "***Beginning catalog {}***".format(cat)
            print ""

            if 'GN' in field:
                # add primary CLEAR pointing to fields
                overlapping_fields_all = overlapping_fields[field]
                overlapping_fields_all.append(field)
                for overlapping_field in overlapping_fields_all:
                    print "***Beginning overlapping 13420 field {}***".format(overlapping_field)
                    print ""
                    if 1 in do_steps:
                        interlace_clear(field=overlapping_field)
                    if 2 in do_steps:
                        model_clear(field=overlapping_field)
                    if 3 in do_steps:
                        tab = Table.read(path_to_REF + cat, format='ascii')
                        extract_clear(field=overlapping_field, tab=tab)

            else:
                if 1 in do_steps:
                    interlace_clear(field=field)
                if 2 in do_steps:
                    model_clear(field=field)

                if 3 in do_steps:
                    tab = Table.read(path_to_REF + cat, format='ascii')
                    extract_clear(field=field, tab=tab)
                
            if 4 in do_steps:
                stack_clear(field=field, tab=tab, cat=cat, catname=catname) 
        


if __name__=='__main__':
    # all_cats = [quiescent_cats, emitters_cats, ivas_cat, plus_cats]
    fields = ['ERSPRIME'] 
    # Steps 3 and 3 should always be done together (the output directory will be messed up otherwise)
    do_steps = [3,4]
    clear_pipeline_main(fields=fields, do_steps=do_steps, cats=all_cats[0])

