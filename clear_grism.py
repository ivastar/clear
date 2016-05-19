import os
import glob
import pyfits
import numpy as np

import threedhst
import unicorn

from astropy.table import Table
from astropy.io import ascii

from find_pointing_start import find_pointing_start
from set_paths import paths

# Define catalogs for S and N.
quiescent_cats = {'S':'UVJ_quiescent_goodss.dat', 'N':'UVJ_quiescent_goodsn.dat'}
steves_cats = {'S':'Steves_source_goodss_w_ids.dat', 'N':'Steves_source_goodsn_w_ids.dat'}


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


    Parameters:
        field : string
            The GOODS field to process. 
    """

    print "Processing field {}".format(field)

    from unicorn.reduce import adriz_blot_from_reference as adriz_blot

    path_to_REF = paths['path_to_ref_files'] + 'REF/'
      
    NGROWX = 200
    NGROWY = 30
    pad = 60

    if 'GS' in field or 'GDS' in field:
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
        
        # Find whether pointing begins with a direct image (0) or grism (1).
        ref_exp = find_pointing_start(pointing)
        print "ref_exp: {}, pointing: {}".format(ref_exp, pointing)

        """
        adriz_blot(pointing=pointing, pad=pad, NGROWX=NGROWX, 
            NGROWY=NGROWY, growx=2, growy=2, auto_offsets=True, 
            ref_exp=ref_exp, ref_image=REF_IMAGE, ref_ext=REF_EXT, 
            ref_filter=REF_FILTER, seg_image=SEG_IMAGE, 
            cat_file=CATALOG, grism='G102')                                                                                  
        unicorn.reduce.interlace_combine(pointing.replace('G102','F105W'), 
            view=False, use_error=True, make_undistorted=False, pad=pad, 
            NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
            growx=2, growy=2, auto_offsets=True, ref_exp=ref_exp)
        unicorn.reduce.interlace_combine(pointing, view=False, 
            use_error=True, make_undistorted=False, pad=pad, 
            NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
            growx=2, growy=2, auto_offsets=True, ref_exp=ref_exp)
        """

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

def extract_clear(field, cat):
    """

    ** Step 3. of Interlace steps. **


    Parameters:
        field : string
            The GOODS field to process. 
        cat : dict 
            Keys either 'S' or 'N' with values of source catalog filenames.
            The catalogs themselves should be ID, RA, DEC.  
    """  
    #ids = [37945,38822,39286,39411]
    path_to_REF = paths['path_to_ref_files'] + 'REF/'

    if cat == 'steves':
        cat = steves_cats
    elif cat == 'quiescent':
        cat = quiescent_cats
    print "Chosen catalogs: {}".format(cat)

    grism = glob.glob(field+'*G102_asn.fits')
    if 'GS' in field:
        tab = Table.read(path_to_REF + cat['S'], format='ascii')
    elif 'GN' in field:
        tab = Table.read(path_to_REF + cat['N'], format='ascii')
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

def stack_clear(field, cat):
    """

    Parameters:
        field : string
            The GOODS field to process. 
        cat : dict 
            Keys either 'S' or 'N' with values of source catalog filenames.
            The catalogs themselves should be ID, RA, DEC.

    Notes:
        How will this step work when trying to stack the CLEAR and the 3DHST 
        visits?
    """
    
    import os
    from unicorn.hudf import Stack2D
 
    path_to_REF = paths['path_to_ref_files'] + 'REF/'
    
    if cat == 'steves':
        cat = steves_cats
    elif cat == 'quiescent':
        cat = quiescent_cats
    print "Chosen catalogs: {}".format(cat)

    grism = glob.glob(field+'*G102_asn.fits')
    if 'GS' in field or 'GDS' in field:
        tab = Table.read(path_to_REF + cat['S'], format='ascii')
    elif 'GN' in field or 'GDN' in field:
        tab = Table.read(path_to_REF + cat['N'], format='ascii')
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
            if (id in model.cat.id) and (not os.path.exists(field+'-G102_{}.2D.fits'.format(id))):
                try:
                    search=field+'-*-*-G102'
                    print '%s*%05d.2D.fits'%(search, id)
                    spec = Stack2D(id=np.int(id), inverse=False, 
                                   scale=[1,99], fcontam=2.,
                                   ref_wave = 1.05e4,
                                   root=field+'-G102', 
                                   search=field+'-*-*-G102', files=None, 
                                   go=True, new_contam=False)
                except:
                    continue

    print "*** stack_clear step complete ***"
            

#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------  

if __name__=='__main__':
    cats =  ['quiescent', 'steves']
    fields = ['GN7']  #['GS1', 'GS2', 'GS3', 'GS5', 'GN4', 'GN5', 'GN7']
    for field in fields:
        print "***Beginning field {}***".format(field)
        print ""
        interlace_clear(field=field)
        #model_clear(field=field)
        #extract_clear(field=field, cat=cats[0])
        #stack_clear(field=field, cat=cats[0])
