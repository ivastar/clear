import os
import glob

import pyfits
import numpy as np

import threedhst
import unicorn

def prep_clear(make_asn = True, check_background = False, run_prep = True):

	import threedhst.prep_flt_astrodrizzle as init
	import unicorn.interlace_test as test

	if make_asn:

		unicorn.candels.make_asn_files()

	if check_background:

		os.chdir('/Users/imomcheva/Work/CLEAR/RAW')
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
	    os.chdir('/Users/imomcheva/Work/CLEAR/PREPARE/')
        
    if remove_satelites:
        
        os.chdir('/Users/imomcheva/Work/CLEAR/RAW/')
        unicorn.prepare.make_IMA_FLT(raw='icxt31r3q_raw.fits', pop_reads=[2])    
        os.chdir('/Users/imomcheva/Work/CLEAR/PREPARE/')
        
	if run_prep:

		os.chdir('/Users/imomcheva/Work/CLEAR/PREPARE/')

		grism_files = glob.glob('GS3-3[01]*G102_asn.fits')
		direct_files = glob.glob('GS3-3[01]*F105W_asn.fits')

		for direct, grism in zip(direct_files, grism_files):
			init.prep_direct_grism_pair(direct_asn=direct, grism_asn=grism, radec='../REF/goodss_3dhst.v4.1.radec.cat',
				raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = None,
				skip_direct=False, ACS=False, align_threshold=6)

def interlace_clear():
	
    from unicorn.reduce import adriz_blot_from_reference as adriz_blot
      
	NGROWX = 400
	NGROWY = 30
	pad = 60

	CATALOG = '/3DHST/Photometry/Release/v4.0/GOODS-S/Detection/goodss_3dhst.v4.0.F125W_conv.cat'
	REF_IMAGE = '../REF/goodss_3dhst.v4.0.F125W_orig_sci.fits'
	REF_EXT = 0
	SEG_IMAGE = '/3DHST/Photometry/Release/v4.0/GOODS-S/Detection/goodss_3dhst.v4.0.F160W_seg.fits'
	REF_FILTER = 'F125W'

	grism = glob.glob('GS3*G102_asn.fits')

	for i in range(len(grism)):
		pointing=grism[i].split('_asn')[0]
        
        if pointing.startswith('GS3-32-173') or pointing.startswith('GS3-33-173'):
            ref_exp = 1
        else:
            ref_exp = 0
    
		adriz_blot(pointing=pointing, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, growx=2, 
			growy=2, auto_offsets=True, ref_exp=ref_exp, ref_image=REF_IMAGE, ref_ext=REF_EXT, 
			ref_filter=REF_FILTER, seg_image=SEG_IMAGE, cat_file=CATALOG, grism='G102')                                                                                  
		unicorn.reduce.interlace_combine(pointing.replace('G102','F105W'), view=False, use_error=True, 
			make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
			growx=2, growy=2, auto_offsets=True, ref_exp=0)
		unicorn.reduce.interlace_combine(pointing, view=False, use_error=True, 
			make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
			growx=2, growy=2, auto_offsets=True, ref_exp=0)
    
def model_clear():
    
    	grism = glob.glob('GS3*G102_asn.fits')
    	
		for i in range(len(grism)):
			root = grism[i].split('-G102')[0]
			#m0 = unicorn.reduce.GrismModel(root='goodss-01')
			#model_list = m0.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
            model = unicorn.reduce.process_GrismModel(root=root, grow_factor=2, growx=2, growy=2, 
            	MAG_LIMIT=24, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, 
            	model_slope=0, direct='F105W', grism='G102', BEAMS=['A', 'B', 'C', 'D','E'], 
            	align_reference=False)
			if not os.path.exists(os.path.basename(model.root) + '-G102_maskbg.dat'):
				model.refine_mask_background(threshold=0.002, grow_mask=14, update=True, resid_threshold=4, 
					clip_left=640, save_figure=True, interlace=True)
    