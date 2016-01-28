import os
import glob

import pyfits
import numpy as np

import threedhst
import unicorn

def figs_gn1(make_asn = True, check_background = True, run_prep = True, run_interlace=True
	make_model=True, ):

	import mywfc3.bg
	import threedhst.prep_flt_astrodrizzle as init
	from unicorn.reduce import adriz_blot_from_reference as adriz_blot
	import unicorn.interlace_test as test

	if make_asn:

		unicorn.candels.make_asn_files()

	if check_background:

		os.chdir('/Users/imomcheva/Work/FIGS/RAW')
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
	    os.chdir('/Users/imomcheva/Work/FIGS/PREPARE/')

	if run_prep:

		os.chdir('/Users/imomcheva/Work/FIGS/PREPARE/')

		grism_files = glob.glob('FIGS-GN1*G102_asn.fits')
		direct_files = glob.glob('FIGS-GN1*F105W_asn.fits')

		for direct, grism in zip(direct_files, grism_files):
			init.prep_direct_grism_pair(direct_asn=direct, grism_asn=grism, radec='../REF/goodsn_3dhst.v4.1.radec.cat',
				raw_path = '../RAW/', mask_grow=8, scattered_light=False, final_scale = None,
				skip_direct=False, ACS=False, align_threshold=6)

	if run_interlace:

		NGROWX = 200
		NGROWY = 30
		pad = 60

		CATALOG = '../REF/GN_aXe_v0.1.cat'
		REF_IMAGE = '../REF/goodsn_all_wfc3_ir_f125w_060mas_v1.0_drz.fits'
		REF_EXT = 0
		SEG_IMAGE = '../REF/finalseg_GN_v0.1.fits'
		REF_FILTER = 'F125W'

		grism = glob.glob('FIGS-GN1*G102_asn.fits')

		for i in range(len(grism)):
			pointing=grism[i].split('_asn')[0]
    
			adriz_blot(pointing=pointing, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, growx=1, 
				growy=1, auto_offsets=True, ref_exp=0, ref_image=REF_IMAGE, ref_ext=REF_EXT, 
				ref_filter=REF_FILTER, seg_image=SEG_IMAGE, cat_file=CATALOG, grism='G102')                                                                                  
    		unicorn.reduce.interlace_combine(pointing.replace('G102','F105W'), view=False, use_error=True, 
    			make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
    			growx=1, growy=1, auto_offsets=True, ref_exp=0)
    		unicorn.reduce.interlace_combine(pointing, view=False, use_error=True, 
    			make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, 
    			growx=1, growy=1, auto_offsets=True, ref_exp=0)

    if make_model:

    	grism = glob.glob('FIGS-GN1*G102_asn.fits')
    	
		for i in range(len(grism)):
			root = grism[i].split('-G102')[0]
			#m0 = unicorn.reduce.GrismModel(root=root)
			#model_list = m0.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
			model = unicorn.reduce.process_GrismModel(root=root, grow_factor=1, growx=1, growy=1, 
				MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, 
				model_slope=0, direct='F105W', grism='G102', BEAMS=['A', 'B', 'C', 'D', 'E'], 
				align_reference=True)
			if not os.path.exists(os.path.basename(model.root) + '-G141_maskbg.dat'):
				model.refine_mask_background(threshold=0.002, grow_mask=14, update=True, resid_threshold=4, 
					clip_left=640, save_figure=True, interlace=True)















