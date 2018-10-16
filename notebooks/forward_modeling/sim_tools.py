__author__ = 'Vince.ec'

from grizli import model
from grizli import multifit
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from astropy import wcs
from astropy.io import fits
from glob import glob
import os

def Extract_BeamCutout(target_id, grism_file, mosaic, seg_map, instruement, catalog):
    flt = model.GrismFLT(grism_file = grism_file ,
                          ref_file = mosaic, seg_file = seg_map,
                            pad=200, ref_ext=0, shrink_segimage=True, force_grism = instrument)
    
    # catalog / semetation image
    ref_cat = Table.read(catalog ,format='ascii')
    seg_cat = flt.blot_catalog(ref_cat,sextractor=False)
    flt.compute_full_model(ids=seg_cat['id'])
    beam = flt.object_dispersers[target_id][2]['A']
    co = model.BeamCutout(flt, beam, conf=flt.conf)
    
    PA = np.round(fits.open(grism_file)[0].header['PA_V3'] , 1)
    
    co.write_fits(root='beams/o{0}'.format(PA), clobber=True)

    ### add EXPTIME to extension 0
    
    
    fits.setval('beams/o{0}_{1}.{2}.A.fits'.format(PA, target_id, instrument), 'EXPTIME', ext=0,
            value=fits.open('beams/o{0}_{1}.{2}.A.fits'.format(PA, target_id, instrument))[1].header['EXPTIME'])   

def Scale_model(data, sigma, model):
    return np.sum(((data * model) / sigma ** 2)) / np.sum((model ** 2 / sigma ** 2))

class Gen_spec(object):
    def __init__(self, beam_102, beam_141, redshift, spec_file_102, spec_file_141, 
                minwv_102 = 7800, maxwv_102 = 11350,
                minwv_141 = 10700, maxwv_141 = 16000):
        self.beam_102 = model.BeamCutout(fits_file = beam_102)
        self.beam_141 = model.BeamCutout(fits_file = beam_141)
        self.redshift = redshift       
        """ 


        """

        self.gal_wv_102, self.gal_fl_102, self.gal_er_102 = np.load(spec_file_102)

        IDX = [U for U in range(len(self.gal_wv_102)) if minwv_102 <= self.gal_wv_102[U] <= maxwv_102]

        self.gal_wv_rf_102 = self.gal_wv_102[IDX] / (1 + self.redshift)
        self.gal_wv_102 = self.gal_wv_102[IDX]
        self.gal_fl_102 = self.gal_fl_102[IDX]
        self.gal_er_102 = self.gal_er_102[IDX]

        self.gal_wv_141, self.gal_fl_141, self.gal_er_141 = np.load(spec_file_141)

        IDX = [U for U in range(len(self.gal_wv_141)) if minwv_141 <= self.gal_wv_141[U] <= maxwv_141]

        self.gal_wv_rf_141 = self.gal_wv_141[IDX] / (1 + self.redshift)
        self.gal_wv_141 = self.gal_wv_141[IDX]
        self.gal_fl_141 = self.gal_fl_141[IDX]
        self.gal_er_141 = self.gal_er_141[IDX]

        ## Get sensitivity function
        flat = self.beam_102.flat_flam.reshape(self.beam_102.beam.sh_beam)
        fwv, ffl, e = self.beam_102.beam.optimal_extract(flat, bin=0)
        
        self.filt_102 = interp1d(fwv, ffl)(self.gal_wv_102)
        
        flat = self.beam_141.flat_flam.reshape(self.beam_141.beam.sh_beam)
        fwv, ffl, e = self.beam_141.beam.optimal_extract(flat, bin=0)
        
        self.filt_141 = interp1d(fwv, ffl)(self.gal_wv_141)
        
    def Sim_spec(self, model_file, model_redshift = 0, dust = 0):
        model_wv, model_fl = np.load(model_file)
        
        if model_redshift ==0:
            model_redshift = self.redshift 
        
        ## Compute the models
        self.beam_102.compute_model(spectrum_1d=[model_wv*(1+model_redshift),model_fl])
        self.beam_141.compute_model(spectrum_1d=[model_wv*(1+model_redshift),model_fl])

        ## Extractions the model (error array here is meaningless)
        w, f, e = self.beam_102.beam.optimal_extract(self.beam_102.model , bin=0)

        ifl = interp1d(w, f)(self.gal_wv_102)
        
        C = Scale_model(self.gal_fl_102, self.gal_er_102, ifl / self.filt_102)

        self.fl_102 = C * ifl / self.filt_102
        
        ## Extractions the model (error array here is meaningless)
        w, f, e = self.beam_141.beam.optimal_extract(self.beam_141.model , bin=0)

        ifl = interp1d(w, f)(self.gal_wv_141)
        
        C = Scale_model(self.gal_fl_141, self.gal_er_141, ifl / self.filt_141)

        self.fl_141 = C * ifl / self.filt_141
        
class Gen_MB_spec(object):
    def __init__(self, beams, redshift, spec_file_102, spec_file_141, minwv_102 = 7800, maxwv_102 = 11350,
                minwv_141 = 10700, maxwv_141 = 16000):
        self.mb = multifit.MultiBeam(beams)
        self.redshift = redshift
        
        """ 


        """
        self.gal_wv_102, self.gal_fl_102, self.gal_er_102 = np.load(spec_file_102)

        IDX = [U for U in range(len(self.gal_wv_102)) if minwv_102 <= self.gal_wv_102[U] <= maxwv_102]

        self.gal_wv_rf_102 = self.gal_wv_102[IDX] / (1 + self.redshift)
        self.gal_wv_102 = self.gal_wv_102[IDX]
        self.gal_fl_102 = self.gal_fl_102[IDX]
        self.gal_er_102 = self.gal_er_102[IDX]

        self.gal_wv_141, self.gal_fl_141, self.gal_er_141 = np.load(spec_file_141)

        IDX = [U for U in range(len(self.gal_wv_141)) if minwv_141 <= self.gal_wv_141[U] <= maxwv_141]

        self.gal_wv_rf_141 = self.gal_wv_141[IDX] / (1 + self.redshift)
        self.gal_wv_141 = self.gal_wv_141[IDX]
        self.gal_fl_141 = self.gal_fl_141[IDX]
        self.gal_er_141 = self.gal_er_141[IDX]
        
        ## Get sensitivity function
        flat = self.mb.optimal_extract(self.mb.flat_flam[self.mb.fit_mask])
        
        self.filt_102 = interp1d(flat['G102']['wave'], flat['G102']['flux'])(self.gal_wv_102)
        self.filt_141 = interp1d(flat['G141']['wave'], flat['G141']['flux'])(self.gal_wv_141)
        
        
    def Sim_spec(self, model_file, model_redshift = 0):
        mwv, mfl = np.load(model_file)
        
        if model_redshift ==0:
            model_redshift = self.redshift 
        
        ## Create model
        spec = self.mb.get_flat_model(spectrum_1d = [mwv *(1 + model_redshift), mfl])

        ## Extract model and flat
        sp = self.mb.optimal_extract(spec)
        
        ifl = interp1d(sp['G102']['wave'], sp['G102']['flux'])(self.gal_wv_102)
        
        C = Scale_model(self.gal_fl_102, self.gal_er_102, ifl / self.filt_102)

        self.fl_102 = C * ifl / self.filt_102
        
        ifl = interp1d(sp['G141']['wave'], sp['G141']['flux'])(self.gal_wv_141)
        
        C = Scale_model(self.gal_fl_141, self.gal_er_141, ifl / self.filt_141)

        self.fl_141 = C * ifl / self.filt_141