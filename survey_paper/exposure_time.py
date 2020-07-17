import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt



grizli_cat_dir = '/Users/rsimons/Dropbox/clear/grizli_extractions_v3.0/grizli_v3.0_cats'
fields = [('GDN', '-'), ('GDS', '--')]


fig, axes = plt.subplots(1, 2, figsize = (8, 4))

alp = 0.4

bns = np.arange(0, 30, 1)
for f, (fld, clr) in enumerate(fields):
    cat = fits.open(grizli_cat_dir + '/%s_lines_grizli_master.fits'%fld)
    t_g102 = cat[1].data['T_G102'] 
    t_g141 = cat[1].data['T_G141']

    axes[f].hist(t_g102/3600., color = 'blue', alpha =alp, bins = bns)
    axes[f].hist(t_g141/3600., color = 'red', alpha = alp, bins = bns)
    axes[f].set_xlabel('Exposure Time (hr)', fontsize = 14)
    axes[f].set_yscale('log')
    if fld == 'GDN': axes[f].annotate('GOODS-N', (0.95, 0.95), xycoords = 'axes fraction', ha = 'right', va = 'top', fontsize = 16)
    if fld == 'GDS': axes[f].annotate('GOODS-S', (0.95, 0.95), xycoords = 'axes fraction', ha = 'right', va = 'top', fontsize = 16)
axes[0].set_ylabel('Number of Objects', fontsize = 14)

alp = 0.7
axes[0].annotate('G102', (0.95, 0.17), xycoords = 'axes fraction', ha = 'right', va = 'top', color = 'blue', alpha = alp, fontsize = 16)
axes[0].annotate('G141', (0.95, 0.10), xycoords = 'axes fraction', ha = 'right', va = 'top', color = 'red', alpha = alp, fontsize = 16)




fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/survey_paper/exposure_time.png', dpi = 400)
