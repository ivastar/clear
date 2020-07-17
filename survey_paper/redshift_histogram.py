import astropy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table, join
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
plt.close('all')
grizli_cat_dir = '/Users/rsimons/Dropbox/clear/catalogs/grizli_v2.1_cats'
tdhst_cat_dir  = '/Users/rsimons/Dropbox/clear/catalogs/'
fields = ['GDN', 'GDS']


fig, axes = plt.subplots(1, 1, figsize = (7, 4))
for f, field in enumerate(fields):
    gcat = fits.open('%s/%s_lines_grizli_master.fits'%(grizli_cat_dir, field))
    if 'N' in field: field_ann = 'GOODS-N'
    if 'S' in field: field_ann = 'GOODS-S'

    axes.hist(np.log10(1 + gcat[1].data['z_MAP']), bins = np.linspace(np.log10(1), np.log10(13), 50), alpha = 0.3, label = field_ann)


z_marks = np.arange(0, 13, 1)

z_ticks = [np.log10(1+z) for z in z_marks]
z_ann   = ['%i'%z for z in z_marks]
axes.set_xticks(z_ticks)
axes.set_xticklabels(z_ann)
axes.legend(loc = 1)
axes.set_xlabel(r'$z_{\textrm{spec, grism}}$', fontsize = 13)
axes.set_ylabel('number of galaxies', fontsize = 13)
fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/survey_paper/zgrizli_histogram.png', dpi = 300)

#stop at z = 4
#one bin for z > 5?
#hatch on one of them






