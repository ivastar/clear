import astropy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table, join
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
plt.close('all')
grizli_cat_dir = '/Users/rsimons/Dropbox/clear/catalogs/grizli_v3.0_cats'
tdhst_cat_dir  = '/Users/rsimons/Dropbox/clear/catalogs/'
fields = ['GDN', 'GDS']


fig, axes = plt.subplots(1, 2, figsize = (8, 4))


fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/survey_paper/z_versus_z.png', dpi = 300)


