import astropy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table, join
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


plt.rcParams['text.usetex'] = True
plt.close('all')
grizli_cat_dir = '/Users/rsimons/Dropbox/clear/grizli_extractions_v3.0/grizli_v3.0_cats'

fig, axes = plt.subplots(2, 3, figsize = (12., 7.))

lam_bins_G102 = [(0.8, 0.9), (0.9, 1.0), (1.0, 1.15)]
lam_bins_G141 = [(1.15, 1.35), (1.35, 1.55), (1.55, 1.75)]

lines = [('Ha', 0.6563), ('OIII', 0.5007), ('Hb', 0.4861), ('OII', 0.3728)]

lam_bins_both = [lam_bins_G102, lam_bins_G141]

for g in np.arange(2):
    lam_bins = lam_bins_both[g]
    if g == 0: 
        t_str = 'G102'
        cm = plt.cm.Blues_r
        cm = truncate_colormap(cm, 0., 0.8)
        (vmn, vmx) = (0., 15.)
    if g == 1: 
        t_str = 'G141'
        cm = plt.cm.Reds_r
        cm = truncate_colormap(cm, 0., 0.8)
        (vmn, vmx) = (0., 6.)


    #cm = plt.cm.viridis

    for f in ['S', 'N']:
        cat = fits.open(grizli_cat_dir + '/' + 'GD%s_lines_grizli_master.fits'%f)
        cat = cat[1].data

        for ll in np.arange(len(lam_bins)):
            ax = axes[g, ll]

            (lmin, lmax) = lam_bins[ll]


            ax.annotate('%.1f '%lmin + r'$<\lambda_{\text{obs}}<$' + ' %.1f'%lmax, (0.95, 0.05), \
                        xycoords = 'axes fraction', ha = 'right', va = 'bottom', fontsize = 16, fontweight = 'bold')


            for line in lines:
                lam_obs = line[-1] * (1+cat['z_MAP'])
                gd = np.where((lam_obs > lmin) & (lam_obs < lmax) & (cat['%s_FLUX'%line[0]]> 0.))[0]
                x = cat['%s_FLUX'%line[0]][gd] * 1.e-17
                y = cat['%s_FLUX'%line[0]][gd]/cat['%s_FLUX_ERR'%line[0]][gd]
                T_exp = cat['T_%s'%t_str][gd]/3600.
                scat = ax.scatter(x, y, c = T_exp,  s = 4., cmap = cm, \
                                 norm = plt.Normalize(vmin = vmn, vmax = vmx), zorder = 10., \
                                 edgecolor = 'black', linewidths = 0.1)

            if (f == 'S'):# & (zz == 3): 
                #cax = fig.add_axes([0.10, 0.93, 0.15, 0.02])
                bbox = (0.05, 0., 0.9, 0.98)
                cax = inset_axes(ax, width="40%", height="5%", loc=2, bbox_to_anchor = bbox, bbox_transform = ax.transAxes)
                cbar = plt.colorbar(scat, cax = cax, orientation = 'horizontal')
                cbar.ax.tick_params(labelsize = 8) 
                
                cbar.set_label('%s Exposure Time (hr)'%t_str, fontsize = 10, labelpad = 2.)
        
for ax in axes.ravel(): 
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(0.1, 200)
    ax.set_xlim(2.e-18, 2.e-15)
    ax.set_yticks([1, 10, 100])
    ax.set_yticklabels(['1', '10', '100'])
    ax.plot([2.e-18, 2.e-15], [0.2, 200], 'k-', zorder = 1, alpha = 0.4)

for ax in axes[1]:
    ax.set_xlabel('Line Flux (erg s$^{-1}$ cm$^{-2}$)')

for ax in axes[:,0]:
    ax.set_ylabel('S/N')


fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/survey_paper/flux_snr_exptime.png', dpi = 300)


