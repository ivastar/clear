import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
from clear_local.utils import tools
import numpy as np
from numpy import *
from astropy.table import join, Table
from astropy.convolution import convolve_fft, Gaussian2DKernel


fig, axes = plt.subplots(2,2, figsize = (6,5.5))
fig2, axes2 = plt.subplots(2,2, figsize = (6,5.5))

mass_min = 8.3
mass_max = 11.
sfr_min  = -2.
sfr_max  = 2.5
z_min = 0.7
z_max = 2.5
nbins = 25.
nbins2 = 40.

sn_limit = 5.

t_min_102 = -1.
t_min_141 = -1.
t_max_102 = 50.
t_max_141 = 50.
lvl_p = 0.68

dic = {}








mass_use = np.array([])
sfr_use  = np.array([])



lines = [('Ha_FLUX','Ha_FLUX_ERR'     , 'red', 6563),\
         ('OIII_FLUX','OIII_FLUX_ERR'   , 'green', 5007),\
         ('Hb_FLUX','Hb_FLUX_ERR'     , 'darkgoldenrod', 4862),\
         ('OII_FLUX','OII_FLUX_ERR' , 'blue', 3727)]#[0:1]


dic_detect = {}

for line in lines:
    dic_detect[line[0]] = {}
    dic_detect[line[0]]['mass_use_detected'] = np.array([])
    dic_detect[line[0]]['sfr_use_detected']  = np.array([])
    dic_detect[line[0]]['mass_use_chance'] = np.array([])
    dic_detect[line[0]]['sfr_use_chance']  = np.array([])


for field in ['s', 'n']:
    eazy_fits = tools.load_eazy_catalog(field_name = field, cat_version = 'v4.5')
    grizli_fits = tools.load_grizli_catalog(field_name = field, cat_version = 'v3.0')

    eazy_cat   = Table(eazy_fits[1].data)
    grizli_cat = Table(grizli_fits[1].data)

    grizli_cat.rename_column('ID', 'id')
    grizli_cat['id'] = grizli_cat['id'].astype('float')
    combined_cat = join(eazy_cat, grizli_cat, keys = 'id')

    mass = combined_cat['mass']
    sfr  = combined_cat['SFR']  
    z    = combined_cat['z_MAP']

    t_g102 = combined_cat['T_G102']
    t_g141 = combined_cat['T_G141']


    crit1 = (log10(mass) > mass_min)   & (log10(mass) < mass_max)
    crit2 = (log10(sfr) > sfr_min)   & (log10(sfr) < sfr_max)
    crit3a = (t_g102/3600. < t_max_102) & (t_g102/3600. > t_min_102)
    crit3b = (t_g141/3600. < t_max_141) & (t_g141/3600. > t_min_141)

    crit3 = crit3a & crit3b

    crit4 = (z > z_min) & (z < z_max)
    use =   np.where(crit1 & crit2 & crit3 & crit4)[0]

    mass_use = np.concatenate((mass_use, log10(mass[use])))
    sfr_use  = np.concatenate((sfr_use, log10(sfr[use])))



    for line in lines:
        detected = np.zeros(len(use))
        sn     = combined_cat[line[0]][use]/combined_cat[line[1]][use]
        use2 = use[sn > sn_limit]
        dic_detect[line[0]]['mass_use_detected'] = np.concatenate((dic_detect[line[0]]['mass_use_detected'], log10(mass[use2])))
        dic_detect[line[0]]['sfr_use_detected']  = np.concatenate((dic_detect[line[0]]['sfr_use_detected'] , log10(sfr [use2])))

        

        lam_102 = [8000., 11500.]
        lam_141 = [10750., 17000.]

        z_102_min = lam_102[0]/line[3] - 1.
        z_102_max = lam_102[1]/line[3] - 1.

        z_141_min = lam_141[0]/line[3] - 1.
        z_141_max = lam_141[1]/line[3] - 1.

        crit_g102 = (t_g102[use]/3600. > 0.) & (z[use] < z_102_max) & (z[use] > z_102_min)
        crit_g141 = (t_g141[use]/3600. > 0.) & (z[use] < z_141_max) & (z[use] > z_141_min)

        use_chance = use[crit_g102 | crit_g141]
        dic_detect[line[0]]['mass_use_chance'] = np.concatenate((dic_detect[line[0]]['mass_use_chance'], log10(mass[use_chance])))
        dic_detect[line[0]]['sfr_use_chance']  = np.concatenate((dic_detect[line[0]]['sfr_use_chance'] , log10(sfr [use_chance])))








dm = (mass_max - mass_min)/nbins
ds = (sfr_max - sfr_min)/nbins



bn_m = np.linspace(mass_min - dm, mass_max + dm, nbins)
bn_s = np.linspace(sfr_min  - ds, sfr_max  + ds,  nbins)
hst_out = np.histogram2d(sfr_use, mass_use, bins = [bn_s, bn_m])



xx = (bn_m[:-1] + bn_m[1:])/2.
yy = (bn_s[:-1] + bn_s[1:])/2.

X, Y = np.meshgrid(xx, yy)

Z          = hst_out[0]
for ax in axes.ravel():
    ax.contourf(X, Y, Z,  levels = np.linspace(-20., np.nanmax(Z), 20.),\
                cmap = 'Greys', zorder = 1.)

for l, line in enumerate(lines):


    mass_use_detected   = dic_detect[line[0]]['mass_use_detected'] 
    sfr_use_detected    = dic_detect[line[0]]['sfr_use_detected']  
    mass_use_chance   = dic_detect[line[0]]['mass_use_chance'] 
    sfr_use_chance    = dic_detect[line[0]]['sfr_use_chance']  

    dm2 = (mass_max - mass_min)/nbins2
    ds2 = (sfr_max - sfr_min)/nbins2


    bn_m2 = np.linspace(mass_min - dm2, mass_max + dm2, nbins2)
    bn_s2 = np.linspace(sfr_min  - ds2, sfr_max  + ds2,  nbins2)
    xx2 = (bn_m2[:-1] + bn_m2[1:])/2.
    yy2 = (bn_s2[:-1] + bn_s2[1:])/2.

    X2, Y2 = np.meshgrid(xx2, yy2)

    hst_out2 = np.histogram2d(sfr_use_chance, mass_use_chance, bins = [bn_s2, bn_m2])
    Z2         = hst_out2[0]


    hst_out_detected    = np.histogram2d(sfr_use_detected, mass_use_detected, bins = [bn_s, bn_m])
    hst_out_detected2    = np.histogram2d(sfr_use_detected, mass_use_detected, bins = [bn_s2, bn_m2])

    Z_detected = hst_out_detected[0]
    Z_detected2 = hst_out_detected2[0]

    Zd_ravel   = Z_detected.ravel()
    Zdect_sort = sort(Zd_ravel)
    Z_cumsum   = cumsum(Zdect_sort)

    lvls = []
    for lvl in [lvl_p]:
        lvls.append(Zdect_sort[argmin(abs(Z_cumsum - Z_cumsum[-1] * (1. - lvl)))])
    

    ax = axes.ravel()[l]
    ax2 = axes2.ravel()[l]
    ax.contour(X, Y, Z_detected,  levels = sort(lvls),\
                colors = line[2], zorder = 3.)

    detection_fraction = Z_detected2/Z2
    detection_fraction[Z2 < 3.] = np.nan
    cmap = plt.cm.viridis
    cmap.set_bad('k')
    im = ax2.imshow(detection_fraction, vmin = 0.0, vmax = 1.0, cmap = cmap, origin = 'lower')

    if line[0] == 'Ha_FLUX'  :  ln_ann = r'H$\alpha$'
    if line[0] == 'Hb_FLUX'  :  ln_ann = r'H$\beta$'
    if line[0] == 'OII_FLUX' :  ln_ann = '[OII]'
    if line[0] == 'OIII_FLUX':  ln_ann = '[OIII]'


    ax.annotate('S/N (%s) $>$ %i'%(ln_ann, sn_limit), (0.95, 0.05), 
                xycoords = 'axes fraction', ha = 'right', va = 'bottom', \
                color = line[2], fontsize = 12)

    ax2.annotate('%s'%(ln_ann), (0.95, 0.05), 
                xycoords = 'axes fraction', ha = 'right', va = 'bottom', \
                color = 'white', fontsize = 12)

axes[0,0].annotate('CLEAR survey\n%.1f $<$ z$_{grism}$ $<$ %.1f'%(z_min, z_max), \
                  (0.05, 0.95),  ha = 'left', va = 'top', xycoords = 'axes fraction',\
                  fontsize = 12)

axes2[0,0].annotate('CLEAR survey\n%.1f $<$ z$_{grism}$ $<$ %.1f'%(z_min, z_max), \
                  (0.05, 0.95),  ha = 'left', va = 'top', xycoords = 'axes fraction',\
                  color = 'white', fontsize = 12)

for ax in axes.ravel():
    ax.set_xlim(mass_min, mass_max)
    ax.set_ylim(sfr_min,  sfr_max)
    ax.set_xlabel(r'$\log$ Stellar Mass (M$_{\odot}$)', fontsize = 10)
    ax.set_ylabel(r'$\log$ Star-Formation Rate (M$_{\odot}$ yr$^{-1}$)', fontsize = 10)


if True:
    #plot tinkering for the detection fraction figure
    cax = fig2.add_axes([0.58, 0.94, 0.15, 0.02])
    cbar = fig2.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label('detection fraction', color = 'white', fontsize = 10, labelpad = 0)
    cbar.ax.xaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), fontsize = 8, color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'xticklabels'), fontsize = 8, color='white')

    xticks = np.arange(8.5, 11.5, 0.5)
    yticks = np.arange(-2, 3, 1)

    xticks_ind = [argmin(abs(tck - bn_m2)) for tck in xticks]
    yticks_ind = [argmin(abs(tck - bn_s2)) for tck in yticks]

    for ax in axes2.ravel():
        ax.set_xticks(xticks_ind)
        ax.set_yticks(yticks_ind)
        ax.set_xticklabels([''])
        ax.set_yticklabels([''])
        xlm_min = argmin(abs(mass_min - bn_m2))
        xlm_max = argmin(abs(mass_max - bn_m2))
        ylm_min = argmin(abs(sfr_min - bn_s2))
        ylm_max = argmin(abs(sfr_max - bn_s2))

        ax.set_xlim(xlm_min, xlm_max)
        ax.set_ylim(ylm_min, ylm_max)


    for ax in axes2[:,0]:
        ax.set_yticklabels(['$%i$'%tck for tck in yticks])    
        ax.set_ylabel(r'$\log$ Star-Formation Rate (M$_{\odot}$ yr$^{-1}$)', fontsize = 10)
    for ax in axes2[1,:]:
        ax.set_xticklabels(['$%.1f$'%tck for tck in xticks])  
        ax.set_xlabel(r'$\log$ Stellar Mass (M$_{\odot}$)', fontsize = 10)


for ax in axes[0]:
    ax.set_xlabel('')
    ax.set_xticklabels([''])

for ax in axes[:,1]:
    ax.set_ylabel('')
    ax.set_yticklabels([''])

fig.subplots_adjust(wspace = 0.05, hspace = 0.05, right = 0.98, top = 0.98, left = 0.10)
fig2.subplots_adjust(wspace = -0.20, hspace = 0.05, right = 1.0, top = 0.98, left = 0.08)
#fig.tight_layout()
#fig2.tight_layout()
fig_name = '/Users/rsimons/Dropbox/clear/figures/survey_paper/mstar_sfr.png'
fig2_name = '/Users/rsimons/Dropbox/clear/figures/survey_paper/mstar_sfr_detection_efficiency.png'
fig.savefig(fig_name, dpi = 600)

fig2.savefig(fig2_name, dpi = 600)


from PIL import Image
images = [Image.open(x) for x in [fig_name, fig2_name]]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

newfig_name = '/Users/rsimons/Dropbox/clear/figures/survey_paper/mstar_sfr_combined.png'
new_im.save(newfig_name)







