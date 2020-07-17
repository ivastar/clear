import astropy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table, join
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
plt.rcParams['xtick.major.size'] = 5.
plt.rcParams['ytick.major.size'] = 5.

plt.rcParams['xtick.minor.size'] = 3.
plt.rcParams['ytick.minor.size'] = 3.
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

plt.close('all')
grizli_cat_dir = '/Users/rsimons/Dropbox/clear/grizli_extractions_v3.0/grizli_v3.0_cats'
tdhst_cat_dir  = '/Users/rsimons/Dropbox/clear/catalogs/'
fields = ['GDN', 'GDS']


fig, axes = plt.subplots(2, 2, figsize = (8, 8))

do_points = False
if do_points:
    postpend = 'with_points'

else:
    postpend  = 'with_contour'


for f, field in enumerate(fields):
    gcat = fits.open('%s/%s_lines_grizli_master.fits'%(grizli_cat_dir, field))
    zcat = fits.open(tdhst_cat_dir + '/goods%s_3dhst.v4.4.cats/Eazy/goods%s_3dhst.v4.4.zout.fits'%(field[-1], field[-1]))
    zspecs     = zcat[1].data['z_spec']
    zphots     = zcat[1].data['z_phot']
    tids       = zcat[1].data['id']
    grizli_ids = gcat[1].data['id']

    dzoz_spec = []
    dzoz_phot = []
    z_grizlis_phot = []
    z_grizlis_spec = []
    z_spec_use = []
    z_phot_use = []
    def lg(z):
        return np.log10(z + 1)


    for ind, di in enumerate(grizli_ids):
        z_grizli   = gcat[1].data['z_MAP'][ind]
        z_spec = float(zspecs[tids == di])
        z_phot = float(zphots[tids == di])
        if (z_grizli < 5) & (z_spec < 5)  & (z_spec > 0.) & (z_grizli  > 0.):
            if do_points: axes[0, f].plot(lg(z_grizli), lg(z_spec), 'k.', zorder = 10, markersize = 4)
            z_grizlis_spec.append(z_grizli)
            z_spec_use.append(z_spec)
            dzoz_spec.append((z_grizli - z_spec)/(1+z_grizli))
        if (z_grizli < 5) & (z_phot < 5)  & (z_phot > 0):
            if do_points: axes[1, f].plot(lg(z_grizli), lg(z_phot), 'k.', zorder = 10, markersize = 4)
            dzoz_phot.append((z_grizli - z_phot)/(1+z_grizli))
            z_phot_use.append(z_phot)
            z_grizlis_phot.append(z_grizli)

    label_fs = 14
    axes[0, f].set_xlabel(r'$z_{\textrm{spec, grism}}$', fontsize = label_fs)
    axes[1, f].set_xlabel(r'$z_{\textrm{spec, grism}}$', fontsize = label_fs)
    axes[0, f].set_ylabel(r'$z_{\textrm{spec, ground}}$', fontsize = label_fs)
    axes[1, f].set_ylabel(r'$z_{\textrm{phot}}$', fontsize = label_fs)


    bn_z = np.linspace(lg(0), lg(5), 125.)
    if not do_points: 
        zt = [np.array(z_spec_use), np.array(z_phot_use)]
        zg = [np.array(z_grizlis_spec), np.array(z_grizlis_phot)]
        for d in np.arange(2):
            hst_out = np.histogram2d(lg(zt[d]), lg(zg[d]), bins = [bn_z, bn_z])
            xx = (bn_z[:-1] + bn_z[1:])/2.
            X, Y = np.meshgrid(xx, xx)
            if False:
                kern = Gaussian2DKernel(0.2)
                Z    =  convolve_fft(hst_out[0], kern)
            Z = hst_out[0]
            axes[d,f].contourf(X, Y, Z, levels = [0.5, 5., 10., 15., 20.], cmap = 'viridis', zorder = 0)

    for ax in axes.ravel():
        minor_ticks      = np.arange(0, 6, 0.2)
        major_ticks      = np.arange(0, 6, 1.)
        major_ticklabels = np.array(['%i'%t for t in major_ticks])


        ax.set_xticks(lg(minor_ticks), minor = True)
        ax.set_yticks(lg(minor_ticks), minor = True)

        ax.set_xticks(lg(major_ticks), minor = False)
        ax.set_yticks(lg(major_ticks), minor = False)

        ax.set_xticklabels(major_ticklabels)
        ax.set_yticklabels(major_ticklabels)


        ax.set_xlim(lg(0),lg(5))
        ax.set_ylim(lg(0),lg(5))
        ax.plot([lg(0),lg(5)], [lg(0),lg(5)], 'r-', alpha = 0.1, zorder = 10)

    if 'N' in field: field_ann = 'GOODS-N'
    if 'S' in field: field_ann = 'GOODS-S'
    #axes[0, f].annotate(field_ann, (0.04, 0.87), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontsize = 25)
    #axes[1, f].annotate(field_ann, (0.04, 0.87), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontsize = 25)
    axes[0, f].set_title(field_ann, fontsize = 25)
    #axes[1, f].set_title(field_ann)

    def nmad(z1, z):
        delz = np.array(z1) - np.array(z)
        return 1.48 * np.median(abs(np.array(delz) - np.median(delz))/(1+np.array(z)))
    nmad_spec = nmad(z_spec_use, z_grizlis_spec)
    nmad_phot = nmad(z_phot_use, z_grizlis_phot)
    print (nmad_spec, nmad_phot)

    dzoz_spec = np.array(dzoz_spec)
    dzoz_phot = np.array(dzoz_phot)
    dzozs = [dzoz_spec, dzoz_phot]
    nmads = [nmad_spec, nmad_phot]

    nout_spec = len(np.where(abs(dzoz_spec - np.median(dzoz_spec)) > 5 * nmad_spec)[0])
    nout_phot = len(np.where(abs(dzoz_phot - np.median(dzoz_phot)) > 5 * nmad_phot)[0])
    print (nout_spec, len(dzoz_spec))
    foutlier_spec = nout_spec/len(dzoz_spec)
    foutlier_phot = nout_phot/len(dzoz_phot)
    ha = 0.03
    va = 0.98

    bbox_params = dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0') 
    axes[0, f].annotate('N = %i\n'%len(dzoz_spec)+r'$\sigma_{\text{NMAD}}$ = %.4f'%(nmad_spec) + '\n' + r'$f_{\text{outlier, $>$5$\sigma$}}$ = %.2f'%(foutlier_spec), \
                        (ha, va), ha = 'left', va = 'top', \
                        xycoords = 'axes fraction', fontsize = 15, \
                        color = 'Grey',  bbox=bbox_params)

    axes[1, f].annotate('N =  %i\n'%len(dzoz_phot)+r'$\sigma_{\text{NMAD}}$ = %.4f'%(nmad_phot) + '\n' + r'$f_{\text{outlier, $>$5$\sigma$}}$ = %.2f'%(foutlier_phot), \
                         (ha, va), ha = 'left', va = 'top', xycoords = 'axes fraction', fontsize = 15, \
                         color = 'Grey', bbox=bbox_params)


    xlms = [0.02, 0.06]    

    #([-0.01, 0.01], np.arange(-0.01, 0.02, 0.005), [ '-0.02', '', '0.0', '', '0.02'], np.linspace(-0.03, 0.03, 100)),
                 
    tck_props = [([-0.01, 0.01], np.arange(-0.01, 0.02, 0.005), [ '-0.02', '', '0.0', '', '0.02'], np.linspace(-0.03, 0.03, 100)),
                 ([-0.09, 0.09], np.arange(-0.08, 0.09, 0.04), [ '-0.08', '', '0.0', '', '0.08'], np.linspace(-0.10, 0.10, 50)),
                ]



    for d in np.arange(2):
        axinset = inset_axes(axes[d, f], width="100%", height="100%", 
                             bbox_to_anchor=(.65, .20, .30, .30),
                             bbox_transform=axes[d, f].transAxes)
        axinset.set_yticks([])
        axinset.set_xlabel(r'$\Delta$z/(1+z)', fontsize = 10)
        bn = tck_props[d][3]
        axinset.hist(dzozs[d], color = 'black', bins = bn)
        axinset.set_xticks(tck_props[d][1])
        axinset.set_xticklabels(tck_props[d][2])
        xlms = tck_props[d][0]
        axinset.set_xlim(xlms[0], xlms[1])
        axinset.set_ylim(axinset.get_ylim()[0], axinset.get_ylim()[1] * 1.2)
        axinset.axvline(x = 0.0, color = 'red', zorder = 1, alpha = 0.3)

        axinset.axvline(x = np.median(dzozs[d]), color = 'grey', linestyle = '-',zorder = 1, alpha = 0.3)
        axinset.axvline(x = np.median(dzozs[d]) -  nmads[d], color = 'grey', linestyle = '--',zorder = 1, alpha = 0.3)
        axinset.axvline(x = np.median(dzozs[d]) +  nmads[d], color = 'grey', linestyle = '--',zorder = 1, alpha = 0.3)


        if d == 0:
            if f == 0: zcheck = np.linspace(0,2.6, 1000)
            else: zcheck = np.linspace(0,5.0, 1000)
            zcheck2 = np.linspace(0,5.0, 1000)
            axes[d, f].plot(lg(zcheck2),  lg((1+zcheck2)*6563./5007. - 1.),'b-.', alpha = 0.55)
            axes[d, f].plot(lg(zcheck),   lg((1+zcheck)*5007./6563. - 1.),'b-.', alpha = 0.55)

            '''
            zcheck = np.linspace(0,0.5, 1000)
            zcheck2 = np.linspace(0,5, 1000)
            axes[d, f].plot(zcheck2, (1+zcheck2)*6563./5007. - 1.,'b--', alpha = 0.55)
            axes[d, f].plot(zcheck, (1+zcheck)*5007./6563. - 1.,'b--', alpha = 0.55)
            '''
            if f == 0:
                zline  = 2.5
                zother = (1+zline)*6563./5007. - 1.
                
                axes[d, f].annotate(r"[OIII]/H$\alpha$ confusion", (lg(zother),lg(zline)), xycoords = 'data', ha = 'center', va = 'center',color = 'blue', fontweight = 'bold', alpha = 0.55, rotation = 45, fontsize = 10)



for d in np.arange(2):
    axes[0,d].set_xticklabels([])
    axes[0,d].set_xlabel('')

fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/survey_paper/z_versus_z_%s.png'%postpend, dpi = 300)


#number of sources
#try a heatmap
#degeneracy lines






