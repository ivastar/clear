import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob
import clear.plots

def get_footprints():

    from mastquery import query, overlaps

    #CLEAR
    parent = query.run_query(box=None, proposal_id=[14227],
        instruments=['WFC3/IR'],filters=['G102'])
    tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01,filters=['G102', 'G141'],
        instruments=['WFC3/IR'],close=False)

    #CANDELS F160W
    parent = query.run_query(box=[53.12916667,-27.815, 20.], proposal_id = [12061, 12062, 11563, 12099, 11359, 12498],
        instruments=['WFC3/IR'],filters=['F160W'])
    parent.write('goodss_candels.fits', format='fits')
    parent = query.run_query(box = [189.20833, 62.23527, 20.], proposal_id = [12445, 12444, 12443,13063],
        instruments=['WFC3/IR'],filters=['F160W'])
    parent.write('goodsn_candels.fits', format='fits')

def clear_footprints_goodsn():

    from clear.plots import polysplit
    from matplotlib import gridspec

    plt.rcParams['lines.linewidth'] = 0.3
    fontsize = 9
    left = 0.32
    right = 0.02
    top = 0.01
    bottom = 0.115
    width = 10
    corner = 'ur'


    g102_clear_color = 'red'
    g102_other_color = 'orange'
    g141_color = 'blue'

    x1, x0 = 188.8639017749491, 189.60688055895648
    y0, y1 = 62.063791549511998, 62.394068625281309
    xticklab = [r'$12^\mathrm{h}37^\mathrm{m}30^\mathrm{s}$', r'$37^\mathrm{m}00^\mathrm{s}$', r'$36^\mathrm{m}30^\mathrm{s}$', r'$36^\mathrm{m}00^\mathrm{s}$']
    xtickv = SkyCoord(ra=['12 37 30', '12 37 00', '12 36 30', '12 36 00'], dec = ['0 0 0', '0 0 0', '0 0 0', '0 0 0'], unit=(u.hourangle, u.deg)).ra.degree
    yticklab = [r'$+62^\circ10^\prime00^{\prime\prime}$', r'$15^\prime00^{\prime\prime}$', r'$20^\prime00^{\prime\prime}$']
    ytickv = SkyCoord(ra = ['0 0 0', '0 0 0', '0 0 0'], dec=['62 10 00','62 15 00', '62 20 00'], unit=(u.hourangle, u.deg)).dec.degree
    clear_files = glob.glob('/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/v2.1/j12*footprint.fits')

    field_jname = {'GN1':'j123716p6222','GN2':'j123652p6215','GN3':'j123700p6219','GN4':'j123712p6216','GN5':'j123640p6219','GN7':'j123632p6213'}

    dx = np.abs(x1-x0)*np.cos(y0/360*2*np.pi)
    dy = (y1-y0)

    fig = plot_init(square=True, xs=width, aspect=dy/dx, fontsize=fontsize, left=left,
        right=right, top=top, bottom=bottom)

    gs = fig.add_gridspec(4,4, wspace=0.05,hspace=0.05)
    ax = fig.add_subplot(gs[1:,0:3])

    g102_clear_polys = []
    g102_other_polys = []
    g141_polys = []

    for file in clear_files:

        tab = Table.read(file)

        for poly in tab[(tab['proposal_id'] != '14227') & (tab['filter'] == 'G102')]['footprint']:
            gx, gy = polysplit(poly)
            ax.fill(gx, gy, alpha=0.05, color=g102_clear_color)
            ax.plot(gx, gy, alpha=0.5, color=g102_clear_color, linestyle=(0, (5, 10)), linewidth = 0.5)
            #g102_other_polys.append(gx, gy)

        for poly in tab[tab['filter'] == 'G141']['footprint']:
            gx, gy = polysplit(poly)
            ax.fill(gx, gy, alpha=0.05, color=g141_color)
            ax.plot(gx, gy, alpha=0.5, color=g141_color)
            #g141_polys.append(gx, gy)

        for poly in tab[tab['proposal_id'] == '14227']['footprint']:
            gx, gy = polysplit(region=poly)
            ax.fill(gx, gy, alpha=0.1, color=g102_clear_color)
            ax.plot(gx, gy, alpha=0.5, color=g102_clear_color)
            #g102_clear_polys.append(gx, gy)

        for label in np.unique(tab[tab['proposal_id'] == '14227']['target']):
            x = np.mean(tab[tab['target'] == label]['ra'])
            y = np.mean(tab[tab['target'] == label]['dec'])
            if x < 0:
                x = 360+x
            ax.text(x,y,label,va='center', ha='center', fontsize=9)

        ll = ['FIGS-GN1','COLFAX','GN-Z10-1']
        offset = [-0.001, 0.001, 0.0]
        for label,off in zip(ll, offset):
            if label in tab['target']:
                x = np.mean(tab[tab['target'] == label]['ra'])
                y = np.mean(tab[tab['target'] == label]['dec'])
                if x < 0:
                    x = 360+x
                #ax.text(x,y+off,label,va='center', ha='center', fontsize=7)



    candels_tab = Table.read('/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/goodsn_candels.fits')
    candels_polys = []
    for poly in candels_tab['footprint']:
        candels_polys.append(polysplit(poly, get_shapely=True))

    candels_polys_union = candels_polys[0]
    for pp in candels_polys[1:]:
        candels_polys_union = candels_polys_union.union(pp)

    cx,cy = candels_polys_union.exterior.xy
    ax.plot(cx,cy, alpha=0.3, color='0.1', linewidth=1)
    ax.fill(cx,cy, alpha=0.1, color='0.7')


    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')

    fsi = '20'
    field_label = 'GOODS-N'

    ax.text(0.95, 0.95,r'$\mathit{%s}$' %(field_label),
    horizontalalignment='right',
    verticalalignment='top',
    transform = ax.transAxes, fontsize=fsi)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    for field, gs_sub in zip(field_jname.keys(),[gs[0,0],gs[0,1],gs[0,2],gs[1,3],gs[2,3],gs[3,3]]):
        tab_sub = Table.read(f'/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/v2.1/{field_jname[field]}_footprint.fits')
        ax_sub = fig.add_subplot(gs_sub)
        ax_sub.axes.get_xaxis().set_ticks([])
        ax_sub.axes.xaxis.set_ticklabels([])
        ax_sub.axes.get_yaxis().set_ticks([])
        ax_sub.axes.yaxis.set_ticklabels([])
        ax_sub.invert_xaxis()
        for poly in tab_sub:
            gx_sub, gy_sub = polysplit(poly['footprint'].replace('ICRS',''))
            if ((poly['proposal_id'] == '14227') and (poly['filter'] == 'G102')):
                ax_sub.fill(gx_sub, gy_sub, alpha=0.1, color=g102_clear_color)
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g102_clear_color)
            if poly['filter'] == 'G141':
                ax_sub.fill(gx_sub, gy_sub, alpha=0.05, color=g141_color)
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g141_color)
            if ((poly['proposal_id'] != '14227') and (poly['filter'] == 'G102')):
                ax_sub.fill(gx_sub, gy_sub, alpha=0.05, color=g102_clear_color, linestyle='dotted')
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g102_clear_color, linestyle=(0, (5, 10)), linewidth = 0.5)
        for label in np.unique(tab_sub['target']):
            x = np.mean(tab_sub[tab_sub['target'] == label]['ra'])
            y = np.mean(tab_sub[tab_sub['target'] == label]['dec'])
            label_text = label.replace('GNGRISM','ANGST-').replace('GDN','BARRO-')
            if x < 0:
                x = 360+x
            if label_text.startswith('ANGST'):
                y +=0.009
            if label_text == 'ANGST-18':
                y += 0.003
            if label_text == 'BARRO-26':
                y -= 0.003
            if label_text == 'COLFAX+WHEELER':
                y -= 0.009
                x += 0.014
            if label_text == 'COLFAX+HAMLIN':
                y -=0.012
            if label_text == 'BARRO-3':
                y += 0.005
            if label_text == 'BARRO-7':
                y -= 0.005
            if (('14227' in tab_sub[tab_sub['target'] == label]['proposal_id']) and ('G102' in tab_sub[tab_sub['target'] == label]['filter'])):
                label_color = 'black'
            if 'G141' in tab_sub[tab_sub['target'] == label]['filter']:
                label_color = g141_color
            if (('14227' not in tab_sub[tab_sub['target'] == label]['proposal_id']) and ('G102' in tab_sub[tab_sub['target'] == label]['filter'])):
                label_color = g102_clear_color
            ax_sub.text(x,y,label_text,va='center', ha='center', fontsize=7.5, color=label_color)
        ax_sub.set_aspect(2.13,'datalim')

    #plt.show(block=False)

    plt.savefig('/Users/imomcheva/WORK/CLEAR/PLOTS/goodsn_clear_footprints.pdf')

def clear_footprints_goodss():

    from matplotlib import gridspec

    plt.rcParams['lines.linewidth'] = 0.3
    fontsize = 10
    left = 0.32
    width = 10
    corner = 'll'

    g102_clear_color = 'red'
    g141_color = 'blue'

    x0, x1 = 53.324005633802822,  52.966197183098595
    y0, y1 = -27.983151830808076, -27.644474431818176
    xticklab = [r'$3^\mathrm{h}33^\mathrm{m}00^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}30^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}00^\mathrm{s}$']
    xtickv = SkyCoord(ra=['3 33 00', '3 32 30', '3 32 00'], dec = ['0 0 0', '0 0 0', '0 0 0'], unit=(u.hourangle, u.deg)).ra.degree
    yticklab = [r'$-27^\circ40^\prime00^{\prime\prime}$', r'$45^\prime00^{\prime\prime}$', r'$-27^\circ50^\prime00^{\prime\prime}$', r'$55^\prime00^{\prime\prime}$']
    ytickv = SkyCoord(ra = ['0 0 0', '0 0 0', '0 0 0', '0 0 0'], dec=['-27 40 00','-27 45 00', '-27 50 00', '-27 55 00'], unit=(u.hourangle, u.deg)).dec.degree
    clear_files = glob.glob('/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/v2.1/j0*footprint.fits')

    field_jname = {'GS1':'j033300m2742','GS2':'j033232m2742','GS3':'j033236m2744','GS4':'j033236m2748','GS5':'j033228m2743','ERSPRIME':'j033216m2743'}

    dx = np.abs(x1-x0)*np.cos(y0/360*2*np.pi)
    dy = (y1-y0)

    fig = plot_init(square=True, xs=width, aspect=dy/dx, fontsize=fontsize, left=left)

    gs = fig.add_gridspec(4,4, wspace=0.05,hspace=0.05)
    ax = fig.add_subplot(gs[1:,0:3])

    g102_clear_polys = []
    g102_other_polys = []
    g141_polys = []

    for file in clear_files:

        tab = Table.read(file)

        for poly in tab[(tab['proposal_id'] != '14227') & (tab['filter'] == 'G102')]['footprint']:
            gx, gy = polysplit(poly)
            ax.fill(gx, gy, alpha=0.05, color=g102_clear_color)
            ax.plot(gx, gy, alpha=0.5, color=g102_clear_color, linestyle=(0, (5, 10)), linewidth = 0.5)
            #g102_other_polys.append(gx, gy)

        for poly in tab[tab['filter'] == 'G141']['footprint']:
            gx, gy = polysplit(poly)
            ax.fill(gx, gy, alpha=0.05, color=g141_color)
            ax.plot(gx, gy, alpha=0.5, color=g141_color)
            #g141_polys.append(gx, gy)

        for poly in tab[tab['proposal_id'] == '14227']['footprint']:
            gx, gy = polysplit(region=poly)
            ax.fill(gx, gy, alpha=0.08, color=g102_clear_color)
            ax.plot(gx, gy, alpha=0.5, color=g102_clear_color)
            #g102_clear_polys.append(gx, gy)

        for label in np.unique(tab[tab['proposal_id'] == '14227']['target']):
            x = np.mean(tab[tab['target'] == label]['ra'])
            y = np.mean(tab[tab['target'] == label]['dec'])
            if x < 0:
                x = 360+x
            ax.text(x,y,label,va='center', ha='center', fontsize=9)

        ll = ['FIGS-GS1','PRIMO','GOODS-SOUTH-36']
        offset = [-0.000, 0.005, -0.003]
        for label,off in zip(ll, offset):
            if label in tab['target']:
                x = np.mean(tab[tab['target'] == label]['ra'])
                y = np.mean(tab[tab['target'] == label]['dec'])
                if x < 0:
                    x = 360+x
                #if label == 'GOODS-SOUTH-36':
                    #ax.text(x,y+off,'HUDF',va='center', ha='center', fontsize=7)
                #else:
                    #ax.text(x,y+off,label,va='center', ha='center', fontsize=7)


    candels_tab = Table.read('/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/goodss_candels.fits')
    candels_polys = []
    for poly in candels_tab['footprint']:
        candels_polys.append(polysplit(poly, get_shapely=True))

    candels_polys_union = candels_polys[0]
    for pp in candels_polys[1:]:
        candels_polys_union = candels_polys_union.union(pp)

    for sub_poly in candels_polys_union.geoms:
        cx,cy = sub_poly.exterior.xy
        ax.plot(cx,cy, alpha=0.3, color='0.1', linewidth=1)
        ax.fill(cx,cy, alpha=0.1, color='0.7')


    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')

    fsi = '20'
    field_label = 'GOODS-S'

    ax.text(0.05, 0.05,r'$\mathit{%s}$' %(field_label),
        horizontalalignment='left',
        verticalalignment='bottom',
        transform = ax.transAxes, fontsize=fsi)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    for field, gs_sub in zip(field_jname.keys(),[gs[0,0],gs[0,1],gs[0,2],gs[1,3],gs[2,3],gs[3,3]]):
        tab_sub = Table.read(f'/Users/imomcheva/WORK/CLEAR/CLEAR/FOOTPRINTS/v2.1/{field_jname[field]}_footprint.fits')
        ax_sub = fig.add_subplot(gs_sub)
        ax_sub.axes.get_xaxis().set_ticks([])
        ax_sub.axes.xaxis.set_ticklabels([])
        ax_sub.axes.get_yaxis().set_ticks([])
        ax_sub.axes.yaxis.set_ticklabels([])
        ax_sub.invert_xaxis()
        for poly in tab_sub:
            gx_sub, gy_sub = polysplit(poly['footprint'].replace('ICRS',''))
            if poly['filter'] == 'G141':
                ax_sub.fill(gx_sub, gy_sub, alpha=0.05, color=g141_color)
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g141_color)
            if ((poly['proposal_id'] != '14227') and (poly['filter'] == 'G102')):
                ax_sub.fill(gx_sub, gy_sub, alpha=0.01, color=g102_clear_color, linestyle='dotted')
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g102_clear_color, linestyle=(0, (5, 10)), linewidth = 0.5)
            if ((poly['proposal_id'] == '14227') and (poly['filter'] == 'G102')):
                ax_sub.fill(gx_sub, gy_sub, alpha=0.1, color=g102_clear_color)
                ax_sub.plot(gx_sub, gy_sub, alpha=0.5, color=g102_clear_color)
        for label in np.unique(tab_sub['target']):
            if label in ['ANY','BUSH1','BUSH3','WFC3-ERSII-IR03','WFC3-ERSII-IR04','WFC3-ERSII-IR09']:
                continue
            else:
                x = np.mean(tab_sub[tab_sub['target'] == label]['ra'])
                y = np.mean(tab_sub[tab_sub['target'] == label]['dec'])
                label_text = label.replace('GOODS-SOUTH','3DHST').replace('GDN','BARRO-')
                if x < 0:
                    x = 360+x
                if label_text == 'PRIMO':
                    y += 0.009
                if label_text == 'FIGS-GS1':
                    y += 0.015
                    x += 0.02
                if label_text == '3DHST-32':
                    label_text = '3DHST-32,34,36,37,38'
                    x += 0.009
                if label_text in ['3DHST-34','3DHST-36','3DHST-37','3DHST-38']:
                    label_text = ''
                if label_text == '3DHST-30':
                    y -= 0.009
                #if label_text.startswith('ANGST'):
                #    y +=0.009
                if (('14227' in tab_sub[tab_sub['target'] == label]['proposal_id']) and ('G102' in tab_sub[tab_sub['target'] == label]['filter'])):
                    label_color = 'black'
                if 'G141' in tab_sub[tab_sub['target'] == label]['filter']:
                    label_color = g141_color
                if (('14227' not in tab_sub[tab_sub['target'] == label]['proposal_id']) and ('G102' in tab_sub[tab_sub['target'] == label]['filter'])):
                    label_color = g102_clear_color
                ax_sub.text(x,y,label_text,va='center', ha='center', fontsize=7.5, color=label_color)
        ax_sub.set_aspect(1.122,'datalim')

    #plt.show(block=False)

    plt.savefig('/Users/imomcheva/WORK/CLEAR/PLOTS/goodss_clear_footprints.pdf')


def grisms_coverage():

    pass

def plot_init(square=True, xs=6, aspect=1, left=0.22, bottom=0.11, \
    right=0.02, top=0.02, wspace=0.2, hspace=0.02, fontsize=10, NO_GUI=False, \
    use_tex=False, invert=False):
    """
    Wrapper for generating a plot window, contains input parameters for setting the
    full window geometry and also handles toggling the GUI/interactive backend.

    NO_GUI should be set to True if your session has no X11 connection.
    """
    import matplotlib
    rc = matplotlib.rcParams

    plt.rcParams['patch.edgecolor'] = 'None'
    plt.rcParams['font.size'] = fontsize

    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['image.interpolation'] = 'nearest'

    if use_tex:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'

    rc['lines.color'] = 'black'
    rc['patch.edgecolor'] = 'black'
    rc['text.color'] = 'black'
    rc['axes.facecolor'] = 'white'
    rc['axes.edgecolor'] = 'black'
    rc['axes.labelcolor'] = 'black'
    rc['xtick.color'] = 'black'
    rc['ytick.color'] = 'black'
    rc['grid.color'] = 'black'
    rc['figure.facecolor'] = 'white'
    rc['figure.edgecolor'] = 'white'
    rc['savefig.facecolor'] = 'white'
    rc['savefig.edgecolor'] = 'white'

    if square:
        lrbt = np.array([left,right,bottom,top])*5./xs
        ys = (1-lrbt[1]-lrbt[0])/(1-lrbt[3]-lrbt[2])*xs*aspect
        lrbt[[2,3]] /= aspect

        fig = plt.figure(figsize=(xs,ys), dpi=100)
        fig.subplots_adjust(left=lrbt[0], bottom=lrbt[2], right=1-lrbt[1],
            top=1-lrbt[3], wspace=wspace, hspace=hspace)

    else:
        fig = plt.figure(figsize=(7,5), dpi=100)
        fig.subplots_adjust(wspace=wspace, hspace=hspace,left=0.10,
                        bottom=0.10,right=0.99,top=0.97)

    return fig
def throughput():

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    import pysynphot as S

    bp = S.ObsBandpass('wfc3,ir,g141')
    xg141, yg141 = bp.wave, bp.throughput
    bp = S.ObsBandpass('wfc3,ir,g102')
    xg102, yg102 = bp.wave, bp.throughput
    bp = S.ObsBandpass('wfc3,ir,f140w')
    xf140, yf140 = bp.wave, bp.throughput
    bp = S.ObsBandpass('wfc3,ir,f105w')
    xf105, yf105 = bp.wave, bp.throughput

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'Serif'
    plt.rcParams['font.serif'] = 'Times'

    plt.ioff()
    fig = plot_init(square=True, xs=8, aspect=1./2.5, left=0.15, bottom=0.1, top=0.095, right=0.01)
    ax = fig.add_subplot(111)

    ax.plot(xg141, yg141, color='black', linewidth=2, alpha=0.5)
    ax.fill(xg141, yg141, color='red', linewidth=2, alpha=0.1)
    ax.plot(xg102, yg102, color='blue', linewidth=2, alpha=0.5)
    ax.fill(xg102, yg102, color='blue', linewidth=2, alpha=0.1)
    ax.plot(xf140, yf140, color='black', linewidth=2, alpha=0.7)
    ax.plot(xf105, yf105, color='black', linewidth=2, alpha=0.7)


    dlam = 30
    zi = 1

    yy = 0.1

    ax.text(1.35e4, 0.29+yy,'G141',rotation=0., color='black', alpha=0.7)
    ax.text(1.35e4, 0.29+yy,'G141',rotation=0., color='red', alpha=0.4)

    ax.text(9.5e3, 0.3,'G102',rotation=17., color='black', alpha=0.7)
    ax.text(9.5e3, 0.3,'G102',rotation=17., color='orange', alpha=0.4)

    ax.text(1.4e4, 0.57,'F140W',color='black', alpha=0.9)
    ax.text(1.01e4, 0.52,'F105W',rotation=10., color='black', alpha=0.9)

    ax.set_xlim(7500, 1.79e4)
    ax.set_ylim(0,0.65)
    ytick = ax.set_yticks([0,0.2,0.4,0.6])
    ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_ylabel('Throughput')
    minorLocator   = MultipleLocator(1000)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator   = MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(minorLocator)

    #fig.savefig('throughput.eps')
    fig.savefig('/Users/imomcheva/WORK/CLEAR/PLOTS/throughput.pdf')

def grisms_coverage_lines():

    lines = [(r'Ly$\alpha$-1215', 1215.4),('NV-1240', 1240.81),('NIV-1487', 1487.0),('CIV-1549', 1549.48),
        ('HeII-1640', 1640.4),
        ('OIII-1663', 1665.85),
        ('NIII-1750', 1750.0),
        ('CIII-1908', 1908.734),
        ('MgII-2799', 2799.117),
        ('NeV-3346', 3346.8),
        ('NeVI-3426', 3426.85),
        ('OII-3727,3729', 3727.092),
        ('NeIII-3867', 3867.5),
        ('H9-3890', 3890.166),
        ('H8-3971', 3971.198),
        ('Hd-4102', 4102.892),
        ('Hg-4341', 4341.692),
        ('OIII-4363', 4364.436),
        (r'H$\beta$-4862', 4862.71),
        ('OIII-4960,5008', 5008.24),
        ('HeI-5877', 5877.2),
        ('OI-6302, 6363', 6302.046),
        (r'H$\alpha$+NII-6563', 6564.61),
        ('SII-6718,6732', 6718.29),
        ('ArIII-7138', 7138.0),
        ('OII-7325', 7322.0),
        ('SIII-9531,9069', 9530.6),
        ('HeI-10830',10830.0),('PaB-12822',12821.6)]

    #lines = lines[::-1]

    filters = [('G102', 8000, 11500, 'blue'),
        ('G141', 10750, 17000, 'red')]

    filters_overlap = [(10750, 11500, 'purple')]

    #plt.rcParams['text.usetex'] = True
    fig = plot_init(square=True, fontsize=10, xs=7)
    ax = fig.add_subplot(111)
    #ax.set_xlabel('z_line', fontsize = 35)​
    ax.fill_between([0,10],[10.5,10.5],[23.5,23.5], color='gray', alpha=0.1)
    ax.annotate('Rest-frame optical', (0.65, 0.75), xycoords = 'axes fraction', color = 'gray', ha = 'left', va = 'center')

    for f, filt in enumerate(filters[::-1]):
        filt_name, filt_lam0, filt_lam1, txt_color  = filt
        ax.annotate(filt_name, (0.75-0.1*f, 0.9), xycoords = 'axes fraction', color = txt_color, alpha=0.5)
        for l, line in enumerate(lines):
            line_lam_rest = line[1]
            z_min_filt = filt_lam0/line_lam_rest - 1.
            z_max_filt = filt_lam1/line_lam_rest - 1.
            xmn = np.log10((1+z_min_filt))
            xmx = np.log10((1+z_max_filt))
            ymn = (l + 0.6)/(len(lines)+1)
            ymx = (l + 1.4)/(len(lines)+1)
            ax.annotate(u'{} \u212B'.format(line[0].split('-')[1]), (0.015, (ymn+ymx)/2), ha = 'left', va = 'center', xycoords = 'axes fraction', color = 'grey')
            ax.axvspan(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, facecolor = txt_color, edgecolor = None, alpha=0.3)

    # y labels and ticks
    ytick_labels = np.array(['{:>12}'.format(line[0].split('-')[0])  for line in lines])
    yticks = np.arange(len(lines))

    ax.set_yticks(yticks)
    ax.set_ylim(-1, len(lines))
    ax.set_yticklabels(ytick_labels)

    ztck_min = 0.
    ztck_max = 9.
    ax.set_xlim(np.log10(1+ztck_min), np.log10(1.+ztck_max))
    zticks = np.arange(ztck_min, ztck_max)
    zticks_str = ['%i'%z for z in zticks]
    zticks_log10_p1 = np.log10(zticks + 1)
    ax.set_xticks(zticks_log10_p1)
    ax.set_xticklabels(zticks_str)

    plt.savefig('/Users/imomcheva/WORK/CLEAR/PLOTS/grisms_coverage.pdf')
    #fig.subplots_adjust(left = 0.20, bottom = 0.15, top = 0.95, right = 0.98)

    #​plt.show(block = False)
    #fig.savefig('/Users/rsimons/Desktop/clear/figures/%s'%figname, dpi = 300)

def polysplit(region='POLYGON 53.168267 -27.72893 53.129988 -27.725711 53.129952122566124 -27.72608704794391 53.129821 -27.726076 53.126152 -27.76452 53.126399821840309 -27.764534699417919', get_shapely=False):

    from shapely.geometry import Polygon

    spl = region.split()

    px = spl[1::2]
    py = spl[2::2]
    px.append(px[0])
    py.append(py[0])
    px, py = np.cast[float](px), np.cast[float](py)
    if px[0] < 0:
        px = 360.+px

    if get_shapely:
        list = []
        for i in range(len(px)):
            list.append((px[i], py[i]))

        poly = Polygon(tuple(list))
        return poly
    else:
        return px, py
