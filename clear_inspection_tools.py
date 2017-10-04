"""
CLEAR grizli inspection tools.

"""

import matplotlib.pyplot as plt


def flt_residuals(grp):
	"""
	"""
	for i in range(len(grp.FLTs)):
	    fig, ax = plt.subplots(figsize=[12,12])
	    ax.imshow(grp.FLTs[i].grism['SCI'] - grp.FLTs[i].model, vmin=-0.05, vmax=0.001,cmap='Greys',
	          interpolation='Nearest', origin='lower')
	    ax.set_title('%s' %(grp.FLTs[i].grism.parent_file))
	    fig.savefig('{}_residuals.png'.format((grp.FLTs[i].grism.parent_file).split('_flt.fits')[0]))
	    plt.close()