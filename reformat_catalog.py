#! /usr/bin/env python

"""
Takes subset of new catalogs, GoodsN_plus.cat and GoodsS_plus.cat, and reformats them to
match the old SourceExtractor-formatted catalogs.

I will actually copy/paste only the new sources into two new files to run this
script over and then paste them into the old catalogs (which will be renamed
GoodsN_plus_2.cat and GoodsS_plus_2.cat).

Author: 

    C.M. Gosmeyer, 1 Nov. 2016
"""

from clear.set_paths import paths
from astropy.io import ascii
import numpy as np
import os



def asciiwrite_prep(column_lists, column_names=[]):
    """
    Parameters:
        column_lists : list of lists/arrays
            The columns to be inserted into file.
        column_names : list of strings
            The names of the columns. Optional. 
            By default will name columns 0 - len(column_lists).

    Returns:
        tab : dictionary
            Columns are values, column names are keys.
        column_names : list of strings
            Names of the columns.
    """
    tab = {}

    if column_names == []:
        column_names = np.arange(len(column_lists))
        column_names = [str(i) for i in column_names]
    #print len(column_lists)
    print column_names

    for column, column_name in zip(column_lists, column_names):
        tab[column_name] = column

    return tab, column_names
        


def reformat_catalog():
    """
    """
    path_to_ref = os.path.join(paths['path_to_ref_files'], 'REF')
    gn_truncated_cat = os.path.join(path_to_ref, 'GoodsN_plus_truncated.cat')
    gs_truncated_cat = os.path.join(path_to_ref, 'GoodsS_plus_truncated.cat')

    catalogs = [gn_truncated_cat, gs_truncated_cat]
    out_catalogs = [os.path.join(path_to_ref, 'GoodsN_plus_reformatted.cat'),
        os.path.join(path_to_ref, 'GoodsS_plus_reformatted.cat')]

    for cat, out_cat in zip(catalogs, out_catalogs):
        data = ascii.read(cat)
        numbers = np.array(data['NUMBER'])
        x_images = np.array(data['X_IMAGE'])
        y_images = np.array(data['Y_IMAGE'])
        x_worlds = np.array(data['X_WORLD'])
        y_worlds = np.array(data['Y_WORLD'])

        # Need fill out the MAG_AUTO and MAG_BEST columsn with 30.0
        # FLAGS with 1
        # And all other columsn with 0 or 99.
        nrows = len(numbers)
        mag_autos = np.array([30.0]*nrows)
        mag_bests = np.array([30.0]*nrows)
        flags = np.array([1]*nrows)
        arrays_0 = np.array([0.0]*nrows)
        arrays_99 = np.array([99.0]*nrows)

        column_lists = [numbers, 
                        x_images, 
                        y_images, 
                        arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, 
                        x_worlds,
                        y_worlds,
                        arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, 
                        arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, arrays_0, 
                        arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, 
                        arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, arrays_99, 
                        arrays_0, arrays_0, 
                        arrays_99, 
                        mag_autos,
                        arrays_0,
                        arrays_99, arrays_99,
                        mag_bests,
                        arrays_0, arrays_0, arrays_0, 
                        flags,
                        arrays_0]
        tab, column_names = asciiwrite_prep(column_lists=column_lists)

        ascii.write(tab, out_cat, names=column_names)



if __name__=="__main__":
    reformat_catalog()