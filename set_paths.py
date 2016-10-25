

#! /usr/bin/env python

"""Sets the paths for the CLEAR project, replacing all those that are 
hardcoded in the following repos:
    clear
    mywfc3
    threedhst
    unicorn


Example:
    To get a path, import the dictionary and get key.

    from set_paths import paths
    path_to_outputs = paths['path_to_outputs']

Authors:

    C.M. Gosmeyer, Feb. 2016

Outputs:

"""

paths = {'path_to_outputs' : '/astro/clear/cgosmeyer/outputs/', \
         'path_to_RAW' : '/astro/clear/cgosmeyer/RAW/', \
         'path_to_ref_files' : '/astro/clear/cgosmeyer/ref_files/', \
         'path_to_PREPARE' : '/astro/clear/cgosmeyer/PREPARE/', \
         'path_to_software' : '/astro/clear/cgosmeyer/software/',
         'path_to_PERSIST' : '/astro/clear/cgosmeyer/PERSIST/',
         'path_to_Extractions' : '/astro/clear/cgosmeyer/Extractions/'}

         # path_to_ref_files contains REF, CONF, Synphot, iref, jref, and templates 
         # ref_files used to be Work, and REF used to be its own directory, not
         # a sub-directory.



