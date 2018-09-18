#! /usr/bin/env python

"""
Script to check that all parameters can be queried correctly.
"""

from psrqpy import QueryATNF
from psrqpy.config import *

# query all parameters for one pulsar (the Crab) one at a time via the webform
query = QueryATNF(psrs=['J0534+2200'], include_refs=True)
t = query.as_table

for i, p in enumerate(PSR_ALL_PARS):
    print('Parameter {}'.format(p))
 
    # check error and reference exist if expected
    if PSR_ALL[p]['err']:
        try:
            err = t[p+'_ERR']
        except IOError:
            raise IOError('No error value found for parameter "{}"'.format(p))

    if PSR_ALL[p]['ref']:
        try:
            ref = t[p+'_REF']
        except IOError:
            raise IOError('No reference found for parameter "{}"'.format(p))
