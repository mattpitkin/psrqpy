#! /usr/bin/env python

"""
Script to check that all parameters can be queried correctly.
"""

from psrqpy import QueryATNF
from psrqpy.config import *

# query all parameters for one pulsar (the Crab) one at a time
for i, p in enumerate(PSR_ALL_PARS):
    print('Parameter: "{}"'.format(p))
    if i == 0:
        query = QueryATNF(params=p, psrs=['J0534+2200'], include_refs=True)
    else:
        # just call class methods rather than creating a new class
        query.generate_query(params=p)
        query.parse_query()
    t = query.table()

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
