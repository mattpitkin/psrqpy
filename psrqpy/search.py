"""
Search query
"""

from __future__ import print_function

import warnings
from collections import OrderedDict
import six

import numpy as np
import requests
from bs4 import BeautifulSoup

from .config import *
from .utils import *

class pulsar(object):
    """
    An object to hold a single pulsar
    """

    def __init__(self, **kwargs):
        """
        Set object attributes from kwargs
        """
        self._raw = kwargs
        for key, value in six.iteritems(kwargs):
            setattr(self, key, value)

    def keys(self):
        return self._raw.keys()

    def items(self):
        return self._raw.items()

class QueryATNF(object):
    """
    Class to generate a query of the ATNF catalogue
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None, sort_attr='jname', sort_order='asc', psrs=None,
                 include_errs=True, include_refs=False, version=None, adsref=False, **kwargs):
        """
        Set up and perform the query of the ATNF catalogue

        :param params: a list of strings with the pulsar parameters to return
        :param condition: a string with conditions for the returned parameters
        :param psrtype: a list of strings, or single string, of conditions on the 'type' of pulsars to return (logical AND will be used for any listed types)
        :param assoc: a condition on the associations of pulsars to return ()
        :param sort_attr: the parameter on which with sort the returned pulsars
        :param sort_ord: the order of the sorting (defaults to ascending)
        :param psrs: a list of pulsar names to get the information for
        :param include_errs: boolean to set whether to include parameter errors
        :param include_refs: boolean to set whether to include parameter references
        :param version: a string with the ATNF version to use (this will default to the current version if set as None)
        :param adsref: boolean to set whether the python 'ads' module can be used to get reference information
        """

        self._include_errs = include_errs
        self._include_refs = include_refs
        self._atnf_version = version
        self._adsref = adsref
        
        self._psr = psrs

        self._references = None # set of pulsar references

        # check parameters are allowed values
        if isinstance(params, list):
            if len(params) == 0:
                raise Exception("No parameters in list")

            for p in params:
                if not isinstance(p, str):
                    raise Exception("Non-string value '{}' found in params list".format(p))

            self._query_params = params
        else:
            if isinstance(params, str):
                self._query_params = [params]
            else:
                raise Exception("'params' must be a list or string")

        for p in list(self._query_params):
            if p not in PSR_ALL_PARS:
                warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                self._query_params.remove(p)
        if len(p) == 0:
            raise Exception("No parameters left in list")

        # set parameter string for query
        pquery = ''
        for p in self._query_params:
            pquery += '&{}={}'.format(p, p)

        # set conditions
        conditionsquery = ''
        if condition is not None:
            conditionsquery = condition

        if psrtype is not None:
            if isinstance(psrtype, list):
                if len(psrtype) == 0:
                    raise Exception("No pulsar types in list")

                for p in psrtype:
                    if not isinstance(p, str):
                        raise Exception("Non-string value '{}' found in params list".format(p))
                self._query_psr_types = psrtype
            else:
                if isinstance(psrtype, str):
                    self._query_psr_types = [psrtype]
                else:
                    raise Exception("'psrtype' must be a list or string")

            for p in list(self._query_psr_types):
                if p.upper() not in PSR_TYPES:
                    warnings.warn("Pulsar type '{}' is not recognised, no type will be required")
                    self._query_psr_types.remove(p)
                else:
                    if len(conditionsquery) == 0:
                        conditionsquery = 'type({})'.format(p.upper())
                    else:
                        conditionsquery += '+&&+type({})'.format(p.upper())

        sortquery = {'sortattr': sort_attr, 'sortorder': sort_order}

        self.query_url = ATNF_URL.format(self.get_version()) + pquery + CONDITION_QUERY.format(conditionsquery) + \
                         PSRNAMES_QUERY.format('') + SORT_QUERY.format(**sortquery) + QUERY_FLUFF

        print(self.query_url)

        # get references is required
        if self._include_refs:
            self._refs = get_references(useads=self._adsref)

        # generate request
        psrrequest = requests.get(self.query_url)

        if psrrequest.status_code != 200:
            raise Exception('Error... their was a problem with the request: status code {}'.format(psrrequest.status_code))

        # parse with BeautifulSoup
        try:
            psrsoup = BeautifulSoup(psrrequest.content, 'html.parser')
        except:
            raise Exception('Error... problem parsing catalogue with BeautifulSoup')

        self.query_output = psrsoup.find('pre').text

        # put the data in an ordered dictionary dictionary
        self.query_output_dict = OrderedDict()
        self.npulsars = 0
        if len(self.query_output) > 0:
            plist = self.query_output.strip().split('\n') # split output string
            
            self.npulsars = len(plist)

            for p in self._query_params:
                if p in PSR_ALL_PARS:
                    self.query_output_dict[p] = np.zeros(self.npulsars, dtype=PSR_ALL[p]['format'])
                    
                    if PSR_ALL[p]['err'] and self._include_errs:
                        self.query_output_dict[p+'_ERR'] = np.zeros(self.npulsars, dtype=PSR_ALL[p]['format'])
                    
                    if PSR_ALL[p]['ref'] and self._include_refs:
                        self.query_output_dict[p+'_REF'] = np.zeros(self.npulsars, dtype='S1024')
                

    def table(self):
        """
        Return an astropy table of the pulsar data
        """
        
        from astropy.table import Table

    def get_version(self):
        """
        Return a string with the ATNF version number, or the default giving in ATNF_VERSION if not found
        """

        if self._atnf_version is None:
            self._atnf_version = get_version()
        
        return self._atnf_version
