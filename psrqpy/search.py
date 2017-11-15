"""
Search query
"""

from __future__ import print_function

import warnings
import six

import requests
from bs4 import BeautifulSoup

from .config import *

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


PSR_ALL = PSR_GENERAL + PSR_TIMING + PSR_BINARY + PSR_DERIVED

class QueryATNF(object):
    """
    Class to generate a query of the ATNF catalogue
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None, sort_attr='jname', sort_order='asc', psrs=None, **kwargs):
        """
        Set up and perform the query of the ATNF catalogue

        :param params: a list of strings with the pulsar parameters to return
        :param condition: a string with conditions for the returned parameters
        :param psrtype: a list of strings, or single string, of conditions on the 'type' of pulsars to return (logical AND will be used for any listed types)
        :param assoc: a condition on the associations of pulsars to return ()
        :param sort_attr: the parameter on which with sort the returned pulsars
        :param sort_ord: the order of the sorting (defaults to ascending)
        :param psrs: a list of pulsar names to get the information for
        """

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
            if p not in PSR_ALL:
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

        self.query_url = ATNF_URL + pquery + CONDITION_QUERY.format(conditionsquery) + \
                         PSRNAMES_QUERY.format('') + SORT_QUERY.format(**sortquery) + QUERY_FLUFF

        print(self.query_url)

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

    def table(self):
        """
        Return an astropy table of the pulsar data
        """
        
        from astropy.table import Table

        