"""
Search query
"""

from __future__ import print_function

import warnings
from collections import OrderedDict
import re
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

        self._psrs = psrs
        self._include_errs = include_errs
        self._include_refs = include_refs
        self._atnf_version = version
        self._atnf_version = self.get_version() # if no version is set this will return the current or default value
        self._adsref = adsref
        self._sort_order = sort_order
        self._sort_attr = sort_attr

        self._refs = None # set of pulsar references
        self._query_output = None

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

        # set conditions
        self._conditions_query = ''
        if condition is not None:
            self._conditions_query = condition

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
                    if not self._conditions_query:
                        self._conditions_query = 'type({})'.format(p.upper())
                    else:
                        self._conditions_query += '+&&+type({})'.format(p.upper())

        # get references is required
        if self._include_refs:
            self._refs = get_references(useads=self._adsref)

        # perform query
        self._query_content = self.generate_query()

        # parse the query with BeautifulSoup into a dictionary
        self._query_output = self.parse_query()

    def generate_query(self, version='', params=None, condition='', sortorder='asc', sortattr='JName', psrnames=None, getephemeris=False):
        """
        Generate a query URL and return the content of the request from that URL. If set the class attributes are
        used for generating the query, otherwise arguments can be given.

        :param version: a string containing the ATNF version
        :param params: a list of parameters to query
        :param condition: the condition string for the query
        :param sortorder: the order for sorting the results
        :param sortattr: the attribute on which to perform the sorting
        :param psrnames: a list of pulsar names to get
        """

        query_dict = {}
        self._atnf_version = self._atnf_version if not version else version
        query_dict['version'] = self._atnf_version

        if params:
            if isinstance(params, basestring):
                params = [params] # convert to list
            else:
                if not isinstance(params, list):
                    raise Exception('Error... input "params" for generate_query() must be a list')
            qparams = list(params)
            for p in params:
                if p not in PSR_ALL_PARS:
                    warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                    qparams.remove(p)
            self._query_params = qparams

        pquery = ''
        for p in self._query_params:
            pquery += '&{}={}'.format(p, p)

        query_dict['params'] = pquery
        self._conditions_query = self._conditions_query if not condition else condition
        query_dict['condition'] = self._conditions_query
        self._sort_order = self._sort_order if sortorder == self._sort_order else sortorder
        query_dict['sortorder'] = self._sort_order
        self._sort_attr = self._sort_attr if sortattr == self._sort_attr else sortattr
        query_dict['sortattr'] = self._sort_attr

        if psrnames:
            if isinstance(psrnames, basestring):
                self._psrs = [psrnames] # convert to list
            else:
                if not isinstance(psrnames, list):
                    raise Exception('Error... input "psrnames" for generate_query() must be a list')
                self._psrs = list(psrnames) # reset self._psrs

        qpulsars = '' # pulsar name query string
        if self._psrs is not None:
            for psr in self._psrs:
                if '+' in psr: # convert '+'s in pulsar names to '%2B' for the query string
                    qpulsars += psr.replace('+', '%2B')
                else:
                    qpulsars += psr
                qpulsars += '+' # seperator between pulsars
            qpulsars = qpulsars.strip('+') # remove the trailing '+'
        query_dict['psrnames'] = qpulsars

        # get pulsar ephemeris rather than table (parsing of this is not implemented yet)
        if getephemeris:
            query_dict['getephemeris'] = 'Get+Ephemeris'
        else:
            query_dict['getephemeris'] = ''

        # generate query URL
        self._query_url = QUERY_URL.format(**query_dict)

        # generate request
        psrrequest = requests.get(self._query_url)

        if psrrequest.status_code != 200:
            raise Exception('Error... their was a problem with the request: status code {}'.format(psrrequest.status_code))

        return psrrequest.content

    def parse_query(self, requestcontent=''):
        """
        Parse the query returned by requests

        :param requestcontent: the content of a request returned by the requests module
        """

        # update request if required
        self._query_content = requestcontent if requestcontent else self._query_content

        # parse through BeautifulSoup
        try:
            psrsoup = BeautifulSoup(self._query_content, 'html.parser')
        except:
            raise Exception('Error... problem parsing catalogue with BeautifulSoup')

        pretags = psrsoup.find_all('pre') # get any <pre> html tags

        # check for any warnings generated by the request
        for pt in pretags:
            if 'WARNING' in pt.text:
                warnings.warn('Request generated warning: "{}"'.format(pt.text), UserWarning)

        # actual table should be in the final <pre> tag
        qoutput = pretags[-1].text

        # put the data in an ordered dictionary dictionary
        self._query_output = OrderedDict()
        self._npulsars = 0
        if qoutput:
            plist = qoutput.strip().split('\n') # split output string

            self._npulsars = len(plist)

            for p in self._query_params:
                if p in PSR_ALL_PARS:
                    self._query_output[p] = np.zeros(self._npulsars, dtype=PSR_ALL[p]['format'])

                    if PSR_ALL[p]['err'] and self._include_errs:
                        self._query_output[p+'_ERR'] = np.zeros(self._npulsars, dtype='f8') # error can only be floats

                    if PSR_ALL[p]['ref'] and self._include_refs:
                        self._query_output[p+'_REF'] = np.zeros(self._npulsars, dtype='S1024')

                        if self._adsref: # also add reference URL for NASA ADS
                            self._query_output[p+'_REFURL'] = np.zeros(self._npulsars, dtype='S1024')

            for idx, line in enumerate(plist):
                pvals = line.split()[1:] # ignore the first entry as it is always in index

                vidx = 0 # index of current value
                for p in self._query_params:
                    if PSR_ALL[p]['format'] == 'f8':
                        if pvals[vidx] == '*':
                            self._query_output[p][idx] = None # put NaN entry in numpy array
                        else:
                            self._query_output[p][idx] = float(pvals[vidx])
                    elif PSR_ALL[p]['format'] == 'i4':
                        if pvals[vidx] == '*':
                            self._query_output[p][idx] = None
                        else:
                            self._query_output[p][idx] = int(pvals[vidx])
                    else:
                        self._query_output[p][idx] = pvals[vidx]
                    vidx += 1

                    # get errors
                    if PSR_ALL[p]['err']:
                        if self._include_errs:
                            if pvals[vidx] == '*':
                                self._query_output[p+'_ERR'][idx] = None
                            else:
                                self._query_output[p+'_ERR'][idx] = float(pvals[vidx])
                        vidx += 1

                    # get references
                    if PSR_ALL[p]['ref']:
                        if self._include_refs:
                            reftag = pvals[vidx]

                            if reftag in self._refs:
                                thisref = self._refs[reftag]
                                refstring = '{authorlist}, {year}, {title}, {journal}, {volume}'
                                refstring2 = re.sub(r'\s+', ' ', refstring.format(**thisref)) # remove any superfluous whitespace
                                self._query_output[p+'_REF'][idx] = ','.join([a for a in refstring2.split(',') if a.strip()]) # remove any superfluous empty ',' seperated values

                                if self._adsref and 'ADS URL' in thisref:
                                    self._query_output[p+'_REFURL'][idx] = thisref['ADS URL'] # remove any superfluous whitespace
                            else:
                                warnings.warn('Reference tag "{}" not found so omitting reference'.format(reftag), UserWarning)
                        vidx += 1

        return self._query_output

    def get_dict(self):
        """
        Return the output dictionary generated from the query
        """

        return self._query_output

    def num_pulsars(self):
        """
        Return the number of pulsars found in with query
        """

        return self._npulsars

    def table(self):
        """
        Return an astropy table of the pulsar data
        """

        from astropy.table import Table

        # make a table from the dictionary
        psrtable = Table(data=self.get_dict())

        # add units to columns
        for p in self._query_params:
            if PSR_ALL[p]['units']:
                psrtable.columns[p].unit = PSR_ALL[p]['units']

                if PSR_ALL[p]['err'] and self._include_errs:
                    psrtable.columns[p+'_ERR'].unit = PSR_ALL[p]['units']

        return psrtable

    def get_version(self):
        """
        Return a string with the ATNF version number, or the default giving in ATNF_VERSION if not found
        """

        if self._atnf_version is None:
            self._atnf_version = get_version()
        
        return self._atnf_version
