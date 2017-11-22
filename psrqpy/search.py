"""
Search query
"""

from __future__ import print_function, division

import warnings
from collections import OrderedDict
import re
import six

import numpy as np
import requests
from bs4 import BeautifulSoup

from .config import *
from .utils import *

class QueryATNF(object):
    """
    Class to generate a query of the ATNF catalogue
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None, bincomp=None,
                 exactmatch=False, sort_attr='jname', sort_order='asc', psrs=None,
                 include_errs=True, include_refs=False, get_ephemeris=False, version=None,
                 adsref=False):
        """
        Set up and perform the query of the ATNF catalogue

        :param params: a list of strings with the pulsar parameters to return
        :param condition: a string with conditions for the returned parameters
        :param psrtype: a list of strings, or single string, of conditions on the 'type' of pulsars to return (logical AND will be used for any listed types)
        :param assoc: a condition on the associations of pulsars to return (logical AND will be used for any listed associations)
        :parsm bincomp: a list of strings, or single string, of conditions on the binary companiion types of pulsars to return (logical AND will be used for any listed associations)
        :param extractmatch: a boolean stating whether assciations and types given as the condition should be an exact match
        :param sort_attr: the parameter on which with sort the returned pulsars
        :param sort_ord: the order of the sorting, either 'asc' or 'desc' (defaults to ascending)
        :param psrs: a list of pulsar names to get the information for
        :param include_errs: boolean to set whether to include parameter errors
        :param include_refs: boolean to set whether to include parameter references
        :param get_ephemeris: boolean to set whether to get pulsar ephemerides (only works if `psrs` have been specified)
        :param version: a string with the ATNF version to use (this will default to the current version if set as None)
        :param adsref: boolean to set whether the python 'ads' module can be used to get reference information
        """

        self._psrs = psrs
        self._include_errs = include_errs
        self._include_refs = include_refs
        self._atnf_version = version
        self._atnf_version = self.get_version # if no version is set this will return the current or default value
        self._adsref = adsref

        # check sort order is either 'asc' or 'desc' (or some synonyms)
        if sort_order.lower() in ['asc', 'up', '^']:
            self._sort_order = 'asc'
        elif sort_order.lower() in ['desc', 'descending', 'v']:
            self._sort_order = 'desc'
        else:
            warnings.warn('Unrecognised sort order "{}", defaulting to "ascending"'.format(sort_order), UserWarning)
            self._sort_order = 'asc'

        self._sort_attr = sort_attr

        self._refs = None # set of pulsar references
        self._query_output = None
        self._get_ephemeris = get_ephemeris

        self._pulsars = None # gets set to a Pulsars object by get_pulsars()

        # check parameters are allowed values
        if isinstance(params, list):
            if len(params) == 0:
                raise Exception("No parameters in list")

            for p in params:
                if not isinstance(p, basestring):
                    raise Exception("Non-string value '{}' found in params list".format(p))

            self._query_params = [p.upper() for p in params] # make sure parameter names are all upper case
        else:
            if isinstance(params, basestring):
                self._query_params = [params.upper()] # make sure parameter is all upper case
            else:
                if self._psrs and self._get_ephemeris: # if getting ephemerides then param can be None 
                    self._query_params = []
                else:
                    raise Exception("'params' must be a list or string")

        for p in list(self._query_params):
            if p not in PSR_ALL_PARS:
                warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                self._query_params.remove(p)
        if len(self._query_params) == 0 and (not self._psrs or not self._get_ephemeris):
            raise Exception("No parameters left in list")

        # set conditions
        self._conditions_query = self.parse_conditions(condition, psrtype=psrtype, assoc=assoc, bincomp=bincomp, exactmatch=exactmatch)

        # get references is required
        if self._include_refs:
            self._refs = get_references(useads=self._adsref)

        # perform query
        self._query_content = self.generate_query()

        # parse the query with BeautifulSoup into a dictionary
        self._query_output = self.parse_query()

    def generate_query(self, version='', params=None, condition='', sortorder='asc', sortattr='JName', psrnames=None, **kwargs):
        """
        Generate a query URL and return the content of the request from that URL. If set the class attributes are
        used for generating the query, otherwise arguments can be given.

        :param version: a string containing the ATNF version
        :param params: a list of parameters to query
        :param condition: the condition string for the query
        :param sortorder: the order for sorting the results
        :param sortattr: the attribute on which to perform the sorting
        :param psrnames: a list of pulsar names to get
        :param get_ephemeris: boolean stating whether to get pulsar ephemerides rather than a table of parameter values (only works if pulsar names are given)
        """

        # get_ephemeris is the only keyword argument at the moment
        for key, value in six.iteritems(kwargs):
            if key == 'get_ephemeris':
                if isinstance(value, bool):
                    self._get_ephemeris = value # overwrite the pre-set class _get_ephemeris value

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
                if p.upper() not in PSR_ALL_PARS:
                    warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                    qparams.remove(p)
            self._query_params = [qp.upper() for qp in qparams] # convert parameter names to all be upper case

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
            if isinstance(self._psrs, basestring):
                self._psrs = [self._psrs] # if a string pulsar name then convert to list

            for psr in self._psrs:
                if '+' in psr: # convert '+'s in pulsar names to '%2B' for the query string
                    qpulsars += psr.replace('+', '%2B')
                else:
                    qpulsars += psr
                qpulsars += '+' # seperator between pulsars
            qpulsars = qpulsars.strip('+') # remove the trailing '+'
        query_dict['psrnames'] = qpulsars

        # get pulsar ephemeris rather than table (parsing of this is not implemented yet)
        query_dict['getephemeris'] = ''
        if self._get_ephemeris:
            if self._psrs is not None:
                query_dict['getephemeris'] = 'Get+Ephemeris'
            else:
                warnings.warn('Cannot get ephemeris if no pulsar names are provided. No ephemerides will be returned.', UserWarning)
                self._get_ephemeris = False

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
        self._bad_pulsars = [] # any requested pulsars that were not found
        for pt in pretags:
            if 'WARNING' in pt.text:
                warnings.warn('Request generated warning: "{}"'.format(pt.text), UserWarning)

                # check if warning was for a specific requested pulsar: given by warning string "WARNING: PSR XXXXXXX not in catalogue"
                if 'PSR' in pt.text:
                    pat = r'WARNING: PSR (?P<psr>\S+) not in catalogue'
                    wvalues = re.search(pat, pt.text).groupdict()

                    if 'psr' in wvalues:
                        self._bad_pulsars.append(wvalues['psr'])
                        # remove any pulsars that weren't found
                        if wvalues['psr'] in self._psrs:
                            del self._psrs[wvalues['psr']]

        # actual table or ephemeris values should be in the final <pre> tag
        qoutput = pretags[-1].text
        self._query_output = OrderedDict()
        self._npulsars = 0
        self._pulsars = None # reset to None in case a previous query had already been performed

        if not self._get_ephemeris: # not getting ephemeris values
            # put the data in an ordered dictionary dictionary
            if qoutput:
                plist = qoutput.strip().split('\n') # split output string

                if self._psrs:
                    if len(self._psrs) != len(plist):
                        raise Exception('Number of pulsars returned is not the same as the number requested')

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
                    # split the line on whitespace or \xa0 using re (if just using split it ignores \xa0,
                    # which may be present for, e.g., empty reference fields, and results in the wrong
                    # number of line entries, also ignore the first entry as it is always in index
                    pvals = [lv.strip() for lv in re.split(r'\s+| \xa0 | \D\xa0', line)][1:] # strip removes '\xa0' now

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
        else: # getting ephemeris
            # split ephemerides for each requested pulsar (they are seperated by '@-----'...)
            if qoutput:
                psrephs = re.split(r'@-+', qoutput)

                if len(psrephs) != len(self._psrs):
                    raise Exception('Number of pulsar ephemerides returned is not the same as the number requested')

                self._npulsars = len(self._psrs)

                # query output in this case is a dictionary of ephemerides
                for psr, psreph in zip(self._psrs, psrephs):
                    self._query_output[psr] = psreph

        return self._query_output

    def get_dict(self):
        """
        Return the output dictionary generated from the query
        """

        return self._query_output

    @property
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

        # add catalogue version to metadata
        psrtable.meta['version'] = self.get_version
        psrtable.meta['ATNF Pulsar Catalogue'] = ATNF_BASE_URL

        return psrtable

    def get_pulsars(self):
        """
        Return the queried pulsars as a Pulsars object, which is a dictionary
        of Pulsar objects. If 'JNAME' or 'NAME' was not in the original query,
        it will be performed again, so that a name is present, which is required
        for a Pulsar object
        """

        if not self._pulsars:
            from .pulsar import Pulsar, Pulsars

            # check if JNAME or NAME was queried
            if 'JNAME' not in self._query_params and 'NAME' not in self._query_params:
                self._query_params.append('JNAME') # add JNAME parameter

                # re-do query
                self._query_content = self.generate_query()

                # parse the query with BeautifulSoup into a dictionary
                self._query_output = self.parse_query()
                nameattr = 'JNAME'
            elif 'JNAME' in self._query_params:
                nameattr = 'JNAME'
            else:
                nameattr = 'NAME'

            self._pulsars = Pulsars()

            # add pulsars one by one
            psrtable = self.table()
            for row in psrtable:
                attrs = {}
                for key in psrtable.colnames:
                    attrs[key] = row[key]

                P = Pulsar(attrs[nameattr], version=self.get_version, **attrs)
                self._pulsars.add_pulsar(P)

        return self._pulsars

    @property
    def get_version(self):
        """
        Return a string with the ATNF version number, or the default giving in
        ATNF_VERSION if not found
        """

        if self._atnf_version is None:
            self._atnf_version = get_version()

        return self._atnf_version

    def parse_conditions(self, condition, psrtype=None, assoc=None, bincomp=None, exactmatch=False):
        """
        Parse a string on conditions, i.e., logical statements on with to perform a search, like
        condition = 'f0 > 2.5 && assoc(GC)'

        :param condition: a string of conditional statements
        :param psrtype: a list of strings, or single string, of conditions on the 'type' of pulsars to return (logical AND will be used for any listed types)
        :param assoc: a list of strings, or single string, of conditions on the associations of pulsars to return (logical AND will be used for any listed associations)
        :parsm bincomp: a list of strings, or single string, of conditions on the binary companiion types of pulsars to return (logical AND will be used for any listed associations)
        :param extractmatch: a boolean stating whether assciations and types given as the condition should be an exact match
        """

        if not condition:
            conditionparse = ''
        else:
            if not isinstance(condition, basestring):
                warnings.warn('Condition "{}" must be a string. No condition being set'.format(condition), UserWarning)
                return ''

            # split condition on >, <, &&, ||, ==, <=, >=, !=, (, ), and whitespace
            splitvals = r'(&&)|(\|\|)|(>=)|>|(<=)|<|\(|\)|(==)|(!=)|!' # perform splitting by substitution and then splitting on whitespace
            condvals = re.sub(splitvals, ' ', condition).split()

            # check values are numbers, parameter, names, assocition names, etc
            for cv in condvals:
                if cv.upper() not in PSR_ALL_PARS + PSR_TYPES + PSR_BINARY_TYPE + PSR_ASSOC_TYPE:
                    # check if it's a number
                    try:
                        float(cv)
                    except ValueError:
                        warnings.warn('Unknown value "{}" in condition string "{}". No condition being set'.format(cv, condition), UserWarning)
                        return ''

            # remove spaces (turn into '+'), and convert values in condition
            conditionparse = condition.strip() # string preceeding and trailing whitespace
            conditionparse = re.sub(r'\s+', '+', conditionparse) # change whitespace to '+'

            # substitute && for %26%26
            conditionparse = re.sub(r'(&&)', '%26%26', conditionparse)

            # substitute || for %7C%7C
            conditionparse = re.sub(r'(\|\|)', '%7C%7C', conditionparse)

            # substitute '==' for %3D%3D
            conditionparse = re.sub(r'(==)', '%3D%3D', conditionparse)

            # substitute '!=' for %21%3D
            conditionparse = re.sub(r'(!=)', '%21%3D', conditionparse)

            # substitute '>=' for >%3D
            conditionparse = re.sub(r'(>=)', '>%3D', conditionparse)

            # substitute '<=' for <%3D
            conditionparse = re.sub(r'(<=)', '>%3D', conditionparse)

        # add on any extra given pulsar types
        if psrtype is not None:
            if isinstance(psrtype, list):
                if len(psrtype) == 0:
                    raise Exception("No pulsar types in list")

                for p in psrtype:
                    if not isinstance(p, basestring):
                        raise Exception("Non-string value '{}' found in pulsar type list".format(p))
                self._query_psr_types = psrtype
            else:
                if isinstance(psrtype, basestring):
                    self._query_psr_types = [psrtype]
                else:
                    raise Exception("'psrtype' must be a list or string")

            for p in list(self._query_psr_types):
                if p.upper() not in PSR_TYPES:
                    warnings.warn("Pulsar type '{}' is not recognised, no type will be required".format(p))
                    self._query_psr_types.remove(p)
                else:
                    if not conditionparse:
                        conditionparse = 'type({})'.format(p.upper())
                    else:
                        conditionparse += '+%26%26+type({})'.format(p.upper())

        # add on any extra given associations
        if assoc is not None:
            if isinstance(assoc, list):
                if len(assoc) == 0:
                    raise Exception("No pulsar types in list")

                for p in assoc:
                    if not isinstance(p, basestring):
                        raise Exception("Non-string value '{}' found in associations list".format(p))
                self._query_assocs = assoc
            else:
                if isinstance(assoc, basestring):
                    self._query_assocs = [assoc]
                else:
                    raise Exception("'assoc' must be a list or string")

            for p in list(self._query_assocs):
                if p.upper() not in PSR_ASSOC_TYPE:
                    warnings.warn("Pulsar association '{}' is not recognised, no type will be required".format(p))
                    self._query_assocs.remove(p)
                else:
                    if not conditionparse:
                        conditionparse = 'assoc({})'.format(p.upper())
                    else:
                        conditionparse += '+%26%26+assoc({})'.format(p.upper())

        # add on any extra given binary companion types
        if bincomp is not None:
            if isinstance(bincomp, list):
                if len(assoc) == 0:
                    raise Exception("No pulsar types in list")

                for p in bincomp:
                    if not isinstance(p, basestring):
                        raise Exception("Non-string value '{}' found in binary companions list".format(p))
                self._query_bincomps = bincomp
            else:
                if isinstance(bincomp, basestring):
                    self._query_bincomps = [bincomp]
                else:
                    raise Exception("'bincomp' must be a list or string")

            for p in list(self._query_bincomps):
                if p.upper() not in PSR_BINARY_TYPE:
                    warnings.warn("Pulsar binary companion '{}' is not recognised, no type will be required".format(p))
                    self._query_bincomps.remove(p)
                else:
                    if not conditionparse:
                        conditionparse = 'bincomp({})'.format(p.upper())
                    else:
                        conditionparse += '+%26%26+bincomp({})'.format(p.upper())

        if exactmatch and conditionparse:
            conditionparse += '&exact_match=match'

        return conditionparse

    def __len__(self):
        """
        Length method returns the number of pulsars
        """

        return self._npulsars

    def __str__(self):
        """
        String method returns string method of astropy table
        """
        
        if self._npulsars > 0:
            return str(self.table())
        else:
            return str(self._query_output) # should be empty dict

    def __repr__(self):
        """
        repr method returns repr method of astropy table
        """
        
        if self._npulsars > 0:
            return repr(self.table())
        else:
            return repr(self._query_output) # should be empty dict
