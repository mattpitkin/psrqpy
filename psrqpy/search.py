"""
The classes defined here are for querying the `ATNF pulsar catalogue
<http://www.atnf.csiro.au/people/pulsar/psrcat/>`_ and viewing the resulting
information.
"""

from __future__ import print_function, division

import warnings
from collections import OrderedDict
import re
import cPickle as pickle
import six

import numpy as np
import requests
from bs4 import BeautifulSoup

from .config import *
from .utils import *

class QueryATNF(object):
    """
    A class to generate a query of the
    `ATNF pulsar catalogue <http://www.atnf.csiro.au/people/pulsar/psrcat/>`_.

    Args:
        params (str, :obj:`list`, required): a list of strings with the pulsar `parameters
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#par_list>`_
            to query. The parameter names are case insensitive.
        condition (str): a string with logical conditions for the returned parameters. The allowed
            format of the condition string is given `here
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_.
            Defaults to None.
        psrtype (:obj:`list`): a list of strings, or single string, of conditions on the `type
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_ of
            pulsars to return (logical AND will be used for any listed types). Defaults to None.
        assoc (:obj:`list`, str): a condition on the associations of pulsars to return (logical AND
            will be used for any listed associations). Currently this can contain either ``GC`` for
            pulsars in globular clusters or ``SNR`` for pulsars with associated supernova remnants.
            Defaults to None.
        bincomp (str, :obj:`list`): a list of strings, or single string, of conditions on the
            `binary
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#bincomp_type>`_
            companion types of pulsars to return (logical AND will be used for any listed
            associations). Defaults to None.
        extractmatch (bool): a boolean stating whether assciations and types given as the condition
            should be an exact match. Defaults to False.
        sort_attr (str): the (case insensitive) parameter name on which with sort the returned
            pulsars. Defaults to ``JName``.
        sort_ord (str): the order of the sorting, can be either ``asc`` or ``desc``. Defaults to
            ascending.
        psrs (:obj:`list`): a list of pulsar names for which to get the requested parameters.
            Defaults to None.
        include_errs (bool): Set if wanting parameter errors to be returned. Defaults to True.
        include_refs (bool): Set if wanting parameter
            `references <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_ref.html>`_ to be
            returned. Defaults to False.
        get_ephemeris (bool): Set if wanting to get pulsar ephemerides (only works if `psrs` have
            been specified). Defaults to False.
        version (str): a string with the ATNF version to use (this will default to the current
            version if set as None)
        adsref (bool): Set if wanting to use an :class:`ads.search.SearchQuery` to get reference
            information. Defaults to False.
        loadfromfile (str): load an instance of :class:`~psrqpy.search.QueryATNF` from the given
            file, rather than performing a new query. Defaults to None.
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None, bincomp=None,
                 exactmatch=False, sort_attr='jname', sort_order='asc', psrs=None,
                 include_errs=True, include_refs=False, get_ephemeris=False, version=None,
                 adsref=False, loadfromfile=None):
        self._psrs = psrs
        self._include_errs = include_errs
        self._include_refs = include_refs
        self._atnf_version = version
        self._atnf_version = self.get_version # if no version is set this will return the current or default value
        self._adsref = adsref

        self._savefile = None # file to save class to
        self._loadfile = None # file class loaded from

        if loadfromfile:
            self.load(loadfromfile)
            return

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
            self._refs = get_references()

        # perform query
        self._query_content = self.generate_query()

        # parse the query with BeautifulSoup into a dictionary
        self._query_output = self.parse_query()

    def save(self, fname):
        """
        Output the :class:`~psrqpy.search.QueryATNF` instance to a pickle file for future loading.

        Args:
            fname (str): the filename to output the pickle to
        """

        try:
            fp = open(fname, 'wb')
            self._savefile = fname
            pickle.dump(self.__dict__, fp, 2)
            fp.close()
            self._savefile = fname
        except IOError:
            raise Exception("Error outputing class to pickle file")
          
    def load(self, fname):
        """
        Load a previously saved pickle of this class.

        Args:
            fname (str): the filename of the pickled object
        """

        try:
            fp = open(fname, 'rb')
            tmp_dict = pickle.load(fp)
            fp.close()          
            self.__dict__.clear() # clear current self
            self.__dict__.update(tmp_dict)
            self._loadfile = fname
        except IOError:
            raise Exception("Error reading in pickle")

    def generate_query(self, version='', params=None, condition='', sortorder='asc',
                       sortattr='JName', psrnames=None, **kwargs):
        """
        Generate a query URL and return the content of the :class:`~requests.Response` from that
        URL. If the required class attributes are set then they are used for generating the query,
        otherwise arguments can be given to override those set when initialising the class.

        Args:
            version (str): a string containing the ATNF version.
            params (list, str): a list of `parameters <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#par_list>`_ to query.
            condition (str): the logical `condition
                <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_
                string for the query.
            sortorder (str): the order for sorting the results.
            sortattr (str): the parameter on which to perform the sorting.
            psrnames (list, str): a list of pulsar names to get parameters for
            get_ephemeris (bool): a boolean stating whether to get pulsar ephemerides rather than
                a table of parameter values (only works if pulsar names are given)

        Returns:
            str: the :attr:`~requests.Response.content` of a :class:`~requests.Response` returned
            by :func:`~requests.get` for the ATNF pulsar catalogue query page.

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

        Args:
            requestcontent (str): The content of a :class:`~requests.Response` returned by
                :func:`~requests.get`

        Returns:
            :class:`collections.OrderedDict`: an ordered dictionary of requested parameter values

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

                                    if self._adsref:
                                        if 'ADS URL' not in thisref: # get ADS reference
                                            try:
                                                import ads
                                            except ImportError:
                                                warnings.warn('Could not import ADS module, so no ADS information will be included', UserWarning)
                                                article = []

                                            try:
                                                article = ads.SearchQuery(year=thisref['year'], first_author=thisref['authors'][0], title=thisref['title'])
                                            except IOError:
                                                warnings.warn('Could not get reference information, so no ADS information will be included', UserWarning)
                                                article = []

                                            article = list(article)
                                            if len(article) > 0:
                                                self._refs[reftag]['ADS URL'] = ADS_URL.format(list(article)[0].bibcode)

                                        self._query_output[p+'_REFURL'][idx] = thisref['ADS URL']
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
        Returns:
            :class:`~collections.OrderedDict`: the output dictionary generated by the query.
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
        Returns:
             :class:`astropy.table.Table`: a table of the pulsar data returned by the query.
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
        Returns:
            :class:`psrqpy.pulsar.Pulsars`: the queried pulsars returned as a
            :class:`~psrqpy.pulsar.Pulsars` object, which is a dictionary of
            :class:`~psrqpy.pulsar.Pulsar` objects. If ``JNAME`` or ``NAME`` was not in the
            original query, it will be performed again, so that a name is present, which is
            required for a :class:`~psrqpy.pulsar.Pulsar` object
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
        Returns:
            str: seturn a string with the ATNF version number, or the default giving in
                :attr:`~psrqpy.config.ATNF_VERSION` if not found.
        """

        if self._atnf_version is None:
            self._atnf_version = get_version()

        return self._atnf_version

    def parse_conditions(self, condition, psrtype=None, assoc=None, bincomp=None, exactmatch=False):
        """
        Parse a string of `conditions
        <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_, i.e.,
        logical statements with which to apply to a catalogue query, e.g.,
        ``condition = 'f0 > 2.5 && assoc(GC)'``, so that they are in the format required for the
        query URL.

        Args:
            condition (str): a string of `conditional <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_
                statements
            psrtype (list, str): a list of strings, or single string, of conditions on the
                `type <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_ of
                pulsars to return (logical AND will be used for any listed types)
            assoc (list, str): a list of strings, or single string, of conditions on the
                associations of pulsars to return (logical AND will be used for any listed
                associations)
            bincomp (list, str): a list of strings, or single string, of conditions on the
                `binary companion <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type>`_
                types of pulsars to return (logical AND will be used for any listed associations)
            exactmatch (bool): a boolean stating whether assciations and types given as the
                condition should be an exact match

        Returns:
            str: a string with the format required for use in :attr:`~psrqpy.config.QUERY_URL`

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
        Returns:
            int: :func:`len` method returns the number of pulsars
        """

        return self._npulsars

    def __str__(self):
        """
        Returns:
            str: :func:`str` method returns the str method of an :class:`astropy.table.Table`.
        """

        if self._npulsars > 0:
            return str(self.table())
        else:
            return str(self._query_output) # should be empty dict

    def __repr__(self):
        """
        Returns:
            str: :func:`repr` method returns the repr method of an :class:`astropy.table.Table`.
        """

        if self._npulsars > 0:
            return repr(self.table())
        else:
            return repr(self._query_output) # should be empty dict

    def ppdot(self, intrinsicpdot=False, excludeGCs=False, showtypes=[], showGCs=False,
              showSNRs=False, markertypes={}, deathline=True, deathmodel='Ip', filldeath=True,
              filldeathtype={}, showtau=True, brakingidx=3, tau=None, showB=True, Bfield=None,
              rcparams={}):
        """
        Draw a lovely period vs period derivative diagram.

        Args:
            intrinsicpdot (bool): use the intrinsic period derivative corrected for the
                `Shklovskii effect <https://en.wikibooks.org/wiki/Pulsars_and_neutron_stars/Pulsar_properties#Pulse_period>`_
                rather than the observed value. Defaults to False.
            excludeGCs (bool): exclude globular cluster pulsars as their period derivatives can be
                contaminated by intra-cluster accelerations. Defaults to False.
            showtypes (list, str): a list of pulsar types to highlight with markers in the plot.
                These can contain any of the following: ``BINARY``, ``HE``, ``NRAD``, ``RRAT``,
                ``XINS``, ``AXP`` or ``SGR``, or ``ALL`` to show all types. Default to showing no
                types.
            showGCs (bool): show markers to denote the pulsars in globular clusters. Defaults to
                False.
            showSNRs (bool): show markers to denote the pulsars with supernova remnants associated
                with them. Defaults to False.
            markertypes (dict): a dictionary of marker styles and colors keyed to the pulsar types
                above
            deathline (bool): draw the pulsar death line. Defaults to True.
            deathmodel (str): the type of death line to draw based on the models in
                :func:`psrqpy.utils.death_line`. Defaults to ``'Ip'``.
            filldeath (bool): set whether to fill the pulsar graveyard under the death line.
                Defaults to True.
            filldeathtype (dict): a dictionary of keyword arguments for the fill style of the
                pulsar graveyard.
            showtau (bool): show lines for a selection of characteritic ages. Defaults to True,
                and shows lines for :math:`10^5` through to :math:`10^9` yrs with steps in powers
                of 10.
            brakingidx (int): a braking index to use for the calculation of the characteristic age
                lines. Defaults to 3 for magnetic dipole radiation.
            tau (list): a list of characteristic ages to show on the plot.
            showB (bool): show lines of constant magnetic field strength. Defaults to True, and
                shows lines for :math:`10^{10}` through to :math:`10^{14}` gauss with steps in
                powers of 10.
            Bfield (list): a list of magnetic field strengths to plot.
            rcparams (dict): a dictionary of :py:obj:`matplotlib.rcParams` setup parameters for the
                plot.

        Returns:
            :class:`matplotlib.figure.Figure`: the figure object
        """

        try:
            import matplotlib as mpl
            from matplotlib import pyplot as pl
        except ImportError:
            raise Exception('Cannot produce P-Pdot plot as Matplotlib is not available')

        # check the we have periods and period derivatives
        nparams = len(self._query_params)
        if 'P0' not in self._query_params:
            self._query_params.append('P0')
        if 'P1' not in self._query_params:
            self._query_params.append('P1')

        # check if we want to use intrinsic period derivatives
        if intrinsicpdot and 'P1_I' not in self._query_params:
            self._query_params.append('P1_I')

        if isinstance(showtypes, basestring):
            nshowtypes = [showtypes]
        else:
            nshowtypes = showtypes

        for stype in nshowtypes:
            if 'ALL' == stype.upper():
                nshowtypes = PSR_TYPES
                del nshowtypes[nshowtypes.index('RADIO')] # remove radio as none are returned as this 
                break
            if 'SGR' == stype.upper(): # synonym for AXP
                nshowtypes[nshowtypes.index(stype)] = 'AXP'

        if nshowtypes and 'TYPE' not in self._query_params:
            self._query_params.append('TYPE')

        if 'BINARY' in nshowtypes and 'BINARY' not in self._query_params:
            self._query_params.append('BINARY')

        if showGCs or showSNRs:
            if 'ASSOC' not in self._query_params:
                self._query_params.append('ASSOC')

        # redo query if required
        if len(self._query_params) != nparams:
            # perform query
            self._query_content = self.generate_query()
            self._query_output = self.parse_query()

        if not self.num_pulsars:
            print("No pulsars found, so no P-Pdot plot has been produced")
            return None

        # set plot parameters
        rcparams['figure.figsize'] = rcparams['figure.figsize'] if 'figure.figsize' in rcparams else (9, 9.5)
        rcparams['figure.dpi'] = rcparams['figure.dpi'] if 'figure.dpi' in rcparams else 250
        rcparams['text.usetex'] = rcparams['text.usetex'] if 'text.usetex' in rcparams else True
        rcparams['axes.linewidth'] = rcparams['axes.linewidth'] if 'axes.linewidth' in rcparams else 0.5
        rcparams['axes.grid'] = rcparams['axes.grid'] if 'axes.grid' in rcparams else False
        rcparams['font.family'] = rcparams['font.family'] if 'font.family' in rcparams else 'sans-serif'
        rcparams['font.sans-serif'] = rcparams['font.sans-serif'] if 'font.sans-serif' in rcparams else 'Avant Garde, Helvetica, Computer Modern Sans serif'
        rcparams['font.size'] = rcparams['font.size'] if 'font.size' in rcparams else 20
        rcparams['legend.fontsize'] = rcparams['legend.fontsize'] if 'legend.fontsize' in rcparams else 16
        rcparams['legend.frameon'] = rcparams['legend.frameon'] if 'legend.frameon' in rcparams else False

        mpl.rcParams.update(rcparams)

        fig, ax = pl.subplots()

        # get astropy table of parameters
        t = self.table()

        # extract periods and period derivatives
        periods = t['P0']
        pdots = t['P1']
        if intrinsicpdot: # use instrinsic period derivatives if requested
            ipdotidx = np.isfinite(t['P1_I'])
            pdots[ipdotidx] = t['P1_I'][ipdotidx]

        # get only finite values
        pidx = (np.isfinite(periods)) & (np.isfinite(pdots))
        periods = periods[pidx]
        pdots = pdots[pidx]

        if 'ASSOC' in self._query_params:
            assocs = t['ASSOC'][pidx]   # associations
        if 'TYPE' in self._query_params:
            types = t['TYPE'][pidx]     # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = t['BINARY'][pidx] # binary pulsars

        # now get only positive pdot values
        pidx = pdots > 0.
        periods = periods[pidx]
        pdots = pdots[pidx]
        if 'ASSOC' in self._query_params:
            assocs = assocs[pidx]   # associations
        if 'TYPE' in self._query_params:
            types = types[pidx]     # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = binaries[pidx] # binary pulsars

        # check whether to exclude globular cluster pulsars that could have contaminated spin-down value
        if excludeGCs:
            nongcidxs = np.flatnonzero(np.char.find(assocs,'GC:')==-1) # use '!=' to find GC indexes
            periods = periods[nongcidxs]
            pdots = pdots[nongcidxs]
            if 'ASSOC' in self._query_params:
                assocs = assocs[nongcidxs]
            if 'TYPE' in self._query_params:
                types = types[nongcidxs]
            if 'BINARY' in nshowtypes:
                binaries = binaries[nongcidxs]

        # plot pulsars
        ax.loglog(periods, pdots, marker='.', color='dimgrey', linestyle='none')
        ax.set_xlabel(r'Period (s)')
        ax.set_ylabel(r'Period Derivative')

        # get limits
        periodlims = [10**np.floor(np.min(np.log10(periods))), 10.*int(np.ceil(np.max(pdots)/10.))]
        pdotlims = [10**np.floor(np.min(np.log10(pdots))), 10**np.ceil(np.max(np.log10(pdots)))]
        ax.set_xlim(periodlims);
        ax.set_ylim(pdotlims);

        if deathline:
            deathpdots = 10**death_line(np.log10(periodlims), linemodel=deathmodel)
            ax.loglog(periodlims, deathpdots, 'k--', linewidth=0.5)

            if filldeath:
                if not filldeathtype:
                    filldeathtype = {}

                filldeathtype['linestyle'] = filldeathtype['linestyle'] if 'linestyle' in filldeathtype else '-'
                filldeathtype['alpha'] = filldeathtype['alpha'] if 'alpha' in filldeathtype else 0.15
                filldeathtype['facecolor'] = filldeathtype['facecolor'] if 'facecolor' in filldeathtype else 'darkorange'
                filldeathtype['hatch'] = filldeathtype['hatch'] if 'hatch' in filldeathtype else ''
                ax.fill_between(periodlims, deathpdots, pdotlims[0], **filldeathtype)

        # add markers for each pulsar type
        if not markertypes:
            markertypes = {}
    
        # check if markers have been defined by the user or not
        markertypes['AXP'] = {'marker': 's', 'markeredgecolor': 'red'} if 'AXP' not in markertypes else markertypes['AXP']
        markertypes['BINARY'] = {'marker': 'o', 'markeredgecolor': 'grey'} if 'BINARY' not in markertypes else markertypes['BINARY']
        markertypes['HE'] = {'marker': 'D', 'markeredgecolor': 'orange'} if 'HE' not in markertypes else markertypes['HE']
        markertypes['RRAT'] = {'marker': 'h', 'markeredgecolor': 'green'} if 'RRAT' not in markertypes else markertypes['RRAT']
        markertypes['NRAD'] = {'marker': 'v', 'markeredgecolor': 'blue'} if 'NRAD' not in markertypes else markertypes['NRAD']
        markertypes['XINS'] = {'marker': '^', 'markeredgecolor': 'magenta'} if 'XINS' not in markertypes else markertypes['XINS']
        markertypes['GC'] = {'marker': '8', 'markeredgecolor': 'cyan'} if 'GC' not in markertypes else markertypes['GC']
        markertypes['SNR'] = {'marker': '*', 'markeredgecolor': 'darkorchid'} if 'SNR' not in markertypes else markertypes['SNR']

        # legend strings for different types
        typelegstring = {}
        typelegstring['AXP'] = r'SGR/AXP'
        typelegstring['NRAD'] = r'"Radio-Quiet"'
        typelegstring['XINS'] = r'Pulsed Thermal X-ray'
        typelegstring['BINARY'] = r'Binary'
        typelegstring['HE'] = r'Radio-IR Emission'
        typelegstring['GC'] = r'Globular Cluster'

        # show globular cluster pulsars
        if showGCs and not excludeGCs:
            nshowtypes.append('GC')

        # show pulsars with associated supernova remnants
        if showSNRs:
            nshowtypes.append('SNR')

        handles = OrderedDict()

        for stype in nshowtypes:
            if stype.upper() in PSR_TYPES + ['GC', 'SNR']:
                thistype = stype.upper()
                if thistype == 'BINARY':
                    # for binaries used the 'BINARY' column in the table
                    typeidx = np.flatnonzero(np.char.find(binaries, '*')==-1)
                elif thistype in ['GC', 'SNR']:
                    typeidx = np.flatnonzero(np.char.find(assocs, thistype)!=-1)
                else:
                    typeidx = np.flatnonzero(np.char.find(types, thistype)!=-1)

                if len(typeidx) == 0:
                    continue

                # default to empty markers with no lines between them
                if 'markerfacecolor' not in markertypes[thistype]:
                    markertypes[thistype]['markerfacecolor'] = 'none'
                if 'linestyle' not in markertypes[thistype]:
                    markertypes[thistype]['linestyle'] = 'none'
                typehandle, = ax.loglog(periods[typeidx], pdots[typeidx], **markertypes[thistype])
                if thistype in typelegstring:
                    handles[typelegstring[thistype]] = typehandle
                else:
                    handles[thistype] = typehandle

                ax.legend(handles.values(), handles.keys(), loc='upper left', numpoints=1);

        # add characteristic age lines
        tlines = OrderedDict()
        if showtau:
            if tau is None:
                taus = [1e5, 1e6, 1e7, 1e8, 1e9] # default characteristic ages
            else:
                taus = tau

            nbrake = brakingidx
            for tauv in taus:
                pdots_tc = age_pdot(periodlims, tau=tauv, braking_idx=nbrake)
                tline, = ax.loglog(periodlims, pdots_tc, 'k-.', linewidth=0.5)
                # check if taus are powers of 10
                taupow = np.floor(np.log10(tauv))
                numv = tauv/10**taupow
                if numv == 1.:
                    tlines[r'$10^{{{0:d}}}\,{{\rm yr}}$'.format(int(taupow))] = tline
                else:
                    tlines[r'${{0:.1f}}!\times\!10^{{{1:d}}}\,{{\rm yr}}$'.format(numv, taupow)] = tline

        # add magnetic field lines
        Blines = OrderedDict()
        if showB:
            if Bfield is None:
                Bs = [1e10, 1e11, 1e12, 1e13, 1e14]
            else:
                Bs = Bfield
            for B in Bs:
                pdots_B = B_field_pdot(periodlims, Bfield=B)
                bline, = ax.loglog(periodlims, pdots_B, 'k:', linewidth=0.5)
                # check if Bs are powers of 10
                Bpow = np.floor(np.log10(B))
                numv = B/10**Bpow
                if numv == 1.:
                    Blines[r'$10^{{{0:d}}}\,{{\rm G}}$'.format(int(Bpow))] = bline
                else:
                    Blines[r'${{0:.1f}}!\times\!10^{{{1:d}}}\,{{\rm G}}$'.format(numv, Bpow)] = bline

        fig.tight_layout()

        # add text for characteristic age lines and magnetic field strength lines
        for l in tlines:
            ttext = label_line(ax, tlines[l], l, color='k', fs=18, frachoffset=0.05)

        for l in Blines:
            ttext = label_line(ax, Blines[l], l, color='k', fs=18, frachoffset=0.90)

        # return the figure
        return fig
