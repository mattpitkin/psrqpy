"""
The classes defined here are for querying the `ATNF pulsar catalogue
<http://www.atnf.csiro.au/people/pulsar/psrcat/>`_ and viewing the resulting
information.
"""

from __future__ import print_function, division

import warnings
from collections import OrderedDict
import re
import six

from six.moves import cPickle as pickle
from six import string_types

import numpy as np
import requests
from bs4 import BeautifulSoup

from astropy.coordinates import SkyCoord
import astropy.units as aunits
from astropy.table import Table

from .config import *
from .utils import *


# set formatting of warnings to not include line number and code (see
# e.g. https://pymotw.com/3/warnings/#formatting)
def warning_format(message, category, filename, lineno, file=None, line=None):
    return '{}: {}'.format(category.__name__, message)


warnings.formatwarning = warning_format


class QueryATNF(object):
    """
    A class to generate a query of the
    `ATNF pulsar catalogue <http://www.atnf.csiro.au/people/pulsar/psrcat/>`_.
    By default this class will download and cache the latest version of the
    catalogue database file, although a query can be generated from the
    catalogue webform interface if requested. The catalogue can be queried for
    specific pulsar parameters and for specific named pulsars. Conditions on
    the parameter can be specified. The results will be stored as an
    :class:`astropy.table.Table`.

    Args:
        params (str, :obj:`list`): a list of strings with the
            pulsar `parameters
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=expert#par_list>`_
            to query. The parameter names are case insensitive. If this is not
            given then all parameters will be returned by default, unless
            querying via the webform in which case only `JNAME` will be
            returned by default.
        condition (str): a string with logical conditions for the returned
            parameters. The allowed format of the condition string is given
            `here
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_.
            Defaults to None.
        psrtype (:obj:`list`): a list of strings, or single string, of
            conditions on the `type
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_
            of pulsars to return (logical AND will be used for any listed
            types). Defaults to None.
        assoc (:obj:`list`, str): a condition on the associations of pulsars to
            return (logical AND will be used for any listed associations).
            Currently this can contain either ``GC`` for pulsars in globular
            clusters or ``SNR`` for pulsars with associated supernova remnants.
            Defaults to None.
        bincomp (str, :obj:`list`): a list of strings, or single string, of
            conditions on the
            `binary
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#bincomp_type>`_
            companion types of pulsars to return (logical AND will be used for
            any listed associations). Defaults to None.
        exactmatch (bool): a boolean stating whether associations and types
            given as the condition should be an exact match. Defaults to False.
        sort_attr (str): the (case insensitive) parameter name on which with
            sort the returned pulsars. Defaults to ``JName``.
        sort_ord (str): the order of the sorting, can be either ``asc`` or
            ``desc``. Defaults to ascending.
        psrs (:obj:`list`): a list of pulsar names for which to get the
            requested parameters. Defaults to None.
        circular_boundary (:obj:`list`, tuple): a list containing three entries
            defining the centre (in right ascension and declination), and
            radius of a circle in which to search for and return pulsars. The
            first entry is the centre point right ascension as a string in
            format 'hh:mm:ss' or a float in radians, the second entry is the
            centre point declination as a string in format 'dd:mm:ss' or a
            float in radians, the final entry is the circle's radius in
            degrees. Alternatively `coord1`, `coord2`, and `radius` can be
            used.
        coord1 (str): a string containing a right ascension in the format
            ('hh:mm:ss') that centres a circular boundary in which to search
            for pulsars (requires coord2 and radius to be set).
        coord2 (str): a string containing a declination in the format
            ('dd:mm:ss') that centres a circular boundary in which to search
            for pulsars (requires coord1 and radius to be set).
        radius (float): the radius (in degrees) of a circular boundary in which
            to search for pulsars (requires coord1 and coord2 to be set).
        include_errs (bool): Set if wanting parameter errors to be returned.
            Defaults to True.
        include_refs (bool): Set if wanting parameter
            `references <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_ref.html>`_
            to be returned. Defaults to False.
        get_ephemeris (bool): Set if wanting to get pulsar ephemerides (only
            works if `psrs` have been specified). Defaults to False.
        adsref (bool): Set if wanting to use an :class:`ads.search.SearchQuery`
            to get reference information. Defaults to False.
        loadfromdb (str): Load a pulsar database file from a given path rather
            than using the ATNF Pulsar Catalogue database. Defaults to None.
        loadquery (str): load an instance of :class:`~psrqpy.search.QueryATNF`
            from the given file, rather than performing a new query. This was
            `loadfromfile` in earlier versions, which still works but has been
            deprecated. Defaults to None.
        cache (bool): Cache the catalogue database file for future use. This is
            ignored if `loadfromdb` is given or the request is via the webform.
            Defaults to True.
        forceupdate (bool): Remove a cached file, so that the catalogue will be
            re-downloaded. This is ignored if `loadfromdb` is given or the
            request is via the webform. Defaults to False.
        webform (bool): Query the catalogue webform rather than downloading the
           database file. Defaults to False.
        version (str): A string with the ATNF version to use. This will only be
            used if querying via the webform and will default to the current
            version if set as None.
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None,
                 bincomp=None, exactmatch=False, sort_attr='jname',
                 sort_order='asc', psrs=None, include_errs=True,
                 include_refs=False, get_ephemeris=False, version=None,
                 adsref=False, loadfromfile=None, loadquery=None,
                 loadfromdb=None, cache=True, forceupdate=False,
                 circular_boundary=None, coord1=None, coord2=None, radius=0.,
                 webform=False):
        if loadfromfile is not None and loadquery is None:
            loadquery = loadfromfile
        if loadquery:
            self.load(loadquery)
            return
        
        self._psrs = psrs
        self._include_errs = include_errs
        self._include_refs = include_refs
        self._atnf_version = version
        self._adsref = adsref
        self._savefile = None  # file to save class to
        self._loadfile = None  # file class loaded from
        self.__table = Table()
        self._condition = condition
        self._exactmatch = exactmatch
        self._sort_order = sort_order
        self._sort_attr = sort_attr.upper()
        self._dbfile = loadfromdb
        self._webform = webform

        if not self._webform:
            # download and cache (if requested) the database file
            try:
                self.__table = get_catalogue(path_to_db=self._dbfile,
                                             cache=cache,
                                             update=forceupdate)
            except IOError:
                raise IOError("Could not get catalogue database file")
            self._atnf_version = self._table.meta['version']
        else:
            # if no version is set this will return the current or default value
            self._atnf_version = self.get_version

        self._refs = None  # set of pulsar references
        self._get_ephemeris = get_ephemeris

        self._pulsars = None  # gets set to a Pulsars object by get_pulsars()

        # conditions for finding pulsars within a circular boundary
        self._coord1 = coord1
        self._coord2 = coord2
        self._radius = radius
        if (isinstance(circular_boundary, list) or
                isinstance(circular_boundary, tuple)):
            if len(circular_boundary) != 3:
                raise Exception("Circular boundary must contain three values")
            self._coord1 = circular_boundary[0]
            self._coord2 = circular_boundary[1]
            self._radius = circular_boundary[2]
        elif self._coord1 is None or self._coord2 is None or self._radius == 0.:
            # if any are not set then we can't define a boundary
            self._coord1 = self._coord2 = ''
            self._radius = 0.
        else:
            if (not isinstance(self._coord1, string_types)
                    or not isinstance(self._coord2, string_types)):
                raise Exception("Circular boundary centre coordinates must "
                                "be strings")
            if not isinstance(self._radius, float) and not isinstance(self._radius, int):
                raise Exception("Circular boundary radius must be a float or "
                                "int")

        # check parameters are allowed values
        if webform:
            self._query_params = ['JNAME']  # query JNAME by default for all queries
        else:
            self._query_params = None
        if isinstance(params, list):
            if len(params) == 0:
                print('No query parameters have been specified, so only '
                      '"JNAME" will be queried')

            for p in params:
                if not isinstance(p, string_types):
                    raise Exception("Non-string value '{}' found in params "
                                    "list".format(p))

            self._query_params += [p.upper() for p in params]
        else:
            if isinstance(params, string_types):
                self._query_params += [params.upper()]  # make sure parameter is all upper case
            elif params is not None:
                if self._psrs and self._get_ephemeris:  # if getting ephemerides then param can be None
                    self._query_params = None
                else:
                    raise Exception("'params' must be a list or string")

        self._coord = None
        if self._coord1 and self._coord2 and self._radius != 0.:
            # set centre coordinate as an astropy SkyCoord
            coord = SkyCoord(self._coord1, self._coord2,
                             unit=(aunits.hourangle, aunits.deg))

            # make sure 'RAJ' and 'DECJ' are queried
            self._query_params.append('RAJ')
            self._query_params.append('DECJ')

        # remove any duplicate
        if self._query_params is not None:
            self._query_params = list(set(self._query_params))

            for p in list(self._query_params):
                if p not in PSR_ALL_PARS:
                    warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                    self._query_params.remove(p)
            if len(self._query_params) == 0 and (not self._psrs or not self._get_ephemeris):
                raise Exception("No parameters left in list")

        # set conditions
        condparse = self.parse_conditions(psrtype=psrtype, assoc=assoc,
                                          bincomp=bincomp)
        if len(condparse) > 0:
            if self._condition is None:
                self._condition = condparse
            else:
                self._condition += condparse

        # get references if required
        if self._include_refs:
            self._refs = get_references()  # CHANGE GET_REFERENCES TO USE FILE IN TARBALL

        if self._webform:
            # perform query
            self.generate_query()

            # parse the query with BeautifulSoup into a dictionary
            self.parse_query()

        # perform sorting
        self.sort()

    def sort(self, sort_attr=None, sort_order=None):
        """
        Sort the generated catalogue :class:`~astropy.table.Table` on a given
        attribute and in either ascending or descending order.
        """

        if sort_attr is None:
            if self._sort_attr is None:
                self._sort_attr = 'JNAME'  # sort by name by default
        else:
            self._sort_attr = sort_attr.upper()

        if self._sort_attr not in self._table.colnames:
            raise KeyError("Sorting by attribute '{}' is not possible as it "
                           "is not in the table".format(self._sort_attr))

        if sort_order is None:
            if self._sort_order is None:
                self._sort_order = 'asc'  # sort ascending by default
        else:
            self._sort_order = sort_order

        # check sort order is either 'asc' or 'desc' (or some synonyms)
        if self._sort_order.lower() in ['asc', 'ascending', 'up', '^']:
            self._sort_order = 'asc'
        elif self._sort_order.lower() in ['desc', 'descending', 'down', 'v']:
            self._sort_order = 'desc'
        else:
            warnings.warn(('Unrecognised sort order "{}", defaulting to'
                           '"ascending"').format(sort_order), UserWarning)
            self._sort_order = 'asc'

        # sort the table
        self.__table.sort(self._sort_attr)

        # reverse table if in descending order
        if self._sort_order == 'desc':
            self.__table.reverse()

    def save(self, fname):
        """
        Output the :class:`~psrqpy.search.QueryATNF` instance to a pickle file
        for future loading.

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
            self.__dict__.clear()  # clear current self
            self.__dict__.update(tmp_dict)
            self._loadfile = fname
        except IOError:
            raise Exception("Error reading in pickle")

    def generate_query(self, version='', params=None, condition='',
                       psrnames=None, coord1='', coord2='', radius=0.,
                       **kwargs):
        """
        Generate a query URL and return the content of the
        :class:`~requests.Response` from that URL. If the required class
        attributes are set then they are used for generating the query,
        otherwise arguments can be given to override those set when
        initialising the class.

        Args:
            version (str): a string containing the ATNF version.
            params (list, str): a list of `parameters <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#par_list>`_
                to query.
            condition (str): the logical `condition
                <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_
                string for the query.
            psrnames (list, str): a list of pulsar names to get parameters for
            coord1 (str): a string containing a right ascension in the format
                ('hh:mm:ss') that centres a circular boundary in which to
                search for pulsars (requires `coord2` and `radius` to be set).
            coord2 (str): a string containing a declination in the format
                ('dd:mm:ss') that centres a circular boundary in which to
                search for pulsars (requires `coord1` and `radius` to be set).
            radius (float): the radius (in degrees) of a circular boundary in
                which to search for pulsars (requires `coord1` and `coord2` to
                be set).
            get_ephemeris (bool): a boolean stating whether to get pulsar
                ephemerides rather than a table of parameter values (only works
                if pulsar names are given)

        """

        # get_ephemeris is the only keyword argument at the moment
        for key, value in six.iteritems(kwargs):
            if key == 'get_ephemeris':
                if isinstance(value, bool):
                    self._get_ephemeris = value  # overwrite the pre-set class _get_ephemeris value

        query_dict = {}
        self._atnf_version = self._atnf_version if not version else version
        query_dict['version'] = self._atnf_version

        if params:
            if isinstance(params, string_types):
                params = [params]  # convert to list
            else:
                if not isinstance(params, list):
                    raise Exception('Error... input "params" for generate_query() must be a list')
            qparams = list(params)
            for p in params:
                if p.upper() not in PSR_ALL_PARS:
                    warnings.warn("Parameter {} not recognised".format(p), UserWarning)
                    qparams.remove(p)
            self._query_params = [qp.upper() for qp in qparams]  # convert parameter names to all be upper case

        pquery = ''
        for p in self._query_params:
            pquery += '&{}={}'.format(p, p)

        query_dict['params'] = pquery
        self._coord1 = self._coord1 if not coord1 else coord1
        self._coord2 = self._coord2 if not coord2 else coord2
        self._radius = self._radius if not radius else radius

        if psrnames:
            if isinstance(psrnames, string_types):
                self._psrs = [psrnames]  # convert to list
            else:
                if not isinstance(psrnames, list):
                    raise Exception('Error... input "psrnames" for generate_query() must be a list')
                self._psrs = list(psrnames)  # reset self._psrs

        qpulsars = ''  # pulsar name query string
        if self._psrs is not None:
            if isinstance(self._psrs, string_types):
                self._psrs = [self._psrs]  # if a string pulsar name then convert to list

            for psr in self._psrs:
                if '+' in psr:  # convert '+'s in pulsar names to '%2B' for the query string
                    qpulsars += psr.replace('+', '%2B')
                else:
                    qpulsars += psr
                qpulsars += '+'  # seperator between pulsars
            qpulsars = qpulsars.strip('+')  # remove the trailing '+'
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

        self._query_content = psrrequest.content

    def parse_query(self, requestcontent=''):
        """
        Parse the query returned by requests.

        Args:
            requestcontent (str): The content of a :class:`~requests.Response` returned by
                :func:`~requests.get`

        """

        # update request if required
        self._query_content = requestcontent if requestcontent else self._query_content

        # parse through BeautifulSoup
        try:
            psrsoup = BeautifulSoup(self._query_content, 'html.parser')
        except RuntimeError:
            raise RuntimeError('Error... problem parsing catalogue with BeautifulSoup')

        pretags = psrsoup.find_all('pre')  # get any <pre> html tags

        if pretags is None:
            # couldn't find anything, or their was a query problem
            raise Exception('Error... problem parsing catalogue for currently requested parameters')

        # check for any warnings generated by the request
        self._bad_pulsars = []  # any requested pulsars that were not found
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
                            del self._psrs[self._psrs.index(wvalues['psr'])]

                            # if there are no pulsars left in the list then return None
                            if len(self._psrs) == 0:
                                print('No requested pulsars were found in the catalogue')
                                query_output = None
                                self._npulsars = 0
                                self._pulsars = None
                                self.__table = Table()  # empty table
                                return

        # actual table or ephemeris values should be in the final <pre> tag
        qoutput = pretags[-1].text
        query_output = []  # list to contain dictionary of pulsars
        self._npulsars = 0
        self._pulsars = None  # reset to None in case a previous query had already been performed

        if not self._get_ephemeris:  # not getting ephemeris values
            # put the data in an ordered dictionary dictionary
            if qoutput:
                plist = qoutput.strip().split('\n')  # split output string

                if self._psrs:
                    if len(self._psrs) != len(plist):
                        raise Exception('Number of pulsars returned is not the same as the number requested')

                self._npulsars = len(plist)

                for p in self._query_params:
                    if p in PSR_ALL_PARS:
                        query_output[p] = []

                        if PSR_ALL[p]['err'] and self._include_errs:
                            query_output[p+'_ERR'] = []

                        if PSR_ALL[p]['ref'] and self._include_refs:
                            query_output[p+'_REF'] = []

                            if self._adsref:  # also add reference URL for NASA ADS
                                query_output[p+'_REFURL'] = []

                query_output.append({})

                for idx, line in enumerate(plist):
                    # split the line on whitespace or \xa0 using re (if just
                    # using split it ignores \xa0, which may be present for,
                    # e.g., empty reference fields, and results in the wrong
                    # number of line entries, also ignore the first entry as it
                    # is always in index
                    pvals = [lv.strip() for lv in re.split(r'\s+| \xa0 | \D\xa0', line)][1:]  # strip removes '\xa0' now

                    vidx = 0  # index of current value
                    for p in self._query_params:
                        if pvals[vidx] != '*':
                            try:
                                query_output[-1][p] = float(pvals[vidx])
                            except ValueError:
                                query_output[-1][p] = pvals[vidx]
                        vidx += 1

                        # get errors
                        if PSR_ALL[p]['err']:
                            if self._include_errs:
                                if pvals[vidx] != '*':
                                    try:
                                        query_output[-1][p+'_ERR'] = float(pvals[vidx])
                                    except ValueError:
                                        raise ValueError("Problem converting error value to float")

                            vidx += 1

                        # get references
                        if PSR_ALL[p]['ref']:
                            if self._include_refs:
                                reftag = pvals[vidx]

                                if reftag in self._refs:
                                    thisref = self._refs[reftag]
                                    refstring = ('{authorlist}, {year}, '
                                                 '{title}, {journal}, '
                                                 '{volume}')
                                    # remove any superfluous whitespace
                                    refstring2 = re.sub(r'\s+', ' ',
                                                   refstring.format(**thisref))
                                    query_output[-1][p+'_REF'] = ','.join([a for a in refstring2.split(',') if a.strip()])  # remove any superfluous empty ',' seperated values

                                    if self._adsref:
                                        if 'ADS URL' not in thisref:  # get ADS reference
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

                                        query_output[-1][p+'_REFURL'] = thisref['ADS URL']
                                else:
                                    if reftag != '*':
                                        warnings.warn('Reference tag "{}" not found so omitting reference'.format(reftag), UserWarning)
                            vidx += 1
        else:  # getting ephemeris
            # split ephemerides for each requested pulsar (they are seperated by '@-----'...)
            if qoutput:
                psrephs = re.split(r'@-+', qoutput)

                if len(psrephs) != len(self._psrs):
                    raise Exception('Number of pulsar ephemerides returned is not the same as the number requested')

                self._npulsars = len(self._psrs)

                # query output in this case is a dictionary of ephemerides
                for psr, psreph in zip(self._psrs, psrephs):
                    query_output.append({})
                    query_output[-1][psr] = psreph

        # set the table
        _ = self.table(query_output, self._query_params)

    def as_array(self):
        """
        Returns:
            :class:`~numpy.ndarray`: the output table as an array.
        """

        return self.table.as_array()

    @property
    def num_pulsars(self):
        """
        Return the number of pulsars found in with query
        """

        return self._npulsars

    @property
    def table(self):
        table = self.__table[self._query_params]

        if self._condtion is not None:
            # apply condition
            table = condition(table, self._condition, self._exactmatch) 

        if (self._coord is not None and 'RAJ' in table.colnames
                and 'DECJ' in table.colnames):
            # apply sky coordinate constraint
            catalog = SkyCoord(table['RAJ'], table['DECJ'],
                               unit=(aunits.hourangle, aunits.deg))

            # get seperations
            d2d = self._coord.separation(catalog)  

            # find seperations within required radius
            catalogmsk = d2d < self._radius*aunits.deg

            table = table[catalogmsk]

        return table

    def table(self, query_list=None, query_params=None, usecondtion=True,
              useseparation=True):
        """
        Set and return an :class:`astropy.table.Table` of the query.

        Args:
            query_list (list): a list of dictionarys of pulsar parameters
                for each pulsar as returned by a query. These are converted
                and set as the query table. If this is None and a table already
                exists then that table will be returned.
            query_params (str, list): a parameter, or list of parameters, to
                return from the query. If this is None then all parameters are
                returned.
            usecondition (bool, str): If True then the condition parsed to the
                :class:`psrqpy.QueryATNF`: class will be used when returning
                the table. If False no condition will be applied to the
                returned table. If a string is given then that will be the
                assumed condition string.
            useseparation (bool): If True and a set of sky coordinates and
                radius around which to return pulsars was set in the
                :class:`psrqpy.QueryATNF`: class then only pulsars within the
                given radius of the sky position will be returned. Otherwise
                all pulsars will be returned.

        Returns:
             :class:`astropy.table.Table`: a table of the pulsar data returned
                 by the query.
        """

        if query_list is not None:
            if not isinstance(query_list, list):
                raise TypeError("Query list is not a list!")

            from pandas import DataFrame

            # add RA and DEC in degs and JNAME/BNAME if necessary
            for i, psr in enumerate(list(query_list)):
                # add RA and DEC in degs
                if np.all([rd in psr.keys() for rd in ['RAJ', 'DECJ']]):
                    if not np.all([rd in psr.keys() for rd in ['RAJD', 'DECJD']]):
                        coord = SkyCoord(psr['RAJ'], psr['DECJ'],
                                         unit=(aunits.hourangle, aunits.deg))
                        query_list[i]['RAJD'] = coord.ra.deg    # right ascension in degrees
                        query_list[i]['DECJD'] = coord.dec.deg  # declination in degrees

                # add 'JNAME', 'BNAME' and 'NAME'
                if 'PSRJ' in psr.keys():
                    if 'JNAME' not in psr.keys():
                        query_list[i]['JNAME'] = psr['PSRJ']
                    
                        if 'NAME' not in psr.keys():
                            query_list[i]['NAME'] = psr['PSRJ']

                if 'PSRB' in psr.keys():
                    query_list[i]['BNAME'] = psr['PSRB']

                    if 'NAME' not in query_list.keys():
                        query_list[i]['NAME'] = psr['PSRB']

            # convert query list to a pandas DataFrame
            try:
                df = DataFrame(query_list)
            except RuntimeError:
                raise RuntimeError("Could not convert list to DataFrame")

            # convert the DataFrame to an astropy table
            self.__table = Table.from_pandas(df)
            
            # add units if known
            for key in PSR_ALL_PARS:
                if key in self.__table.colnames:
                    if PSR_ALL[key]['units']:
                        self.__table.columns[key].unit = PSR_ALL[key]['units']

                    if PSR_ALL[key]['err'] and key+'_ERR' in self.__table.colnames:
                        self.__table.columns[key+'_ERR'].unit = PSR_ALL[key]['units']

            # add catalogue version to metadata
            self.__table.meta['version'] = self.get_version
            self.__table.meta['ATNF Pulsar Catalogue'] = ATNF_BASE_URL

        if isinstance(self.__table, Table):
            if query_params is None:
                return self.__table
            else:
                if isinstance(query_params, string_types):
                    query_params = [query_params]
                elif not isinstance(query_params, list):
                    raise TypeError("query_params must be a string or list.")

                # convert to numpy array
                query_params = np.array(query_params)

                # check parameters are in table
                intab = np.array([par in self.__table.colnames for par in query_params])

                if not np.all(intab):
                    warnings.warn("Not all request parameters '{}' were in"
                                  " the table".format(query_params[~intab].tolist()))

                if not np.any(intab):
                    warnings.warn("No requested parameters were in the "
                                  "table")

                # return only requested parameters
                table = self.__table[query_params[intab].tolist()]

                # return given the condition
                expression = None
                if usecondition == True and isinstance(self._condition, string_types):
                    expression = self._condition
                elif isinstance(usecondition, string_types):
                    expression = usecondition

                if expression is not None:
                    table = condition(table, expression, self._exactmatch)

                if (useseparation and self._coord is not None and 'RAJ' in
                        table.colnames and 'DECJ' in table.colnames):
                    # apply sky coordinate constraint
                    catalog = SkyCoord(table['RAJ'], table['DECJ'],
                                       unit=(aunits.hourangle, aunits.deg))

                    # get seperations
                    d2d = self._coord.separation(catalog)  

                    # find seperations within required radius
                    catalogmsk = d2d < self._radius*aunits.deg

                    table = table[catalogmsk]

                return table
        else:
            raise TypeError("Table is not an astropy Table!")

    def get_pulsars(self):
        """
        Returns:
            :class:`psrqpy.pulsar.Pulsars`: the queried pulsars returned as a
            :class:`~psrqpy.pulsar.Pulsars` object, which is a dictionary of
            :class:`~psrqpy.pulsar.Pulsar` objects.
        """

        if not self._pulsars:
            from .pulsar import Pulsar, Pulsars

            self._pulsars = Pulsars()

            # add pulsars one by one
            psrtable = self.table
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

    def parse_conditions(self, psrtype=None, assoc=None, bincomp=None):
        """
        Parse a string of `conditions
        <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_,
        i.e., logical statements with which to apply to a catalogue query,
        e.g., ``condition = 'f0 > 2.5 && assoc(GC)'``, so that they are in the
        format required for the query URL.

        Args:
            psrtype (list, str): a list of strings, or single string, of
                conditions on the
                `type <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_
                of pulsars to return (logical AND will be used for any listed
                types)
            assoc (list, str): a list of strings, or single string, of
                conditions on the associations of pulsars to return (logical
                AND will be used for any listed associations)
            bincomp (list, str): a list of strings, or single string, of
                conditions on the
                `binary companion <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type>`_
                types of pulsars to return (logical AND will be used for any
                listed associations)
            exactmatch (bool): a boolean stating whether assciations and types
                given as the condition should be an exact match

        Returns:
            str: a string with the format required for use in
                :attr:`~psrqpy.config.QUERY_URL`

        """

        conditionparse = ''

        # add on any extra given pulsar types
        if psrtype is not None:
            if isinstance(psrtype, list):
                if len(psrtype) == 0:
                    raise Exception("No pulsar types in list")

                for p in psrtype:
                    if not isinstance(p, string_types):
                        raise Exception("Non-string value '{}' found in pulsar type list".format(p))
                self._query_psr_types = psrtype
            else:
                if isinstance(psrtype, string_types):
                    self._query_psr_types = [psrtype]
                else:
                    raise Exception("'psrtype' must be a list or string")

            for p in list(self._query_psr_types):
                if p.upper() not in PSR_TYPES:
                    warnings.warn("Pulsar type '{}' is not recognised, no type will be required".format(p))
                    self._query_psr_types.remove(p)
                else:
                    if len(conditionparse) == 0:
                        conditionparse = 'type({})'.format(p.upper())
                    else:
                        conditionparse += ' && type({})'.format(p.upper())

        # add on any extra given associations
        if assoc is not None:
            if isinstance(assoc, list):
                if len(assoc) == 0:
                    raise Exception("No pulsar types in list")

                for p in assoc:
                    if not isinstance(p, string_types):
                        raise Exception("Non-string value '{}' found in associations list".format(p))
                self._query_assocs = assoc
            else:
                if isinstance(assoc, string_types):
                    self._query_assocs = [assoc]
                else:
                    raise Exception("'assoc' must be a list or string")

            for p in list(self._query_assocs):
                if p.upper() not in PSR_ASSOC_TYPE:
                    warnings.warn("Pulsar association '{}' is not recognised, no type will be required".format(p))
                    self._query_assocs.remove(p)
                else:
                    if len(conditionparse) == 0:
                        conditionparse = 'assoc({})'.format(p.upper())
                    else:
                        conditionparse += ' && assoc({})'.format(p.upper())

        # add on any extra given binary companion types
        if bincomp is not None:
            if isinstance(bincomp, list):
                if len(assoc) == 0:
                    raise Exception("No pulsar types in list")

                for p in bincomp:
                    if not isinstance(p, string_types):
                        raise Exception("Non-string value '{}' found in binary companions list".format(p))
                self._query_bincomps = bincomp
            else:
                if isinstance(bincomp, string_types):
                    self._query_bincomps = [bincomp]
                else:
                    raise Exception("'bincomp' must be a list or string")

            for p in list(self._query_bincomps):
                if p.upper() not in PSR_BINARY_TYPE:
                    warnings.warn("Pulsar binary companion '{}' is not recognised, no type will be required".format(p))
                    self._query_bincomps.remove(p)
                else:
                    if len(conditionparse) == 0:
                        conditionparse = 'bincomp({})'.format(p.upper())
                    else:
                        conditionparse += ' && bincomp({})'.format(p.upper())

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

        return str(self.table)

    def __repr__(self):
        """
        Returns:
            str: :func:`repr` method returns the repr method of an :class:`astropy.table.Table`.
        """

        return repr(self.table)

    def ppdot(self, intrinsicpdot=False, excludeGCs=False, showtypes=[],
              showGCs=False, showSNRs=False, markertypes={}, deathline=True,
              deathmodel='Ip', filldeath=True, filldeathtype={}, showtau=True,
              brakingidx=3, tau=None, showB=True, Bfield=None, pdotlims=None,
              periodlims=None, usecondition=True, rcparams={}):
        """
        Draw a lovely period vs period derivative diagram.

        Args:
            intrinsicpdot (bool): use the intrinsic period derivative corrected
                for the `Shklovskii effect <https://en.wikibooks.org/wiki/Pulsars_and_neutron_stars/Pulsar_properties#Pulse_period>`_
                rather than the observed value. Defaults to False.
            excludeGCs (bool): exclude globular cluster pulsars as their period
                derivatives can be contaminated by intra-cluster accelerations.
                Defaults to False.
            showtypes (list, str): a list of pulsar types to highlight with
                markers in the plot. These can contain any of the following:
                ``BINARY``, ``HE``, ``NRAD``, ``RRAT``, ``XINS``, ``AXP`` or
                ``SGR``, or ``ALL`` to show all types. Default to showing no
                types.
            showGCs (bool): show markers to denote the pulsars in globular
                clusters. Defaults to False.
            showSNRs (bool): show markers to denote the pulsars with supernova
                remnants associated with them. Defaults to False.
            markertypes (dict): a dictionary of marker styles and colors keyed
                to the pulsar types above
            deathline (bool): draw the pulsar death line. Defaults to True.
            deathmodel (str): the type of death line to draw based on the
                models in :func:`psrqpy.utils.death_line`. Defaults to
                ``'Ip'``.
            filldeath (bool): set whether to fill the pulsar graveyard under
                the death line. Defaults to True.
            filldeathtype (dict): a dictionary of keyword arguments for the
                fill style of the pulsar graveyard.
            showtau (bool): show lines for a selection of characteritic ages.
                Defaults to True, and shows lines for :math:`10^5` through to
                :math:`10^9` yrs with steps in powers of 10.
            brakingidx (int): a braking index to use for the calculation of the
                characteristic age lines. Defaults to 3 for magnetic dipole
                radiation.
            tau (list): a list of characteristic ages to show on the plot.
            showB (bool): show lines of constant magnetic field strength.
                Defaults to True, and shows lines for :math:`10^{10}` through
                to :math:`10^{14}` gauss with steps in powers of 10.
            Bfield (list): a list of magnetic field strengths to plot.
            periodlims (array_like): the [min, max] period limits to plot with
            pdotlims (array_like): the [min, max] pdot limits to plot with
            usecondition (bool): if True create the P-Pdot diagram only with
                pulsars that conform the the original query condition values.
                Defaults to True.
            rcparams (dict): a dictionary of :py:obj:`matplotlib.rcParams`
                setup parameters for the plot.

        Returns:
            :class:`matplotlib.figure.Figure`: the figure object
        """

        try:
            import matplotlib as mpl
            from matplotlib import pyplot as pl
        except ImportError:
            raise Exception('Cannot produce P-Pdot plot as Matplotlib is not '
                            'available')

        if self._webform:
            raise Exception("Please repeat query with 'webform=False'")

        # get table containing all required parameters
        table = self.table(condition=usecondition,
                           query_params=['P0', 'P1', 'P1_I', 'ASSOC',
                                         'BINARY', 'TYPE'])

        if len(table) == 0:
            print("No pulsars found, so no P-Pdot plot has been produced")
            return None

        if isinstance(showtypes, string_types):
            nshowtypes = [showtypes]
        else:
            nshowtypes = showtypes

        for stype in list(nshowtypes):
            if 'ALL' == stype.upper():
                nshowtypes = list(PSR_TYPES)
                # remove radio as none are returned as this
                del nshowtypes[nshowtypes.index('RADIO')]
                break
            elif stype.upper() not in list(PSR_TYPES):
                warnings.warn('"TYPE" {} is not recognised, so will not be '
                              'included'.format(stype))
                del nshowtypes[nshowtypes.index(stype)]
            if 'SGR' == stype.upper():  # synonym for AXP
                nshowtypes[nshowtypes.index(stype)] = 'AXP'

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

        # extract periods and period derivatives
        periods = table['P0']
        pdots = table['P1']
        if intrinsicpdot:  # use instrinsic period derivatives if requested
            ipdotidx = np.isfinite(table['P1_I'])
            pdots[ipdotidx] = table['P1_I'][ipdotidx]

        # get only finite values
        pidx = (np.isfinite(periods)) & (np.isfinite(pdots))
        periods = periods[pidx]
        pdots = pdots[pidx]

        if 'ASSOC' in self._query_params:
            assocs = table['ASSOC'][pidx]     # associations
        if 'TYPE' in self._query_params:
            types = table['TYPE'][pidx]       # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = table['BINARY'][pidx]  # binary pulsars

        # now get only positive pdot values
        pidx = pdots > 0.
        periods = periods[pidx]
        pdots = pdots[pidx]
        if 'ASSOC' in self._query_params:
            assocs = assocs[pidx]      # associations
        if 'TYPE' in self._query_params:
            types = types[pidx]        # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = binaries[pidx]  # binary pulsars

        # check whether to exclude globular cluster pulsars that could have
        # contaminated spin-down value
        if excludeGCs:
            # use '!=' to find GC indexes
            nongcidxs = np.flatnonzero(np.char.find(assocs, 'GC:') == -1)
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
        if periodlims is None:
            periodlims = [10**np.floor(np.min(np.log10(periods))), 10.*int(np.ceil(np.max(pdots)/10.))]
        if pdotlims is None:
            pdotlims = [10**np.floor(np.min(np.log10(pdots))), 10**np.ceil(np.max(np.log10(pdots)))]
        ax.set_xlim(periodlims)
        ax.set_ylim(pdotlims)

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
        typelegstring['SNR'] = r'SNR'
        typelegstring['RRAT'] = r'RRAT'

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
                    typeidx = np.flatnonzero(np.char.find(binaries, '*') == -1)
                elif thistype in ['GC', 'SNR']:
                    typeidx = np.flatnonzero(np.char.find(assocs, thistype) != -1)
                else:
                    typeidx = np.flatnonzero(np.char.find(types, thistype) != -1)

                if len(typeidx) == 0:
                    continue

                # default to empty markers with no lines between them
                if 'markerfacecolor' not in markertypes[thistype]:
                    markertypes[thistype]['markerfacecolor'] = 'none'
                if 'linestyle' not in markertypes[thistype]:
                    markertypes[thistype]['linestyle'] = 'none'
                typehandle, = ax.loglog(periods[typeidx], pdots[typeidx],
                                        label=typelegstring[thistype],
                                        **markertypes[thistype])
                if thistype in typelegstring:
                    handles[typelegstring[thistype]] = typehandle
                else:
                    handles[thistype] = typehandle

                ax.legend(handles.values(), handles.keys(), loc='upper left', numpoints=1)

        # add characteristic age lines
        tlines = OrderedDict()
        if showtau:
            if tau is None:
                taus = [1e5, 1e6, 1e7, 1e8, 1e9]  # default characteristic ages
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
