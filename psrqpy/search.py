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

from pandas import DataFrame
from copy import deepcopy

from .config import *
from .utils import get_version, condition


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
    catalogue database file. The catalogue can be queried for specificpulsar
    parameters and for specific named pulsars. Conditions on the parameter can
    be specified. The results will be stored as a :class:`pandas.DataFrame`,
    but can also be accessed as an :class:`astropy.table.Table`.

    Args:
        params (str, :obj:`list`): a list of strings with the
            pulsar `parameters
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=expert#par_list>`_
            to query. The parameter names are case insensitive. If this is not
            given then all parameters will be returned by default.
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
        adsref (bool): Set if wanting to use an :class:`ads.search.SearchQuery`
            to get reference information. Defaults to False.
        loadfromdb (str): Load a pulsar database file from a given path rather
            than using the ATNF Pulsar Catalogue database. Defaults to None.
        loadquery (str): load an instance of :class:`~psrqpy.search.QueryATNF`
            from the given file, rather than performing a new query. This was
            `loadfromfile` in earlier versions, which still works but has been
            deprecated. Defaults to None.
        cache (bool): Cache the catalogue database file for future use. This is
            ignored if `loadfromdb` is given. Defaults to True.
        checkupdate (bool): If True then check whether a cached catalogue file
            has an update available, and re-download if there is an update.
            Defaults to False.
        frompandas (:class:`pandas.DataFrame`): create a new
            :class:`psrqpy.QueryATNF` object from an existing
            :class:`pandas.DataFrame`.
        fromtable (:class:`astropy.table.Table`): create a new
            :class:`psrqpy.QueryATNF` object from an existing
            :class:`astropy.table.Table`.
    """

    def __init__(self, params=None, condition=None, psrtype=None, assoc=None,
                 bincomp=None, exactmatch=False, sort_attr='jname',
                 sort_order='asc', psrs=None, include_errs=True,
                 include_refs=False, adsref=False, loadfromfile=None,
                 loadquery=None, loadfromdb=None, cache=True,
                 checkupdate=False, circular_boundary=None, coord1=None,
                 coord2=None, radius=0., frompandas=None, fromtable=None):
        if loadfromfile is not None and loadquery is None:
            loadquery = loadfromfile
        if loadquery:
            self.load(loadquery)
            return

        self.__dataframe = DataFrame()
        self._include_errs = include_errs
        self._include_refs = include_refs
        self._adsref = adsref
        self._savefile = None  # file to save class to
        self._loadfile = None  # file class loaded from
        self.condition = condition
        self.exactmatch = exactmatch
        self.psrs = psrs
        self._sort_order = sort_order
        self._sort_attr = sort_attr.upper()

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

        self._coord = None
        if self._coord1 and self._coord2 and self._radius != 0.:
            if params is None:
                params = []
            
            # set centre coordinate as an astropy SkyCoord
            coord = SkyCoord(self._coord1, self._coord2,
                             unit=(aunits.hourangle, aunits.deg))

        # set conditions
        condparse = self.parse_conditions(psrtype=psrtype, assoc=assoc,
                                          bincomp=bincomp)
        if len(condparse) > 0:
            if self._condition is None:
                self._condition = condparse
            else:
                self._condition += condparse

        self.query_params = params
        self._refs = None  # set of pulsar references
        self._pulsars = None  # gets set to a Pulsars object by get_pulsars()

        # store passed pandas DataFrame
        if isinstance(frompandas, DataFrame):
            self.__dataframe = frompandas.copy()
            return

        # store passed astropy Table
        if isinstance(fromtable, Table):
            self.__dataframe = fromtable.to_pandas()
            return

        # download and cache (if requested) the database file
        try:
            _ = self.get_catalogue(path_to_db=loadfromdb, cache=cache,
                                   update=checkupdate)
        except IOError:
            raise IOError("Could not get catalogue database file")

        # get references if required
        if self._include_refs:
            self.get_references(adsref, cache=self._cache)

        # perform requested sorting
        _ = self.sort(inplace=True)

    def get_references(self, useads=False, cache=True):

        from .utils import get_references

        self._useads = useads
        self._refs = None
        self._adsrefs = None

        if self._useads:
            self._refs, self._adsrefs = get_references(self._useads,
                                                       cache=cache)
        else:
            self._refs = get_references(False, cache=cache)

    def get_catalogue(self, path_to_db=None, cache=True, update=False,
                      overwrite=True):
        """
        Call the :func:`psrqpy.utils.get_catalogue` function to download the
        ATNF Pulsar Catalogue, or load a given catalogue path.

        Args:
            path_to_db (str): if the path to a local version of the database
            file is given then that will be read in rather than attempting to
            download the file (defaults to None).
        cache (bool): cache the downloaded ATNF Pulsar Catalogue file. Defaults
            to True. This is ignored if `path_to_db` is given.
        update (bool): if True the ATNF Pulsar Catalogue will be
            re-downloaded and cached if there has been a change compared to the
            currently cached version. This is ignored if `path_to_db` is given.
        overwrite (bool): if True the returned catalogue will overwrite the
            catalogue currently contained within the :class:`~psrqpy.QueryATNF`
            class. If False then a new :class:`~psrqpy.QueryATNF` copy of the
            catalogue will be returned.

        Returns:
            :class:`psrqpy.QueryATNF`: a table containing the catalogue.
        """
        from .utils import get_catalogue

        try:
            dbtable = get_catalogue(path_to_db=path_to_db, cache=cache,
                                    update=update, pandas=True)
        except RuntimeError:
            raise RuntimeError("Problem getting catalogue")

        if not overwrite:
            newcatalogue = QueryATNF(params=self.query_params,
                                     condition=self.condition,
                                     exactmatch=self.exactmatch,
                                     sort_attr=self._sort_attr,
                                     sort_order=self._sort_order,
                                     psrs=self.psrs,
                                     include_errs=self._include_errs,
                                     include_refs=self._include_refs,
                                     adsref=self._adsref, cache=False,
                                     coord1=self._coord1, coord2=self._cord2,
                                     radius=self._radius,
                                     frompandas=dbtable)
            return newcatalogue

        # update current catalogue
        self.__dataframe = DataFrame(dbtable)
        self._dbfile = path_to_db
        self._checkupdate = update
        self._cache = cache
        self._atnf_version = dbtable.version

        # calculate derived parameters
        self.set_derived()

        return self

    @property
    def columns(self):
        """
        Return the table column names.
        """

        return self.__dataframe.columns

    def sort(self, sort_attr='JNAME', sort_order='asc', inplace=False):
        """
        Sort the generated catalogue :class:`~astropy.table.Table` on a given
        attribute and in either ascending or descending order.

        Args:
            sort_attr (str): The parameter on which to perform the sorting of
                the query output. Defaults to 'JNAME'.
            sort_order (str): Set to 'asc' to sort the parameter values in
                ascending order, or 'desc' to sort in descending order.
                Defaults to ascending.
            inplace (bool): If True, and sorting the class' internal
                :class:`~pandas.DataFrame`, then the sorting will be done
                in place without returning a copy of the table, otherwise
                a sorted copy of the table will be returned.

        Returns:
            :class:`~pandas.DataFrame`: a table containing the sorted
                catalogue.
        """

        self._sort_attr = sort_attr.upper()

        if self._sort_attr not in self.columns:
            raise KeyError("Sorting by attribute '{}' is not possible as it "
                           "is not in the table".format(self._sort_attr))

        self._sort_order = sort_order

        # check sort order is either 'asc' or 'desc' (or some synonyms)
        if self._sort_order.lower() in ['asc', 'ascending', 'up', '^']:
            self._sort_order = 'asc'
        elif self._sort_order.lower() in ['desc', 'descending', 'down', 'v']:
            self._sort_order = 'desc'
        else:
            warnings.warn("Unrecognised sort order '{}', defaulting to "
                          "'ascending'".format(sort_order), UserWarning)
            self._sort_order = 'asc'

        sortorder = True if self._sort_order == 'asc' else False

        if inplace:
            # sort the stored dataframe
            _ = self.__dataframe.sort_values(self._sort_attr,
                                             ascending=sortorder,
                                             inplace=inplace)
            return self.__dataframe
        else:
            return self.__dataframe.sort_values(self._sort_attr,
                                                ascending=sortorder)

    def __getitem__(self, key):
        if key not in self.pandas.columns:
            raise KeyError("Key '{}' not in queried results".format(key))

        # return astropy table column
        return self.table[key]

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

    def as_array(self):
        """
        Returns:
            :class:`~numpy.ndarray`: the output table as an array.
        """

        return self.table.as_array()

    def __len__(self):
        return len(self.pandas)

    @property
    def psrs(self):
        """
        Return the name(s) of particular pulsars asked for in the query.
        """

        return self._psrs

    @psrs.setter
    def psrs(self, psrnames=None):
        """
        Set a list of names of pulsars to be returned by the query.

        Args:
            psrnames (str, list): a list of names, or a single name, of pulsars
                to be returned by the query.
        """

        # set the pulsar name list
        if psrnames is None:
            self._psrs = None
        else:
            if isinstance(psrnames, string_types):
                self._psrs = [psrnames]
            elif isinstance(psrnames, list):
                self._psrs = psrnames
            elif isinstance(psrnames, np.ndarray):
                if psrnames.dtype == np.str or psrnames.dtype == np.unicode:
                    self._psrs = psrnames.tolist()
                else:
                    raise TypeError("psrnames must be a list of strings")
            else:
                raise TypeError("psrnames must be a list of strings")

    @property
    def num_pulsars(self):
        """
        Return the number of pulsars found in with query
        """

        return len(self)

    @property
    def table(self):
        # convert to astropy table
        thistable = Table.from_pandas(self.pandas)

        if (self._coord is not None and 'RAJ' in thistable.colnames
                and 'DECJ' in thistable.colnames):
            # apply sky coordinate constraint
            catalog = SkyCoord(thistable['RAJ'], thistable['DECJ'],
                               unit=(aunits.hourangle, aunits.deg))

            # get seperations
            d2d = self._coord.separation(catalog)  

            # find seperations within required radius
            catalogmsk = d2d < self._radius*aunits.deg

            thistable = thistable[catalogmsk]

        # add units if known
        for key in PSR_ALL_PARS:
            if key in thistable.colnames:
                if PSR_ALL[key]['units']:
                    thistable.columns[key].unit = PSR_ALL[key]['units']

                    if (PSR_ALL[key]['err'] and 
                            key+'_ERR' in thistable.colnames):
                        thistable.columns[key+'_ERR'].unit = PSR_ALL[key]['units']

        # add catalogue version to metadata
        thistable.meta['version'] = self.get_version
        thistable.meta['ATNF Pulsar Catalogue'] = ATNF_BASE_URL

        return thistable

    @property
    def empty(self):
        """
        Return True if the :class:`pandas.DataFrame` containing the catalogue
        is empty.
        """

        return self.__dataframe.empty

    def query_table(self, query_params=None, usecondition=True,
                    useseparation=True):
        """
        Return an :class:`astropy.table.Table` from the query with new
        parameters or conditions if given.

        Args:
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

        if not self.empty:  # convert to Table if DataFrame is not empty
            if query_params is None:
                query_params = self.columns
            elif isinstance(query_params, string_types):
                query_params = [query_params]
            elif not isinstance(query_params, list):
                raise TypeError("query_params must be a string or list.")

            # convert to numpy array
            query_params = np.array(query_params)

            # check parameters are in table
            intab = np.array([par in self.columns for par in query_params])

            if not np.all(intab):
                warnings.warn("Not all request parameters '{}' were in the "
                              "table".format(query_params[~intab].tolist()))

            if not np.any(intab):
                warnings.warn("No requested parameters were in the table")

            # return given the condition
            expression = None
            if usecondition == True and isinstance(self._condition, string_types):
                expression = self._condition
            elif isinstance(usecondition, string_types):
                expression = usecondition

            # sort table
            dftable = self.sort(self._sort_attr, self._sort_order)
            if expression is not None:
                # apply conditions
                dftable = condition(dftable, expression, self._exactmatch)
            
            # return only requested parameters and convert to table
            table = Table.from_pandas(dftable[query_params[intab].tolist()])

            # add units if known
            for key in PSR_ALL_PARS:
                if key in table.colnames:
                    if PSR_ALL[key]['units']:
                        table.columns[key].unit = PSR_ALL[key]['units']

                    if PSR_ALL[key]['err'] and key+'_ERR' in table.colnames:
                        table.columns[key+'_ERR'].unit = PSR_ALL[key]['units']

            # add catalogue version to metadata
            table.meta['version'] = self.get_version
            table.meta['ATNF Pulsar Catalogue'] = ATNF_BASE_URL

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
            # return an empty table
            return Table()

    @property
    def condition(self):
        """
        Return the string of logical conditions applied to the pulsars.
        """

        return self._condition

    @condition.setter
    def condition(self, expression):
        """
        Set the logical condition string to apply to queried pulsars.

        Args:
            expression (str): A string containing logical expressions to apply
                to queried pulsars.
        """

        if not isinstance(expression, string_types) and expression is not None:
            raise TypeError("Condition must be a string")

        self._condition = expression

    @property
    def exactmatch(self):
        """
        Return the boolean stating whether certain conditions should apply an
        exact match.
        """

        return self._exactmatch

    @exactmatch.setter
    def exactmatch(self, match):
        """
        Set whether to apply an exact match criterion for certain conditions.

        Args:
            match (bool): A boolean stating whether or not to apply an exact
                match.
        """

        if not isinstance(match, bool):
            if isinstance(match, int):
                if match != 0 and match != 1:
                    raise TypeError("Exact match requires boolean")
            else:
                raise TypeError("Exact match requires boolean")

        self._exactmatch = bool(match)

    @property
    def query_params(self):
        """
        Return the parameters required for the query.
        """

        return self._query_params

    @query_params.setter
    def query_params(self, params):
        """
        Set the parameters with which to query from the catalogue.

        Args:
            params (list, str): A list of parameter names to query from the
               catalogue.
        """

        self._query_params = None

        if isinstance(params, list):
            if len(params) == 0:
                print('No query parameters have been specified')

            for p in params:
                if not isinstance(p, string_types):
                    raise TypeError("Non-string value '{}' found in params "
                                    "list".format(p))

            self._query_params = [p.upper() for p in params]
        elif isinstance(params, string_types):
            # make sure parameter is all upper case
            self._query_params = [params.upper()]
        elif params is not None:
            raise TypeError("'params' must be a list or string")

        # remove any duplicate
        if self._query_params is not None:
            self._query_params = list(set(self._query_params))
 
            for p in list(self._query_params):
                if p not in PSR_ALL_PARS:
                    warnings.warn("Parameter '{}' not recognised.".format(p),
                                  UserWarning)

    @property
    def catalogue(self):
        """
        Return a copy of the entire stored :class:`~pandas.DataFrame` catalogue
        without any sorting or conditions applied.  
        """

        return self.__dataframe.copy()

    @property
    def dataframe(self):
        """
        Return the query table as a :class:`pandas.DataFrame`.
        """

        return self.pandas

    @property
    def pandas(self):
        """
        Return the query table as a :class:`pandas.DataFrame`.
        """

        # get only required parameters and sort
        dftable = self.sort(self._sort_attr, self._sort_order)

        if self._condition is not None:
            # apply condition
            dftable = condition(dftable, self._condition, self._exactmatch) 

        # return only requested pulsars
        if self.psrs is not None:
            jnames = np.zeros(len(dftable), dtype=np.bool)
            if 'JNAME' in dftable.columns:
                jnames = np.array([psr in self.psrs
                                   for psr in dftable['JNAME']])

            bnames = np.zeros(len(dftable), dtype=np.bool)
            if 'BNAME' in dftable.columns:
                bnames = np.array([psr in self.psrs
                                   for psr in dftable['BNAME']])

            if np.any(jnames) and np.any(bnames):
                allnames = jnames | bnames
            elif np.any(jnames):
                allnames = jnames
            elif np.any(bnames):
                allnames = bnames
            else:
                warnings.warn("No requested pulsars '{}' were "
                              "found.".format(self.psrs), UserWarning)

            dftable = dftable[allnames]

        # return only the required query parameters
        if isinstance(self.query_params, list):
            retpars = list(self.query_params)  # return parameters
            
            for par in self.query_params:
                if par in PSR_ALL_PARS:
                    if PSR_ALL[par]['err'] and self._include_errs:
                        retpars.append(par+'_ERR')

                    if PSR_ALL[par]['ref'] and self._include_refs:
                        retpars.append(par+'_REF')

                        if self._useads and self._adsref is not None:
                            retpars.append(p+'_REFURL')

            retpars = list(set(retpars))  # remove duplicates

            dftable = dftable[retpars]

        # convert reference tags to reference strings
        if self._include_refs and isinstance(self._refs, dict):
            for par in dftable.columns:
                if par[-4:] == '_REF':
                    for i in dftable[par].index.values:
                        reftag = dftable[par][i]
                        if reftag in self._refs:
                            dftable.loc[i, par] = self._refs[reftag]

                            if self._useads and reftag in self._adsref:
                                dftable[par+'_REFURL'] = self._adsref[reftag]

        # reset the indices to zero in the dataframe
        return dftable.reset_index(drop=True)

    def set_derived(self):
        """
        Compute any derived parameters and add them to the class.

        These calculations are based on those in the `readCatalogue.c` and
        `defineParameters.c` files from the `PSRCAT`
        `code <http://www.atnf.csiro.au/research/pulsar/psrcat/download.html>`_.
        """

        self.define_dist()      # define the DIST and DIST1 parameters
        self.derived_p0()       # derive P0 from F0 if not given
        self.derived_f0()       # derive F0 from P0 if not given
        self.derived_p1()       # derive P1 from F1 if not given
        self.derived_f1()       # derive F1 from P1 if not given
        self.derived_pb()       # derive binary period from FB0
        self.derived_pbdot()    # derive Pbdot from FB1
        self.derived_fb0()      # derive orbital frequency from period
        self.derived_fb1()      # derive FB1 from PBDOT
        self.derived_age()      # characteristic age
        self.derived_bsurf()    # surace magnetic field
        self.derived_b_lc()     # magnetic field at light cylinder
        self.derived_edot()     # spin-down luminosity
        self.derived_edotd2()   # spin-down flux at Sun
        self.derived_pmtot()    # total proper motion
        self.derived_vtrans()   # transverse velocity
        self.derived_p1_i()     # instrinsic period derivative
        self.derived_age_i()    # intrinsic age
        self.derived_bsurf_i()  # intrinsic Bsurf
        self.derived_edot_i()   # intrinsic luminosity
        self.derived_flux()     # radio flux

    def define_dist(self):
        """
        Set the `DIST` and `DIST1` parameters using other values.
        """

        if not np.all([p in self.__dataframe.columns for p in ['PX', 'PX_ERR', 'DIST_A', 'DIST_AMN', 'DIST_AMX', 'DIST_DM', 'DIST_DM1']]):
            warnings.warn("Could not set distances.",
                          UserWarning)
            return

        PX = self.__dataframe['PX']
        PXERR = self.__dataframe['PX_ERR']
        DIST_A = self.__dataframe['DIST_A']
        DIST_AMN = self.__dataframe['DIST_AMN']
        DIST_AMX = self.__dataframe['DIST_AMX']
        DIST_DM = self.__dataframe['DIST_DM']
        DIST_DM1 = self.__dataframe['DIST_DM1']

        # DIST defaults to DM distance
        DIST = DIST_DM.copy()

        # DIST1 defaults to DM1 distance
        DIST1 = DIST_DM1.copy()

        ONEAU = 149597870.  # AU in km (from psrcat.h)
        ONEPC = 30.857e12   # 1 pc in km (from psrcat.h)

        idxpx = np.isfinite(PX) & np.isfinite(PXERR)

        # set distances using parallax if parallax has greater than 3 sigma significance
        pxsigma = np.zeros(len(PX))
        pxsigma[idxpx] = np.abs(PX[idxpx])/PXERR[idxpx]

        # use DIST_A if available
        idxdista = np.isfinite(DIST_A)

        DIST[idxdista] = DIST_A[idxdista]
        DIST1[idxdista] = DIST_A[idxdista]

        # indexes of parallaxes with greater than 3 sigma significance
        idxpxgt3 = (pxsigma > 3.) & ~np.isfinite(DIST_A)

        DIST[idxpxgt3] = (ONEAU/ONEPC)*(60.*60.*180)/(PX[idxpxgt3]*np.pi)
        DIST1[idxpxgt3] = (ONEAU/ONEPC)*(60.*60.*180)/(PX[idxpxgt3]*np.pi) 

        # if dist_amn and dist_amx exist and dist_dm lies within boundary
        # then use dist_dm else use the closest limit to dist_dm
        # if dist_dm is not defined then use (dism_amn + dist_amx)/2
        idxdist = np.isfinite(DIST) & ~idxpxgt3
        idxdist1 = np.isfinite(DIST1) & ~idxpxgt3

        idxa = ~((DIST <= DIST_AMX) & (DIST >= DIST_AMN)) 
        idxa1 = ~((DIST1 <= DIST_AMX) & (DIST1 >= DIST_AMN))

        DIST[idxa & idxdist & (DIST >= DIST_AMX)] = DIST_AMX[idxa & idxdist & (DIST >= DIST_AMX)]
        DIST1[idxa1 & idxdist1 & (DIST1 >= DIST_AMX)] = DIST_AMX[idxa1 & idxdist1 & (DIST1 >= DIST_AMX)]

        DIST[idxa & idxdist & (DIST < DIST_AMX)] = DIST_AMN[idxa & idxdist & (DIST < DIST_AMX)]
        DIST1[idxa1 & idxdist1 & (DIST1 < DIST_AMX)] = DIST_AMN[idxa1 & idxdist1 & (DIST1 < DIST_AMX)]

        idxdist = ~np.isfinite(DIST) & ~idxpxgt3 & np.isfinite(DIST_AMN) & np.isfinite(DIST_AMX)
        idxdist1 = ~np.isfinite(DIST) & ~idxpxgt3 & np.isfinite(DIST_AMN) & np.isfinite(DIST_AMX)

        DIST[idxdist] = 0.5*(DIST_AMN[idxdist] + DIST_AMX[idxdist])
        DIST1[idxdist1] = 0.5*(DIST_AMN[idxdist1] + DIST_AMX[idxdist1])

        self.__dataframe['DIST'] = DIST
        self.__dataframe['DIST1'] = DIST1 

    def derived_p0(self):
        """
        Calculate the period from the frequency in cases where period is not
        given.
        """

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'F0']]):
            warnings.warn("Could not set periods.",
                          UserWarning)
            return

        F0 = self.__dataframe['F0']
        P0 = self.__dataframe['P0']

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'F0']]):
            warnings.warn("Could not set periods.", UserWarning)
            return

        F0 = self.__dataframe['F0']
        P0 = self.__dataframe['P0']

        # find indices where P0 needs to be set from F0
        idxp0 = ~np.isfinite(P0) & np.isfinite(F0)
        P0new = P0.copy()
        P0new[idxp0] = 1./F0[idxp0]
        self.__dataframe.update(P0new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['P0_REF', 'F0_REF']]):
            P0REFnew = self.__dataframe['P0_REF'].copy()
            F0REF = self.__dataframe['F0_REF']
            P0REFnew[idxp0] = F0REF[idxp0]
            self.__dataframe.update(P0REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['P0_ERR', 'F0_ERR']]):
            P0ERRnew = self.__dataframe['P0_ERR'].copy()
            F0ERR = self.__dataframe['F0_ERR']
            P0ERRnew[idxp0] = F0ERR[idxp0]*P0[idxp0]**2
            self.__dataframe.update(P0ERRnew)

    def derived_f0(self):
        """
        Calculate the frequency from the period in cases where frequency is not
        given.
        """

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'F0']]):
            warnings.warn("Could not set periods.", UserWarning)
            return

        F0 = self.__dataframe['F0']
        P0 = self.__dataframe['P0']

        # find indices where F0 needs to be set from P0
        idxf0 = np.isfinite(P0) & ~np.isfinite(F0)
        F0new = F0.copy()
        F0new[idxf0] = 1./P0[idxf0]
        self.__dataframe.update(F0new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['P0_REF', 'F0_REF']]):
            F0REFnew = self.__dataframe['F0_REF'].copy()
            P0REF = self.__dataframe['P0_REF']
            F0REFnew[idxf0] = P0REF[idxf0]
            self.__dataframe.update(F0REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['P0_ERR', 'F0_ERR']]):
            F0ERRnew = self.__dataframe['F0_ERR'].copy()
            P0ERR = self.__dataframe['P0_ERR']
            F0ERRnew[idxf0] = P0ERR[idxf0]*F0[idxf0]**2
            self.__dataframe.update(F0ERRnew)

    def derived_p1(self):
        """
        Calculate the period derivative from the frequency derivative in cases
        where period derivative is not given.
        """

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'F0', 'F1', 'P1']]):
            warnings.warn("Could not set period derivatives.",
                          UserWarning)
            return

        F0 = self.__dataframe['F0']
        P0 = self.__dataframe['P0']
        F1 = self.__dataframe['F1']
        P1 = self.__dataframe['P1']

        # find indices where P0 needs to be set from F0
        idxp1 = ~np.isfinite(P1) & np.isfinite(F1)
        P1new = P1.copy()
        P1new[idxp1] = -(P0[idxp1]**2)*F1[idxp1]
        self.__dataframe.update(P1new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['P1_REF', 'F1_REF']]):
            P1REFnew = self.__dataframe['P1_REF'].copy()
            F1REF = self.__dataframe['F1_REF']
            P1REFnew[idxp1] = F1REF[idxp1]
            self.__dataframe.update(P1REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['P0_ERR', 'F0_ERR', 'F1_ERR', 'P1_ERR']]):
            P1ERRnew = self.__dataframe['P0_ERR'].copy()
            F1ERR = self.__dataframe['F1_ERR']
            F0ERR = self.__dataframe['F0_ERR']
            P1ERRnew[idxp1] = np.sqrt((P0[idxp1]**2*F1ERR[idxp1])**2 +
						          (2.0*P0[idxp1]**3*F1[idxp1]*F0ERR[idxp1])**2)
            self.__dataframe.update(P1ERRnew)

    def derived_f1(self):
        """
        Calculate the frequency derivative from the period derivative in cases
        where frequency derivative is not given.
        """

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'F0', 'F1', 'P1']]):
            warnings.warn("Could not set period derivatives.",
                          UserWarning)
            return

        F0 = self.__dataframe['F0']
        P0 = self.__dataframe['P0']
        F1 = self.__dataframe['F1']
        P1 = self.__dataframe['P1']

        # find indices where P0 needs to be set from F0
        idxf1 = np.isfinite(P1) & ~np.isfinite(F1)
        F1new = F1.copy()
        F1new[idxf1] = -(F0[idxf1]**2)*P1[idxf1]
        self.__dataframe.update(F1new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['P1_REF', 'F1_REF']]):
            F1REFnew = self.__dataframe['F1_REF'].copy()
            P1REF = self.__dataframe['P1_REF']
            F1REFnew[idxf1] = P1REF[idxf1]
            self.__dataframe.update(F1REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['P0_ERR', 'F0_ERR', 'F1_ERR', 'P1_ERR']]):
            F1ERRnew = self.__dataframe['F0_ERR'].copy()
            P1ERR = self.__dataframe['P1_ERR']
            P0ERR = self.__dataframe['P0_ERR']
            F1ERRnew[idxf1] = np.sqrt((F0[idxf1]**2*P1ERR[idxf1])**2 +
						          (2.0*F0[idxf1]**3*P1[idxf1]*P0ERR[idxf1])**2)
            self.__dataframe.update(F1ERRnew)

    def derived_pb(self):
        """
        Calculate binary orbital period from orbital frequency.
        """

        if not np.all([p in self.__dataframe.columns for p in ['PB', 'FB0']]):
            warnings.warn("Could not set orbital period.",
                          UserWarning)
            return

        FB0 = self.__dataframe['FB0']
        PBnew = self.__dataframe['PB'].copy()

        idxpb = ~np.isfinite(PBnew) & np.isfinite(FB0)
        PBnew[idxpb] = 1./(FB0[idxpb]*86400.)
        self.__dataframe.update(PBnew)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['PB_REF', 'FB0_REF']]):
            PBREFnew = self.__dataframe['PB_REF'].copy()
            FB0REF = self.__dataframe['FB0_REF']
            PBREFnew[idxpb] = FB0REF[idxpb]
            self.__dataframe.update(PBREFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['PB_ERR', 'FB0_ERR']]):
            PBERRnew = self.__dataframe['PB_ERR'].copy()
            FB0ERR = self.__dataframe['FB0_ERR']
            PBERRnew[idxpb] = FB0ERR[idxpb]*PBnew[idxpb]**2*86400.
            self.__dataframe.update(PBERRnew)

    def derived_pbdot(self):
        """
        Calculate binary orbital period derivative from orbital frequency
        derivative.
        """

        if not np.all([p in self.__dataframe.columns for p in ['PBDOT', 'FB1', 'PB']]):
            warnings.warn("Could not set orbital period derivative.",
                          UserWarning)
            return

        FB1 = self.__dataframe['FB1']
        PB = self.__dataframe['PB']
        PBDOTnew = self.__dataframe['PBDOT'].copy()

        idxpbdot = ~np.isfinite(PBDOTnew) & np.isfinite(FB1)
        PBDOTnew[idxpbdot] = -(PB[idxpbdot]**2*FB1[idxpbdot])
        self.__dataframe.update(PBDOTnew)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['PBDOT_REF', 'FB1_REF']]):
            PBDOTREFnew = self.__dataframe['PBDOT_REF'].copy()
            FB1REF = self.__dataframe['FB1_REF']
            PBDOTREFnew[idxpbdot] = FB1REF[idxpbdot]
            self.__dataframe.update(PBDOTREFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['PBDOT_ERR', 'FB1_ERR', 'FB0_ERR']]):
            PBDOTERRnew = self.__dataframe['PBDOT_ERR'].copy()
            FB1ERR = self.__dataframe['FB1_ERR']
            FB0ERR = self.__dataframe['FB0_ERR']
            PBDOTERRnew[idxpbdot] = np.sqrt((PB[idxpbdot]**2*FB1ERR[idxpbdot])**2 +
						          (2.0*PB[idxpbdot]**3*FB1[idxpbdot]*FB0ERR[idxpbdot])**2)
            self.__dataframe.update(PBDOTERRnew)

    def derived_fb0(self):
        """
        Calculate orbital frequency from orbital period.
        """

        if not np.all([p in self.__dataframe.columns for p in ['PB', 'FB0']]):
            warnings.warn("Could not set orbital frequency.",
                          UserWarning)
            return

        PB = self.__dataframe['PB']
        FB0new = self.__dataframe['FB0'].copy()

        idxfb0 = ~np.isfinite(FB0new) & np.isfinite(PB)
        FB0new[idxfb0] = 1./(PB[idxfb0]*86400.)
        self.__dataframe.update(FB0new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['PB_REF', 'FB0_REF']]):
            FB0REFnew = self.__dataframe['FB0_REF'].copy()
            PBREF = self.__dataframe['PB_REF']
            FB0REFnew[idxfb0] = PBREF[idxfb0]
            self.__dataframe.update(FB0REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['PB_ERR', 'FB0_ERR']]):
            FB0ERRnew = self.__dataframe['FB0_ERR'].copy()
            PBERR = self.__dataframe['PB_ERR']
            FB0ERRnew[idxfb0] = PBERR[idxfb0]*(FB0new[idxfb0]**2)*86400.
            self.__dataframe.update(FB0ERRnew)

    def derived_fb1(self):
        """
        Calculate theorbital frequency derivative from the binary orbital
        period derivative.
        """

        if not np.all([p in self.__dataframe.columns for p in ['PBDOT', 'FB1', 'FB0']]):
            warnings.warn("Could not set orbital period derivative.",
                          UserWarning)
            return

        PBDOT = self.__dataframe['PBDOT']
        FB0 = self.__dataframe['FB0']
        FB1new = self.__dataframe['FB1'].copy()

        idxfb1 = ~np.isfinite(FB1new) & np.isfinite(PBDOT)
        FB1new[idxfb1] = -(FB0[idxfb1]**2*PBDOT[idxfb1])
        self.__dataframe.update(FB1new)

        # set the references
        if np.all([p in self.__dataframe.columns for p in ['PBDOT_REF', 'FB1_REF']]):
            FB1REFnew = self.__dataframe['FB1_REF'].copy()
            PBDOTREF = self.__dataframe['PBDOT_REF']
            FB1REFnew[idxfb1] = PBDOTREF[idxfb1]
            self.__dataframe.update(FB1REFnew)

        # set the errors
        if np.all([p in self.__dataframe.columns for p in ['PBDOT_ERR', 'FB1_ERR', 'PB_ERR']]):
            FB1ERRnew = self.__dataframe['FB1_ERR'].copy()
            PBDOTERR = self.__dataframe['PBDOT_ERR']
            PBERR = self.__dataframe['PB_ERR']
            FB1ERRnew[idxfb1] = np.sqrt((FB0[idxfb1]**2*PBDOTERR[idxfb1])**2 +
						          (2.0*FB0[idxfb1]**3*PBDOT[idxfb1]*PBERR[idxfb1]*86400.)**2)
            self.__dataframe.update(FB1ERRnew)

    def derived_p1_i(self):
        """
        Calculate the intrinsic period derivative.
        """

        if 'P1_I' in self.__dataframe.columns:
            return

        if 'VTRANS' not in self.__dataframe.columns:
            self.derived_vtrans()

        if not np.all([p in self.__dataframe.columns for p in ['VTRANS', 'P0', 'P1', 'DIST']]):
            warnings.warn("Could not set intrinsic period derivative.",
                          UserWarning)
            return

        # get required parameters
        VTRANS = self.__dataframe['VTRANS']
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']
        DIST = self.__dataframe['DIST']

        p1i = ((P1/1.0e-15) - VTRANS**2*1.0e10*P0/(DIST*3.086e6)/2.9979e10)*1.0e-15
        p1i[~np.isfinite(P1) | ~np.isfinite(P0) | ~np.isfinite(VTRANS) | ~np.isfinite(DIST)] = np.nan
        self.__dataframe['P1_I'] = p1i

    def derived_age(self):
        """
        Calculate the characteristic age.
        """

        # check if AGE is already defined
        if 'AGE' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1']]):
            warnings.warn("Could not set characteristic age.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']

        age = 0.5 * P0 / P1 / (60.0 * 60.0 * 24.0 * 365.25)
        age[(P1 < 0) | ~np.isfinite(P1) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['AGE'] = age

    def derived_age_i(self):
        """
        Calculate the characteristic age, dervied from period and intrinsic
        period derivative.
        """

        if 'AGE_I' in self.__dataframe.columns:
            return

        if 'P1_I' not in self.__dataframe.columns:
            self.derived_p1_i()

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1_I']]):
            warnings.warn("Could not set characteristic age.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1_I = self.__dataframe['P1_I']

        age_i = 0.5 * (P0 / P1_I) / (60.0 * 60.0 * 24.0 * 365.25)
        age_i[(P1_I < 0) | ~np.isfinite(P1_I) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['AGE_I'] = age_i

    def derived_bsurf(self):
        """
        Calculate the surface magnetic field strength.
        """

        if 'BSURF' in self.__dataframe.columns:
            return
        
        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1']]):
            warnings.warn("Could not set surface magnetic field.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']

        bsurf = 3.2e19 * np.sqrt(np.abs(P0 * P1))
        bsurf[(P1 < 0) | ~np.isfinite(P1) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['BSURF'] = bsurf

    def derived_bsurf_i(self):
        """
        Calculate the surface magnetic field strength, dervied from period and
        intrinsic period derivative.
        """

        if 'BSURF_I' in self.__dataframe.columns:
            return

        if 'P1_I' not in self.__dataframe.columns:
            self.derived_p1_i()

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1_I']]):
            warnings.warn("Could not set surface magnetic field.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1_I = self.__dataframe['P1_I']

        bsurf_i = 3.2e19 * np.sqrt(np.abs(P0 * P1_I))
        bsurf_i[(P1_I < 0) | ~np.isfinite(P1_I) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['BSURF_I'] = bsurf_i

    def derived_b_lc(self):
        """
        Calculate the magnetic field strength at the light cylinder.
        """

        if 'B_LC' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1']]):
            warnings.warn("Could not set light cylinder magnetic field.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']

        blc = 3.0e8*np.sqrt(np.abs(P1))*np.abs(P0)**(-5./2.)
        blc[(P1 < 0) | ~np.isfinite(P1) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['B_LC'] = blc

    def derived_edot(self):
        """
        Calculate the spin-down luminosity.
        """

        if 'EDOT' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1']]):
            warnings.warn("Could not set spin-down luminosity.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']

        edot = 4.0 * np.pi**2 * 1e45 * P1 / P0**3
        edot[(P1 < 0) | ~np.isfinite(P1) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['EDOT'] = edot

    def derived_edot_i(self):
        """
        Calculate the spin-down luminosity, dervied from period and intrinsic
        period derivative.
        """

        if 'EDOT_I' in self.__dataframe.columns:
            return

        if 'P1_I' not in self.__dataframe.columns:
            self.derived_p1_i()

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1_I']]):
            warnings.warn("Could not set spin-down luminosity.",
                          UserWarning)
            return

        # get period and period derivative
        P0 = self.__dataframe['P0']
        P1_I = self.__dataframe['P1_I']

        edot_i = 4.0 * np.pi**2 * 1e45 * P1_I / P0**3
        edot_i[(P1_I < 0) | ~np.isfinite(P1_I) | ~np.isfinite(P0)] = np.nan
        self.__dataframe['EDOT_I'] = edot_i

    def derived_edotd2(self):
        """
        Calculate the spin-down luminosity flux at the Sun.
        """

        if 'EDOTD2' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['P0', 'P1', 'DIST']]):
            warnings.warn("Could not set spin-down luminosity flux.",
                          UserWarning)
            return

        # get period, period derivative and distance
        P0 = self.__dataframe['P0']
        P1 = self.__dataframe['P1']
        DIST = self.__dataframe['DIST']

        edotd2 = 4.0 * np.pi**2 * 1e45 * (P1 / P0**3) / DIST**2
        edotd2[(P1 < 0) | ~np.isfinite(P1) | ~np.isfinite(P0) | ~np.isfinite(DIST)] = np.nan
        self.__dataframe['EDOTD2'] = edotd2

    def derived_pmtot(self):
        """
        Calculate the total proper motion and error.
        """
        
        if 'PMTOT' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['PMRA', 'PMDEC', 'PMELONG', 'PMELAT']]):
            warnings.warn("Could not set total proper motion.",
                          UserWarning)
            return

        # get PMRA and PMDEC
        PMRA = self.__dataframe['PMRA'].copy()
        PMDEC = self.__dataframe['PMDEC'].copy()
        PMELONG = self.__dataframe['PMELONG']
        PMELAT = self.__dataframe['PMELAT']

        # use PM ELONG or ELAT if no RA and DEC
        useelong = ~np.isfinite(PMRA) & np.isfinite(PMELONG)
        useelat = ~np.isfinite(PMDEC) & np.isfinite(PMELAT)
        PMRA[useelong] = PMELONG[useelong]
        PMDEC[useelat] = PMELAT[useelat]

        pmtot = np.sqrt(PMRA**2+PMDEC**2)
        self.__dataframe['PMTOT'] = pmtot

        # get the error
        if not np.all([p in self.__dataframe.columns for p in ['PMRA_ERR', 'PMDEC_ERR', 'PMELONG_ERR', 'PMELAT_ERR']]):
            return

        PMRA_ERR = self.__dataframe['PMRA_ERR'].copy()
        PMDEC_ERR = self.__dataframe['PMDEC_ERR'].copy()
        PMELONG_ERR = self.__dataframe['PMELONG_ERR']
        PMELAT_ERR = self.__dataframe['PMELAT_ERR']
        PMDEC_ERR[useelong] = PMELONG_ERR[useelong]
        PMRA_ERR[useelat] = PMELAT_ERR[useelat]

        pmtoterr = np.sqrt(((PMRA*PMRA_ERR)**2+(PMDEC*PMDEC_ERR)**2)/(PMRA**2 + PMDEC**2))
        self.__dataframe['PMTOT_ERR'] = pmtoterr

    def derived_vtrans(self):
        """
        Calculate the transverse velocity.
        """

        if 'VTRANS' in self.__dataframe.columns:
            return

        if 'PMTOT' not in self.__dataframe.columns:
            self.derived_pmtot()

        if not np.all([p in self.__dataframe.columns for p in ['PMTOT', 'DIST']]):
            warnings.warn("Could not set transverse velocity.",
                          UserWarning)
            return

        PMTOT = self.__dataframe['PMTOT']
        DIST = self.__dataframe['DIST']

        vtrans = (PMTOT/(1000.0*3600.0*180.0*np.pi*365.25*86400.0))*3.086e16*DIST
        vtrans[~np.isfinite(PMTOT) | ~np.isfinite(DIST)] = np.nan
        self.__dataframe['VTRANS'] = vtrans

    def derived_flux(self):
        """
        Calculate spectral index between 400 and 1400 MHz and radio
        flux at 400 and 1400 MHz.
        """

        if 'SI414' in self.__dataframe.columns:
            return

        if not np.all([p in self.__dataframe.columns for p in ['S1400', 'S400']]):
            warnings.warn("Could not set spectral index.",
                          UserWarning)
            return

        S1400 = self.__dataframe['S1400']
        S400 = self.__dataframe['S400']

        SI414 = -np.log10(S400/S1400)/(np.log10(400.0/1400.0))
        SI414[~np.isfinite(S1400) | ~np.isfinite(S400)] = np.nan
        self.__dataframe['SI414'] = SI414

        # need distance for flux
        if 'DIST' not in self.__dataframe.columns:
            self.define_dist()

        if 'DIST' not in self.__dataframe.columns:
            return

        DIST = self.__dataframe['DIST']

        R_LUM = S400 * DIST**2
        R_LUM[~np.isfinite(S400) | ~np.isfinite(DIST)] = np.nan
        self.__dataframe['R_LUM'] = R_LUM

        R_LUM14 = S1400 * DIST**2
        R_LUM14[~np.isfinite(S1400) | ~np.isfinite(DIST)] = np.nan
        self.__dataframe['R_LUM14'] = R_LUM14

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

                P = Pulsar(attrs['JNAME'], **attrs)
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

        return len(self.pandas)

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

        from .utils import death_line, label_line

        # get table containing all required parameters
        table = self.query_table(condition=usecondition,
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

        if 'ASSOC' in self.query_params:
            assocs = table['ASSOC'][pidx]     # associations
        if 'TYPE' in self.query_params:
            types = table['TYPE'][pidx]       # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = table['BINARY'][pidx]  # binary pulsars

        # now get only positive pdot values
        pidx = pdots > 0.
        periods = periods[pidx]
        pdots = pdots[pidx]
        if 'ASSOC' in self.query_params:
            assocs = assocs[pidx]      # associations
        if 'TYPE' in self.query_params:
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
            if 'ASSOC' in self.query_params:
                assocs = assocs[nongcidxs]
            if 'TYPE' in self.query_params:
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
