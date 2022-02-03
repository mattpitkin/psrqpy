"""
The classes defined here are for querying the `ATNF pulsar catalogue
<http://www.atnf.csiro.au/people/pulsar/psrcat/>`_ and viewing the resulting
information.
"""

import warnings
import re

import pickle

import numpy as np
import astropy
from astropy.coordinates import SkyCoord, ICRS, Galactic, BarycentricMeanEcliptic
from astropy.table.column import MaskedColumn, Column
import astropy.units as aunits
from astropy.constants import c, GM_sun
from astropy.table import Table
from packaging import version

from pandas import DataFrame, Series
from copy import deepcopy

from .config import ATNF_BASE_URL, PSR_ALL, PSR_ALL_PARS, PSR_TYPE, PSR_ASSOC_TYPE, PSR_BINARY_TYPE
from .utils import condition, age_pdot, B_field_pdot, h0_to_q22, q22_to_ellipticity


# set default astropy galactocentric frame values
# (https://docs.astropy.org/en/latest/coordinates/galactocentric.html)
if version.parse(astropy.__version__) >= version.parse("4.0"):
    from astropy.coordinates import galactocentric_frame_defaults
    _ = galactocentric_frame_defaults.set("pre-v4.0")


class QueryATNF(object):
    """
    A class to generate a query of the
    `ATNF pulsar catalogue <http://www.atnf.csiro.au/people/pulsar/psrcat/>`_.
    By default, this class will download and cache the latest version of the
    catalogue database file. The catalogue can be queried for specificpulsar
    parameters and for specific named pulsars. Conditions on the parameter can
    be specified. The results will be stored as a :class:`pandas.DataFrame`,
    but can also be accessed as an :class:`astropy.table.Table`.

    Args:
        params (str, :obj:`list`): a list of strings with the
            pulsar `parameters
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=expert#par_list>`_
            to query. The parameter names are case insensitive. If this is not
            given, then all parameters will be returned by default.
        condition (str): a string with logical conditions for the returned
            parameters. The allowed format of the condition string is given
            `here
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#condition>`_.
            Defaults to None.
        psrtype (str, :obj:`list`): a list of strings, or single string, of
            conditions on the `type
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_
            of pulsars to return (logical AND will be used for any listed
            types). Defaults to None.
        assoc (:obj:`list`, str): a condition on the associations of pulsars to
            return (logical AND will be used for any listed associations).
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
            format 'hh:mm:ss', the second entry is the centre point declination
            as a string in format 'dd:mm:ss', and the final entry is the
            circle's radius in degrees. This condition will only be applied if
            viewing the results as an :class:`astropy.table.Table`.
            Alternatively, `coord1`, `coord2`, and `radius` can be used.
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
        include_refs (bool): Set if wanting to include references tags in the
            output tables. Defaults to False.
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
        version (str): the version string (without the leading "v") of the ATNF
            catalogue version to download. This defaults to "latest" to get the
            most up-to-date version.
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
                 coord2=None, radius=0., frompandas=None, fromtable=None,
                 version="latest"):
        if loadfromfile is not None and loadquery is None:
            loadquery = loadfromfile
        if loadquery:
            self.load(loadquery)
            return

        self.__dataframe = DataFrame()
        self.include_errs = include_errs
        self._include_refs = include_refs
        self._savefile = None  # file to save class to
        self._loadfile = None  # file class loaded from
        self.condition = condition
        self.exactmatch = exactmatch
        self.psrs = psrs
        self._sort_order = sort_order
        self.sort_key = sort_attr.upper()
        self._useads = adsref

        # conditions for finding pulsars within a circular boundary
        self._coord1 = coord1
        self._coord2 = coord2
        self._radius = radius
        if isinstance(circular_boundary, (list, tuple)):
            if len(circular_boundary) != 3:
                raise ValueError("Circular boundary must contain three values")
            self._coord1 = circular_boundary[0]
            self._coord2 = circular_boundary[1]
            self._radius = circular_boundary[2]
        elif self._coord1 is None or self._coord2 is None or self._radius == 0.:
            # if any are not set then we can't define a boundary
            self._coord1 = self._coord2 = ''
            self._radius = 0.
        else:
            if (not isinstance(self._coord1, str)
                    or not isinstance(self._coord2, str)):
                raise ValueError("Circular boundary centre coordinates must "
                                 "be strings")
            if not isinstance(self._radius, float) and not isinstance(self._radius, int):
                raise ValueError("Circular boundary radius must be a float or "
                                 "int")

        self._coord = None
        if self._coord1 and self._coord2 and self._radius != 0.:
            if params is None:
                params = []

            # set centre coordinate as an astropy SkyCoord
            self._coord = SkyCoord(self._coord1, self._coord2,
                                   unit=(aunits.hourangle, aunits.deg))

        # set conditions
        condparse = self.parse_conditions(psrtype=psrtype, assoc=assoc,
                                          bincomp=bincomp)
        if len(condparse) > 0:
            if self.condition is None:
                self.condition = condparse
            else:
                self._condition += ' && {}'.format(condparse)

        self.query_params = params
        self._refs = None  # set of pulsar references
        self._pulsars = None  # gets set to a Pulsars object by get_pulsars()

        # store passed pandas DataFrame
        if isinstance(frompandas, DataFrame):
            self.__dataframe = frompandas.copy()

            # set version to None if not defined
            if not hasattr(self.__dataframe, 'version'):
                self.__dataframe.version = None
            return

        # store passed astropy Table
        if isinstance(fromtable, Table):
            self.__dataframe = fromtable.to_pandas()

            # set version if available
            if 'version' in fromtable.meta:
                self.__dataframe.version = fromtable.meta['version']
            else:
                self.__dataframe.version = None
            return

        # download and cache (if requested) the database file
        try:
            _ = self.get_catalogue(path_to_db=loadfromdb, cache=cache,
                                   update=checkupdate, version=version)
        except IOError:
            raise IOError("Could not get catalogue database file")

        # perform requested sorting
        _ = self.sort(inplace=True)

    def get_references(self, useads=False, cache=True):
        """
        Get a dictionary of short reference tags keyed to the full reference
        string. If requested also get a dictionary of reference tags keyed to
        NASA ADS URLs for the given reference. This uses the function
        :func:`psrqpy.utils.get_references`.

        Args:
            useads (bool): Set this to True to get the NASA ADS reference
                URLs. Defaults to False.
            cache (bool): The flag sets whether or not to use a pre-cached
                database of references. Defaults to True.
        """

        from .utils import get_references

        self._refs = None
        self._adsrefs = None
        useadst = useads or self._useads

        if useadst:
            self._refs, self._adsrefs = get_references(useadst,
                                                       cache=cache)
        else:
            self._refs = get_references(False, cache=cache)

    def parse_ref(self, refs, useads=False):
        """
        This function takes a short format reference string from the ATNF
        Pulsar Catalogue and returns the full format reference. It can
        also return a NASA ADS URL if requested and present.

        Args:
            refs (str, array_like): a single short reference string, or
                an array of reference strings.
            useads (bool): Set whether or not to also return a NASA ADS
                reference URL if present.

        Returns:
            array_like: a single full reference string, or an array of full
            reference strings. If NASA ADS URLs are requested, each return
            value will be a tuple containing the reference string and URL.
        """

        useadst = useads or self._useads

        if self._refs is None:
            self.get_references(useads=useadst)
        elif self._adsrefs is None and useadst:
            self.get_references(useads=useadst)

        singleref = False
        if isinstance(refs, str):
            singleref = True
            refs = [refs]

        # check refs is an iterable object
        if not hasattr(refs, '__iter__'):
            raise ValueError("Reference tags must be a string or array like")

        refstrs = []
        for ref in refs:
            if isinstance(ref, str):
                if ref in self._refs:
                    if useadst:
                        if ref in self._adsrefs:
                            refstrs.append((self._refs[ref], self._adsrefs[ref]))
                        else:
                            refstrs.append((self._refs[ref], None))
                    else:
                        refstrs.append(self._refs[ref])
                else:
                    if useadst:
                        refstrs.append((None, None))
                    else:
                        refstrs.append(None)
            else:
                if useadst:
                    refstrs.append((None, None))
                else:
                    refstrs.append(None)

        # just return a single value if only one input
        if singleref:
            return refstrs[0]

        return refstrs

    def get_catalogue(self, path_to_db=None, cache=True, update=False,
                      overwrite=True, version="latest"):
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
            version (str): the version string (without the leading "v") of the ATNF
                catalogue version to download. This defaults to "latest" to get the
                most up-to-date version.

        Returns:
            :class:`psrqpy.QueryATNF`: a table containing the catalogue.
        """
        from .utils import get_catalogue

        try:
            dbtable = get_catalogue(path_to_db=path_to_db, cache=cache,
                                    update=update, pandas=True,
                                    version=version)
        except Exception as e:
            raise RuntimeError("Problem getting catalogue: {}".format(str(e)))

        if not overwrite:
            newcatalogue = QueryATNF(params=self.query_params,
                                     condition=self.condition,
                                     exactmatch=self.exactmatch,
                                     sort_attr=self._sort_attr,
                                     sort_order=self._sort_order,
                                     psrs=self.psrs,
                                     include_errs=self._include_errs,
                                     include_refs=self._include_refs,
                                     adsref=self._useads, cache=False,
                                     coord1=self._coord1, coord2=self._coord2,
                                     radius=self._radius,
                                     frompandas=dbtable)
            return newcatalogue

        # update current catalogue
        self.__dataframe = DataFrame(dbtable)
        self.__dataframe.version = dbtable.version
        self._dbfile = path_to_db
        self._checkupdate = update
        self._cache = cache

        # calculate derived parameters
        self.set_derived()
        self.parse_types()

        return self

    @property
    def columns(self):
        """
        Return the table column names.
        """

        return self.__dataframe.columns

    def update(self, column, name=None, overwrite=False):
        """
        Update a column in the internal :class:`pandas.DataFrame` table using
        :meth:`pandas.DataFrame.update`. If the column does not exist, it will
        be added to the table.

        Args:
            column (:class:`pandas.Series`): a named column of values.
            name (str): the name of the column (required if `column` is not a
                :class:`pandas.Series`, or to overwrite the current column
                name)
            overwrite (bool): set whether to overwrite non-NA values or not if
               the column already exists. Defaults to False, so non-NA values
               will not be overwritten.
        """

        # get column name to update/add
        if name is not None:
            # use input `name` by default
            colname = name
        else:
            try:
                colname = column.name
            except AttributeError:
                colname = None

        if colname is not None:
            if colname in self.columns:
                if not isinstance(column, Series):
                    column = Series(column, name=colname)

                if column.dtype == self.__dataframe[colname].dtype:
                    self.__dataframe.update(column, overwrite=overwrite)
                else:
                    raise ValueError("Could not update table with "
                                     "supplied column")
            else:
                try:
                    self.catalogue[colname] = column
                except Exception as e:
                    raise ValueError("Could not add supplied columns to "
                                     "table: {}".format(str(e)))
        else:
            raise ValueError("No column name given")

    @property
    def sort_key(self):
        return self._sort_attr

    @sort_key.setter
    def sort_key(self, value):
        """
        Set the parameter to sort on.
        """

        if not isinstance(value, str):
            raise ValueError("Sort parameter must be a string")

        self._sort_attr = value

    def sort(self, sort_attr='JNAME', sort_order='asc', inplace=False):
        """
        Sort the generated catalogue :class:`~pandas.DataFrame` on a given
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

        if sort_attr is not None:
            self.sort_key = sort_attr.upper()

        if self.sort_key not in self.columns:
            raise KeyError("Sorting by attribute '{}' is not possible as it "
                           "is not in the table".format(self.sort_key))

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
            _ = self.__dataframe.sort_values(self.sort_key,
                                             ascending=sortorder,
                                             inplace=inplace)
            return self.__dataframe
        else:
            return self.__dataframe.sort_values(self.sort_key,
                                                ascending=sortorder)

    def __getitem__(self, key):
        if key in self.pandas.columns:
            # return astropy table column
            return self.table[key]
        else:
            psrrow = self.get_pulsar(key)

            if psrrow is None:
                raise KeyError("Key '{}' not in queried results".format(key))
            else:
                return psrrow

    def __getstate__(self):
        """
        Define to allow pickling of whole object.
        See, e.g., https://stackoverflow.com/a/2050357/1862861.
        """

        # Pulsars() object can cause pickling issues (in Python 2.7), so have
        # workaround
        from .pulsar import Pulsars
        if isinstance(self._pulsars, Pulsars):
            del self._pulsars
            self._pulsars = True

        # save ATNF version information from DataFrame separately
        self._atnf_version = self.catalogue.version

        return self.__dict__

    def __setstate__(self, d):
        """
        Define to allow pickling.
        """

        self.__dict__.update(d)

        # restore ATNF version info to catalogue
        self.__dataframe.version = self._atnf_version

        # Pulsars() object can cause pickling issues (in Python 2.7), so have
        # workaround
        if isinstance(self._pulsars, bool):
            from .pulsar import Pulsars

            if self._pulsars:
                # get the Pulsars() object
                self._pulsars = None
                _ = self.get_pulsars()

    def save(self, fname):
        """
        Output the :class:`~psrqpy.search.QueryATNF` instance to a pickle file
        for future loading.

        Args:
            fname (str): the filename to output the pickle to
        """

        try:
            fp = open(fname, 'wb')
            pickle.dump(self, fp, 2)
            fp.close()
            self._savefile = fname
        except IOError:
            raise IOError("Error outputing class to pickle file")

    def load(self, fname):
        """
        Load a previously saved pickle of this class.

        Args:
            fname (str): the filename of the pickled object
        """

        try:
            fp = open(fname, 'rb')
            tmpdict = pickle.load(fp)
            fp.close()
            self.__dict__.clear()  # clear current self
            self.__dict__.update(tmpdict.__dict__)
            self._loadfile = fname
        except IOError:
            raise IOError("Error reading in pickle")

    def as_array(self):
        """
        Returns:
            :class:`~numpy.ndarray`: the output table as an array.
        """

        return self.table.as_array()

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
            if isinstance(psrnames, str):
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
        """
        Return a :class:`astropy.table.Table` based on the query.
        """

        # convert to astropy table
        thistable = Table.from_pandas(self.pandas)

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

        return thistable

    @property
    def catalogue_table(self):
        """
        Return the full catalogue as a :class:`astropy.table.Table` without
        any query conditions applied.

        Note: in this returned table any references will not be converted into
        actual reference strings, but will still be the ATNF Pulsar Catalogue
        tags.
        """

        # convert catalogue to astropy table
        thistable = Table.from_pandas(self.catalogue)

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

        return thistable

    @property
    def empty(self):
        """
        Return True if the :class:`pandas.DataFrame` containing the catalogue
        is empty.
        """

        return self.__dataframe.empty

    def query_table(self, query_params=None, usecondition=True,
                    usepsrs=True, useseparation=True):
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
            usepsrs (bool, str): If True then the list of pulsars parsed to the
                :class:`psrqpy.QueryATNF`: class will be used when returning
                the table. If False then all pulsars in the catalogue will be
                returned. If a string, or list of strings, is given then that
                will be assumed to be a set of pulsar names to return.
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
            elif isinstance(query_params, str):
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
            if usecondition is True and isinstance(self.condition, str):
                expression = self.condition
            elif isinstance(usecondition, str):
                expression = usecondition

            # sort table
            dftable = self.sort(self.sort_key, self._sort_order)
            if expression is not None:
                # apply conditions
                dftable = condition(dftable, expression, self._exactmatch)

            # return only requested pulsars
            if usepsrs and self.psrs is not None:
                if not isinstance(usepsrs, bool):
                    # copy original pulsars
                    tmppsrs = deepcopy(self.psrs)

                    # set new pulsars
                    self.psrs = usepsrs

                jnames = np.zeros(len(dftable), dtype=bool)
                if "JNAME" in dftable.columns:
                    jnames = np.array([psr in self.psrs
                                      for psr in dftable["JNAME"]])

                bnames = np.zeros(len(dftable), dtype=bool)
                if "BNAME" in dftable.columns:
                    bnames = np.array([psr in self.psrs
                                       for psr in dftable["BNAME"]])

                if np.any(jnames) and np.any(bnames):
                    allnames = jnames | bnames
                elif np.any(jnames):
                    allnames = jnames
                elif np.any(bnames):
                    allnames = bnames
                else:
                    raise ValueError("No requested pulsars '{}' were "
                                     "found.".format(self.psrs), UserWarning)

                dftable = dftable[allnames]

                # reset original pulsar query
                if not isinstance(usepsrs, bool):
                    self.psrs = tmppsrs

            # return only requested parameters and convert to table
            table = Table.from_pandas(dftable[query_params[intab].tolist()])

            # add units if known
            for key in PSR_ALL_PARS:
                if key in table.colnames:
                    if PSR_ALL[key]["units"]:
                        table.columns[key].unit = PSR_ALL[key]["units"]

                    if PSR_ALL[key]["err"] and key+"_ERR" in table.colnames:
                        table.columns[key+"_ERR"].unit = PSR_ALL[key]["units"]

            # add catalogue version to metadata
            table.meta["version"] = self.get_version
            table.meta["ATNF Pulsar Catalogue"] = ATNF_BASE_URL

            if (useseparation and self._coord is not None and "RAJ" in
                    table.colnames and "DECJ" in table.colnames):
                # apply sky coordinate constraint
                catalog = SkyCoord(table["RAJ"], table["DECJ"],
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

        if not isinstance(expression, str) and expression is not None:
            raise TypeError("Condition must be a string")

        self._condition = expression

    @property
    def include_errs(self):
        """
        Return a boolean stating whether errors are to be included.
        """

        return self._include_errs

    @include_errs.setter
    def include_errs(self, inclerr):
        """
        Set whether to include errors with queried pulsars.

        Args:
            inclerr (bool): Set to True to include errors.
        """

        if isinstance(inclerr, (bool, int)):
            self._include_errs = bool(inclerr)
        else:
            TypeError("Flag must be boolean")

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
                if not isinstance(p, str):
                    raise TypeError("Non-string value '{}' found in params "
                                    "list".format(p))

            self._query_params = [p.upper() for p in params]
        elif isinstance(params, str):
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
        Return the entire stored :class:`~pandas.DataFrame` catalogue
        without any sorting or conditions applied.
        """

        return self.__dataframe

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
        dftable = self.sort(self.sort_key, self._sort_order)

        if (self._coord is not None and 'RAJD' in dftable.columns
                and 'DECJD' in dftable.columns):
            # apply sky coordinate constraint
            catalog = SkyCoord(dftable['RAJD'], dftable['DECJD'],
                               unit=(aunits.deg, aunits.deg))

            # get seperations
            d2d = self._coord.separation(catalog)

            # find seperations within required radius
            catalogmsk = d2d < self._radius*aunits.deg

            dftable = dftable[catalogmsk]

        if self._condition is not None:
            # apply condition
            dftable = condition(dftable, self._condition, self._exactmatch)

        # return only requested pulsars
        if self.psrs is not None:
            jnames = np.zeros(len(dftable), dtype=bool)
            if 'JNAME' in dftable.columns:
                jnames = np.array([psr in self.psrs
                                   for psr in dftable['JNAME']])

            bnames = np.zeros(len(dftable), dtype=bool)
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
                return DataFrame()  # empty dataframe

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

            retpars = list(set(retpars))  # remove duplicates

            dftable = dftable[retpars]

        # reset the indices to zero in the dataframe
        return dftable.reset_index(drop=True)

    def parse_types(self):
        """
        Parse information in 'ASSOC', 'TYPE', and 'BINCOMP', as described in
        `<http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_.
        """

        self.parse_assoc()      # parse the association parameter
        self.parse_type()       # parse the type parameter
        self.parse_bincomp()    # parse the binary companion parameter

    def parse_assoc(self):
        """
        Parse default string representing source associations, extracting (first) value
        and reference. Multiple values and references currently not supported.
        """

        if 'ASSOC' not in self.columns:
            return

        # Save original parameter in new column
        ASSOCorig = self.catalogue['ASSOC'].copy()
        ASSOCorig.name = 'ASSOC_ORIG'
        self.update(ASSOCorig, name='ASSOC_ORIG')

        ASSOCnew = self.catalogue['ASSOC'].copy()
        idxassoc = ~ASSOCnew.isna()

        # Set references first
        if 'ASSOC_REF' not in self.columns:
            ASSOCREFnew = Series(np.full(self.catalogue_len, '', dtype='U64'),
                                 name='ASSOC_REF')
        else:
            ASSOCREFnew = self.catalogue['ASSOC_REF'].copy()

        ASSOCREFnew[idxassoc] = ASSOCnew[idxassoc].apply(lambda x:
                                                         re.split(r'\]', re.split(r'\[', x)[1])[0]
                                                         if len(re.split(r'\[', x)) > 1 else np.nan)

        # Set values
        ASSOCnew[idxassoc] = ASSOCnew[idxassoc].apply(lambda x: re.split(r'\[|,|\(|:', x)[0])

        self.update(ASSOCnew, overwrite=True)
        self.update(ASSOCREFnew, overwrite=True)

    def parse_type(self):
        """
        Parse default string representing source type, extracting (first) value
        and reference. Multiple values and references currently not supported.
        """

        if 'TYPE' not in self.columns:
            return

        # Save original parameter in new column
        TYPEorig = self.catalogue['TYPE'].copy()
        TYPEorig.name = 'TYPE_ORIG'
        self.update(TYPEorig, 'TYPE_ORIG')

        TYPEnew = self.catalogue['TYPE'].copy()
        idxtype = ~TYPEnew.isna()

        # Set references first
        if 'TYPE_REF' not in self.columns:
            TYPEREFnew = Series(np.full(self.catalogue_len, '', dtype='U64'),
                                name='TYPE_REF')
        else:
            TYPEREFnew = self.catalogue['TYPE_REF'].copy()

        TYPEREFnew[idxtype] = TYPEnew[idxtype].apply(lambda x:
                                                     re.split(r'\]', re.split(r'\[', x)[1])[0]
                                                     if len(re.split(r'\[', x)) > 1 else np.nan)

        # Set values
        TYPEnew[idxtype] = TYPEnew[idxtype].apply(lambda x: re.split(r'\[|,|\(|:', x)[0])

        self.update(TYPEnew, overwrite=True)
        self.update(TYPEREFnew, overwrite=True)

    def parse_bincomp(self):
        """
        Parse default string representing source companion type, extracting (first) value
        and reference. Multiple values and references currently not supported.
        """

        if 'BINCOMP' not in self.columns:
            return

        # Save original parameter in new column
        BINCOMPorig = self.catalogue['BINCOMP'].copy()
        BINCOMPorig.name = 'BINCOMP_ORIG'
        self.update(BINCOMPorig)

        BINCOMPnew = self.catalogue['BINCOMP'].copy()
        idxbincomp = ~BINCOMPnew.isna()

        # Set references first
        if 'BINCOMP_REF' not in self.columns:
            BINCOMPREFnew = Series(np.full(self.catalogue_len, '', dtype='U64'),
                                   name='BINCOMP_REF')
        else:
            BINCOMPREFnew = self.catalogue['BINCOMP_REF'].copy()

        BINCOMPREFnew[idxbincomp] = BINCOMPnew[idxbincomp]\
            .apply(lambda x: re.split(r'\]', re.split(r'\[', x)[1])[0]
                   if len(re.split(r'\[', x)) > 1 else np.nan)

        # Set values
        BINCOMPnew[idxbincomp] = BINCOMPnew[idxbincomp]\
            .apply(lambda x: re.split(r'\[|,|\(|:', x)[0])

        self.update(BINCOMPnew, overwrite=True)
        self.update(BINCOMPREFnew, overwrite=True)

    def set_derived(self):
        """
        Compute any derived parameters and add them to the class.

        These calculations are based on those in the `readCatalogue.c` and
        `defineParameters.c` files from the `PSRCAT`
        `code <http://www.atnf.csiro.au/research/pulsar/psrcat/download.html>`_.
        """

        self.define_dist()         # define the DIST and DIST1 parameters
        self.derived_ecliptic()  # derive the ecliptic coordinates if not given
        self.derived_equatorial()  # derive equatorial coords from ecliptic
        self.define_galactic()     # define the galactic coordinates
        self.derived_p0()          # derive P0 from F0 if not given
        self.derived_f0()          # derive F0 from P0 if not given
        self.derived_p1()          # derive P1 from F1 if not given
        self.derived_f1()          # derive F1 from P1 if not given
        self.derived_pb()          # derive binary period from FB0
        self.derived_pbdot()       # derive Pbdot from FB1
        self.derived_fb0()         # derive orbital frequency from period
        self.derived_fb1()         # derive FB1 from PBDOT
        self.derived_age()         # characteristic age
        self.derived_bsurf()       # surface magnetic field
        self.derived_b_lc()        # magnetic field at light cylinder
        self.derived_edot()        # spin-down luminosity
        self.derived_edotd2()      # spin-down flux at Sun
        self.derived_pmtot()       # total proper motion
        self.derived_vtrans()      # transverse velocity
        self.derived_p1_i()        # instrinsic period derivative
        self.derived_age_i()       # intrinsic age
        self.derived_bsurf_i()     # intrinsic Bsurf
        self.derived_edot_i()      # intrinsic luminosity
        self.derived_flux()        # radio flux
        self.derived_binary()      # derived binary parameters
        self.derived_gw_h0_spindown_limit()  # derive GW parameters

    def define_dist(self):
        """
        Set the `DIST` and `DIST1` parameters using other values.
        """

        if 'PX' in self.columns:
            PX = self.catalogue['PX']
        else:
            PX = np.full(self.catalogue_len, np.nan)

        if 'PX_ERR' in self.columns:
            PXERR = self.catalogue['PX_ERR']
        else:
            PXERR = np.full(self.catalogue_len, np.nan)

        if 'DIST_A' in self.columns:
            DIST_A = self.catalogue['DIST_A']
        else:
            DIST_A = np.full(self.catalogue_len, np.nan)

        if 'DIST_AMN' in self.columns:
            DIST_AMN = self.catalogue['DIST_AMN']
        else:
            DIST_AMN = np.full(self.catalogue_len, np.nan)

        if 'DIST_AMX' in self.columns:
            DIST_AMX = self.catalogue['DIST_AMX']
        else:
            DIST_AMX = np.full(self.catalogue_len, np.nan)

        if 'DIST_DM' in self.columns:
            DIST_DM = self.catalogue['DIST_DM']
        else:
            DIST_DM = np.full(self.catalogue_len, np.nan)

        if 'DIST_DM1' in self.columns:
            DIST_DM1 = self.catalogue['DIST_DM1']
        else:
            DIST_DM1 = np.full(self.catalogue_len, np.nan)

        # DIST defaults to DM distance
        DIST = DIST_DM.copy()

        # DIST1 defaults to DM1 distance
        DIST1 = DIST_DM1.copy()

        ONEAU = 149597870.  # AU in km (from psrcat.h)
        ONEPC = 30.857e12   # 1 pc in km (from psrcat.h)

        idxpx = np.isfinite(PX) & np.isfinite(PXERR)

        # set distances using parallax if parallax has greater than 3 sigma
        # significance (and as long as parallax is positive)
        pxsigma = np.zeros(self.catalogue_len)
        pxsigma[idxpx] = PX[idxpx]/PXERR[idxpx]

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

        # ignore warnings from comparisons that involve NaN
        with np.errstate(invalid='ignore'):
            idxa = ~((DIST <= DIST_AMX) & (DIST >= DIST_AMN))
            idxa1 = ~((DIST1 <= DIST_AMX) & (DIST1 >= DIST_AMN))

            DIST[idxa & idxdist & (DIST >= DIST_AMX)] = DIST_AMX[idxa & idxdist &
                                                                 (DIST >= DIST_AMX)]
            DIST1[idxa1 & idxdist1 & (DIST1 >= DIST_AMX)] = DIST_AMX[(idxa1 & idxdist1
                                                                      & (DIST1 >= DIST_AMX))]

            DIST[idxa & idxdist & (DIST < DIST_AMX)] = DIST_AMN[idxa & idxdist & (DIST < DIST_AMX)]
            DIST1[idxa1 & idxdist1 & (DIST1 < DIST_AMX)] = DIST_AMN[(idxa1 & idxdist1
                                                                     & (DIST1 < DIST_AMX))]

        idxdist = (~np.isfinite(DIST) & ~idxpxgt3 &
                   np.isfinite(DIST_AMN) & np.isfinite(DIST_AMX))
        idxdist1 = (~np.isfinite(DIST) & ~idxpxgt3 &
                    np.isfinite(DIST_AMN) & np.isfinite(DIST_AMX))

        DIST[idxdist] = 0.5*(DIST_AMN[idxdist] + DIST_AMX[idxdist])
        DIST1[idxdist1] = 0.5*(DIST_AMN[idxdist1] + DIST_AMX[idxdist1])

        self.update(DIST, name='DIST')
        self.update(DIST1, name='DIST1')

    def derived_equatorial(self):
        """
        Calculate equatorial coordinates if only ecliptic coordinates are
        given. Errors on ecliptical coordinates are converted into
        equavalent errors on equatorial coordinates using a different
        algorithm to that used in `psrcat`.
        """

        reqpar = ['ELONG', 'ELAT']
        if not np.all([p in self.columns for p in reqpar]):
            return

        ELONG = self.catalogue['ELONG']
        ELAT = self.catalogue['ELAT']
        RAJDnew = np.full(self.catalogue_len, np.nan)
        DECJDnew = np.full(self.catalogue_len, np.nan)
        RAJnew = np.full(self.catalogue_len, '', dtype='U32')
        DECJnew = np.full(self.catalogue_len, '', dtype='U32')

        idx = np.isfinite(ELONG) & np.isfinite(ELAT)

        # get sky coordinates
        sc = BarycentricMeanEcliptic(
            ELONG.values[idx]*aunits.deg,
            ELAT.values[idx]*aunits.deg
        ).transform_to(ICRS())

        RAJDnew[idx] = sc.ra.value
        DECJDnew[idx] = sc.dec.value
        RAJnew[idx] = sc.ra.to('hourangle').to_string(sep=':', pad=True)
        DECJnew[idx] = sc.dec.to_string(sep=':', pad=True, alwayssign=True)

        self.update(RAJDnew, name='RAJD')
        self.update(DECJDnew, name='DECJD')
        self.update(RAJnew, name='RAJ')
        self.update(DECJnew, name='DECJ')

        reqpar = ['ELONG_ERR', 'ELAT_ERR']
        if np.all([p in self.columns for p in reqpar]):
            ELONG_ERR = self.catalogue['ELONG_ERR'].values.copy()
            ELAT_ERR = self.catalogue['ELAT_ERR'].values.copy()
            RAJD_ERRnew = np.full(self.catalogue_len, np.nan)
            DECJD_ERRnew = np.full(self.catalogue_len, np.nan)
            RAJ_ERRnew = np.full(self.catalogue_len, np.nan)
            DECJ_ERRnew = np.full(self.catalogue_len, np.nan)

            idxerr = idx & np.isfinite(ELONG_ERR) & np.isfinite(ELAT_ERR)

            # get position angle towards Northern ecliptic pole and rotate
            # ecliptical error ellipse to equatorial coordinates
            # (Jean Meeus, Astronomal Algorithms, 2nd edition, p. 100)
            ecl = 23.4392911 * aunits.deg
            l, b = ELONG.values[idxerr] * aunits.deg, ELAT.values[idxerr] * aunits.deg
            el, eb = ELONG_ERR[idxerr], ELAT_ERR[idxerr]
            elcosb = el * np.cos(b)
            q = np.arctan2(
                np.cos(l) * np.tan(ecl),
                np.sin(b) * np.sin(l) * np.tan(ecl) - np.cos(b)
            )
            cq, sq = np.cos(q), np.sin(q)
            eracosdec, edec = np.abs(elcosb * cq - eb * sq), np.abs(elcosb * sq + eb * cq)
            era = eracosdec / np.cos(np.deg2rad(DECJDnew[idxerr]))

            RAJD_ERRnew[idxerr] = era
            DECJD_ERRnew[idxerr] = edec
            RAJ_ERRnew[idxerr] = 3600 * era / 15
            DECJ_ERRnew[idxerr] = 3600 * edec

            self.update(RAJD_ERRnew, name='RAJD_ERR')
            self.update(DECJD_ERRnew, name='DECJD_ERR')
            self.update(RAJ_ERRnew, name='RAJ_ERR')
            self.update(DECJ_ERRnew, name='DECJ_ERR')

        # set references
        if 'ELONG_REF' in self.columns:
            RAJREFnew = np.full(self.catalogue_len, '', dtype='U32')
            DECJREFnew = np.full(self.catalogue_len, '', dtype='U32')
            ELONGREF = self.catalogue['ELONG_REF']

            DECJREFnew[idx] = ELONGREF[idx]
            RAJREFnew[idx] = ELONGREF[idx]

            self.update(DECJREFnew, name='RAJ_REF')
            self.update(RAJREFnew, name='DECJ_REF')

        # get PMRA and PMDEC if not given
        reqpar = ['PMELONG', 'PMELAT']
        if np.all([p in self.columns for p in reqpar]):
            PMELONG = self.catalogue['PMELONG']
            PMELAT = self.catalogue['PMELAT']
            PMRAnew = np.full(self.catalogue_len, np.nan)
            PMDECnew = np.full(self.catalogue_len, np.nan)

            idx = idx & np.isfinite(PMELONG) & np.isfinite(PMELAT)

            sc = BarycentricMeanEcliptic(
                ELONG[idx].values*aunits.deg,
                ELAT[idx].values*aunits.deg,
                pm_lon_coslat=PMELONG[idx].values*aunits.mas/aunits.yr,
                pm_lat=PMELAT[idx].values*aunits.mas/aunits.yr
                ).transform_to(ICRS())

            PMRAnew[idx] = sc.pm_ra_cosdec.value
            PMDECnew[idx] = sc.pm_dec.value

            self.update(PMRAnew, name='PMRA')
            self.update(PMDECnew, name='PMDEC')

    def derived_ecliptic(self):
        """
        Calculate the ecliptic coordinates, and proper motions, from the
        right ascension and declination if they are not already given.
        The ecliptic used here is the astropy's `BarycentricMeanEcliptic
        <http://docs.astropy.org/en/stable/api/astropy.coordinates.BarycentricMeanEcliptic.html>`_
        (note that this does not include nutation unlike the
        `BarycentricTrueEcliptic <http://docs.astropy.org/en/stable/api/astropy.coordinates.BarycentricTrueEcliptic.html>`_
        for which a bug fix was added in `astropy 3.2 <http://docs.astropy.org/en/v3.2.1/changelog.html#id12>`_),
        which may not exactly match that used in `psrcat`.
        """

        reqpar = ['RAJD', 'DECJD']
        if not np.all([p in self.columns for p in reqpar]):
            return

        RAJD = self.catalogue['RAJD']
        DECJD = self.catalogue['DECJD']
        ELONGnew = np.full(self.catalogue_len, np.nan)
        ELATnew = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(RAJD) & np.isfinite(DECJD)

        # get sky coordinates
        sc = SkyCoord(RAJD[idx].values*aunits.deg,
                      DECJD[idx].values*aunits.deg)

        ELONGnew[idx] = sc.barycentricmeanecliptic.lon.value
        ELATnew[idx] = sc.barycentricmeanecliptic.lat.value

        self.update(ELONGnew, name='ELONG')
        self.update(ELATnew, name='ELAT')

        # get references
        refpar = ['RAJ_REF', 'DECJ_REF']
        if np.all([p in self.columns for p in refpar]):
            RAJREF = self.catalogue['RAJ_REF']
            DECJREF = self.catalogue['DECJ_REF']
            ELONGREFnew = np.full(self.catalogue_len, '', dtype='U32')
            ELATREFnew = np.full(self.catalogue_len, '', dtype='U32')

            ELONGREFnew[idx] = RAJREF[idx]
            ELATREFnew[idx] = DECJREF[idx]

            self.update(ELONGREFnew, name='ELONG_REF')
            self.update(ELATREFnew, name='ELAT_REF')

        # get PMELONG and PMELAT if not given
        reqpar = ['PMRA', 'PMDEC']
        if np.all([p in self.columns for p in reqpar]):
            PMELONGnew = np.full(self.catalogue_len, np.nan)
            PMELATnew = np.full(self.catalogue_len, np.nan)
            PMRA = self.catalogue['PMRA']
            PMDEC = self.catalogue['PMDEC']

            idx = idx & np.isfinite(PMRA) & np.isfinite(PMDEC)

            sc = ICRS(
                RAJD[idx].values*aunits.deg,
                DECJD[idx].values*aunits.deg,
                pm_ra_cosdec=PMRA[idx].values*aunits.mas/aunits.yr,
                pm_dec=PMDEC[idx].values*aunits.mas/aunits.yr
            ).transform_to(BarycentricMeanEcliptic())

            PMELONGnew[idx] = sc.pm_lon_coslat.value
            PMELATnew[idx] = sc.pm_lat.value

            self.update(PMELONGnew, name='PMELONG')
            self.update(PMELATnew, name='PMELAT')

    def define_galactic(self):
        """
        Calculate the galactic longitude, latitude and position.

        .. note::
            The cartesian galactic coordinates returned by this function *do
            not* match those returned by the ATNF Pulsar Catalogue and the
            ``psrcat`` software. They are defined using the conventions in the
            :class:`astropy.coordinates.Galactocentric` class. This uses a
            Galactic centre distance of 8.3 kpc compared to 8.5 kpc in
            ``psrcat`` and rotated 90 degrees anticlockwise compared to
            ``psrcat``.

            The Galactic coordinate proper motions returned by this function
            *do not* match those returned by the ATNF Pulsar Catalogue and the
            ``psrcat`` software. The values returned here convert the observed
            proper motions in right ascension and declination (or elliptic
            longitude and latitude) into equivalent values in the Galactic
            coordinate system (via the :class:`astropy.coordinates.Galactic`
            class). However, the values returned by the ATNF Pulsar Catalogue
            and the ``psrcat`` software are in the Galactic cooridinate system,
            but additionally have the local solar system velocity and Galactic
            rotation of the pulsar removed from them as described in Section 3
            of `Harrison, Lyne & Anderson (1993) <https://ui.adsabs.harvard.edu/?#abs/1993MNRAS.261..113H>`_.
        """

        galpars = ['GL', 'GB', 'ZZ', 'XX', 'YY', 'DMSINB']
        if np.all([p in self.columns for p in galpars]):
            return

        reqpars = ['RAJD', 'DECJD']
        if not np.all([p in self.columns for p in reqpars]):
            return

        # get distance if required
        if 'DIST' not in self.columns:
            self.define_dist()

            if 'DIST' not in self.columns:
                return

        RAJD = self.catalogue['RAJD'].values.copy()
        DECJD = self.catalogue['DECJD'].values.copy()
        DIST = self.catalogue['DIST'].values.copy()

        GL = np.full(self.catalogue_len, np.nan)
        GB = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(RAJD) & np.isfinite(DECJD) & np.isfinite(DIST)

        # get sky coordinates
        sc = SkyCoord(RAJD[idx]*aunits.deg, DECJD[idx]*aunits.deg,
                      DIST[idx]*aunits.kpc)

        GL[idx] = sc.galactic.l.value
        GB[idx] = sc.galactic.b.value

        # set galactic longitude and latitude
        self.update(GL, name='GL')
        self.update(GB, name='GB')

        XX = np.full(self.catalogue_len, np.nan)
        YY = np.full(self.catalogue_len, np.nan)
        ZZ = np.full(self.catalogue_len, np.nan)

        # set galactocentric cartesian position (these seem to have a
        # different orientation (rotated 90 deg anticlockwise) to that
        # defined in the ATNF catalogue, and using a slightly different
        # distance to the galactic centre 8.3 kpc in astropy and 8.5 in psrcat)
        XX[idx] = sc.galactocentric.cartesian.x.value
        YY[idx] = sc.galactocentric.cartesian.y.value
        ZZ[idx] = sc.galactocentric.cartesian.z.value

        self.update(XX, name='XX')
        self.update(YY, name='YY')
        self.update(ZZ, name='ZZ')

        # set DMSINB
        if 'DM' in self.columns:
            DM = self.catalogue['DM']
            GB = self.catalogue['GB']

            DMSINB = np.full(self.catalogue_len, np.nan)

            idx = np.isfinite(GB) & np.isfinite(DM)
            DMSINB[idx] = DM[idx]*np.sin(np.deg2rad(GB[idx]))

            self.update(DMSINB, name='DMSINB')

        # galactic proper motion (in the Local Standard of Rest)
        if not np.all([p in self.columns for p in ['PMB', 'PML']]):
            reqpar = ['PMRA', 'PMDEC']
            if np.all([p in self.columns for p in reqpar]):
                PMRA = self.catalogue['PMRA'].values.copy()
                PMDEC = self.catalogue['PMDEC'].values.copy()

                PMB = np.full(self.catalogue_len, np.nan)
                PML = np.full(self.catalogue_len, np.nan)

                idx = (np.isfinite(PMRA) & np.isfinite(PMDEC) &
                       np.isfinite(RAJD) & np.isfinite(DECJD) &
                       np.isfinite(DIST))

                sc = ICRS(RAJD[idx]*aunits.deg, DECJD[idx]*aunits.deg,
                          distance=DIST[idx]*aunits.kpc,
                          pm_ra_cosdec=PMRA[idx]*aunits.mas/aunits.yr,
                          pm_dec=PMDEC[idx]*aunits.mas/aunits.yr
                          ).transform_to(Galactic())

                PMB[idx] = sc.pm_b.value
                PML[idx] = sc.pm_l_cosb.value

                self.update(PMB, name='PMB')
                self.update(PML, name='PML')

    def derived_binary(self):
        """
        Calculate derived binary system parameters.
        """

        MASS_PSR = 1.35  # canonical pulsar mass (solar masses)

        # derive mass function
        reqpars = ['A1', 'PB']
        if np.all([p in self.columns for p in reqpars]):
            A1 = self.catalogue['A1'].values.copy()*c.value  # convert to m
            PB = self.catalogue['PB'].values.copy()*86400.   # convert to sec

            idx = np.isfinite(A1) & np.isfinite(PB)

            MASSFN = np.full(self.catalogue_len, np.nan)
            MASSFN[idx] = (4.*np.pi**2/GM_sun.value)*A1[idx]**3/(PB[idx]**2)

            self.update(MASSFN, name='MASSFN')

            # derive minimum, median and 90% UL for mass
            MINMASS = np.full(self.catalogue_len, np.nan)
            MEDMASS = np.full(self.catalogue_len, np.nan)
            UPRMASS = np.full(self.catalogue_len, np.nan)
            from scipy.optimize import newton

            def solfunc(m2, sini, mf, m1):
                return (m1 + m2)**2 - (m2*sini)**3/mf

            SINI_MIN = 1.0  # inclination for minimum mass
            SINI_MED = 0.866025403  # inclination of 60deg for median mass
            SINI_90 = 0.438371146   # inclination for 90% UL mass
            for i, mf in enumerate(MASSFN):
                if ~np.isfinite(mf):
                    continue

                try:
                    MINMASS[i] = newton(solfunc, MASS_PSR,
                                        args=(SINI_MIN, mf, MASS_PSR),
                                        maxiter=1000)
                except RuntimeError:
                    MINMASS[i] = np.nan
                try:
                    MEDMASS[i] = newton(solfunc, MASS_PSR,
                                        args=(SINI_MED, mf, MASS_PSR),
                                        maxiter=1000)
                except RuntimeError:
                    MEDMASS[i] = np.nan
                try:
                    UPRMASS[i] = newton(solfunc, MASS_PSR,
                                        args=(SINI_90, mf, MASS_PSR),
                                        maxiter=1000)
                except RuntimeError:
                    UPRMASS[i] = np.nan

            self.update(MINMASS, name='MINMASS')
            self.update(MEDMASS, name='MEDMASS')
            self.update(UPRMASS, name='UPRMASS')

            # add uncertainty on mass function
            reqpars = ['A1_ERR', 'PB_ERR']
            if np.all([p in self.columns for p in reqpars]):
                # convert to metres
                A1ERR = self.catalogue['A1_ERR'].values.copy()*c.value
                # convert to seconds
                PBERR = self.catalogue['PB_ERR'].values.copy()*86400.

                idx = (np.isfinite(MASSFN) & np.isfinite(A1ERR) &
                       np.isfinite(PBERR))

                MASSFN_ERR = np.full(self.catalogue_len, np.nan)
                MASSFN_ERR[idx] = (
                    MASSFN[idx] * np.sqrt((3.*A1ERR[idx]/A1[idx])**2 +
                                          (2.*PBERR[idx]/PB[idx])**2)
                    )

                self.update(MASSFN_ERR, name='MASSFN_ERR')

        # derive eccentricity from EPS1 and EPS2
        reqpars = ['EPS1', 'EPS2']
        if np.all([p in self.columns for p in reqpars]):
            EPS1 = self.catalogue['EPS1'].values.copy()
            EPS2 = self.catalogue['EPS2'].values.copy()
            ECCnew = np.full(self.catalogue_len, np.nan)
            OMnew = np.full(self.catalogue_len, np.nan)

            idx = np.isfinite(EPS1) & np.isfinite(EPS2)

            # set eccentricities
            ECCnew[idx] = np.sqrt(EPS1[idx]**2+EPS2[idx]**2)

            self.update(ECCnew, name='ECC')

            # set angle of peristron
            idxn = idx & (ECCnew != 0.)
            OMnew[idxn] = np.arctan2(EPS1[idxn],
                                     EPS2[idxn])*180./np.pi

            with np.errstate(invalid='ignore'):
                OMnew = np.mod(OMnew+360., 360.)  # make sure angles are positive

            self.update(OMnew, name='OM')

            # set errors
            reqpars = ['EPS1_ERR', 'EPS2_ERR']
            if np.all([p in self.columns for p in reqpars]):
                EPS1ERR = self.catalogue['EPS1_ERR'].values.copy()
                EPS2ERR = self.catalogue['EPS2_ERR'].values.copy()
                ECCERRnew = np.full(self.catalogue_len, np.nan)
                OMERRnew = np.full(self.catalogue_len, np.nan)

                idxn = idx & (np.isfinite(EPS1ERR) & np.isfinite(EPS2ERR) &
                              (ECCnew != 0.))

                OMERRnew[idxn] = (
                    np.sqrt((EPS2[idxn]*EPS1ERR[idxn])**2
                            + (EPS1[idxn]*EPS2ERR[idxn])**2) /
                    (ECCnew[idxn])**2
                    )*180.0/np.pi
                self.update(OMERRnew, name='OM_ERR')

                ECCERRnew[idxn] = (
                    np.sqrt((EPS1[idxn]*EPS1ERR[idxn])**2
                            + (EPS2[idxn]*EPS2ERR[idxn])**2) / ECCnew[idxn]
                    )
                self.update(ECCERRnew, name='ECC_ERR')

        # derive EPS1 and EPS2 from ECC and OM
        reqpars = ['ECC', 'OM']
        if np.all([p in self.columns for p in reqpars]):
            ECC = self.catalogue['ECC'].values.copy()
            OM = self.catalogue['OM'].values.copy()
            EPS1new = np.full(self.catalogue_len, np.nan)
            EPS2new = np.full(self.catalogue_len, np.nan)

            idx = np.isfinite(ECC) & np.isfinite(OM)

            EPS1new[idx] = ECC[idx] * np.sin(OM[idx])
            EPS2new[idx] = ECC[idx] * np.cos(OM[idx])

            self.update(EPS1new, name='EPS1')
            self.update(EPS2new, name='EPS2')

            # set errors
            reqpars = ['ECC_ERR', 'OM_ERR']
            if np.all([p in self.columns for p in reqpars]):
                ECCERR = self.catalogue['ECC_ERR'].values.copy()
                OMERR = self.catalogue['OM_ERR'].values.copy()
                EPS1ERRnew = np.full(self.catalogue_len, np.nan)
                EPS2ERRnew = np.full(self.catalogue_len, np.nan)

                idxn = idx & np.isfinite(ECCERR) & np.isfinite(OMERR)

                EPS1ERRnew[idxn] = (
                    np.abs(EPS1new[idxn]) *
                    np.sqrt((ECCERR[idxn]/ECC[idxn])**2 +
                            (np.abs(np.cos(np.deg2rad(OM[idxn]))) * np.deg2rad(OMERR[idxn]) /
                             np.abs(np.sin(np.deg2rad(OM[idxn]))))**2)
                    )
                EPS2ERRnew[idxn] = (
                    np.abs(EPS2new[idxn]) *
                    np.sqrt((ECCERR[idxn]/ECC[idxn])**2 +
                            (np.abs(np.sin(np.deg2rad(OM[idxn]))) * np.deg2rad(OMERR[idxn]) /
                             np.abs(np.cos(np.deg2rad(OM[idxn]))))**2)
                    )

                self.update(EPS1ERRnew, name='EPS1_ERR')
                self.update(EPS2ERRnew, name='EPS2_ERR')

        # derive MINOMDOT
        reqpars = ['ECC', 'PB', 'MINMASS']
        if np.all([p in self.columns for p in reqpars]):
            MINMASS = self.catalogue['MINMASS'].values.copy()
            PB = self.catalogue['PB'].values.copy()*86400.
            ECC = self.catalogue['ECC'].values.copy()

            MINOMDOT = np.full(self.catalogue_len, np.nan)

            idx = np.isfinite(MINMASS) & np.isfinite(PB) & np.isfinite(ECC)

            MINOMDOT[idx] = (3.*(2.*np.pi/PB[idx])**(5./3.) *
                             ((MASS_PSR+MINMASS[idx]) *
                              4.925490946e-6)**(2./3.) /
                             (1.-ECC[idx]**2))
            MINOMDOT[idx] = np.rad2deg(MINOMDOT[idx])*86400.*365.25

            self.update(MINOMDOT, name='MINOMDOT')

    def derived_p0(self):
        """
        Calculate the period from the frequency in cases where period is not
        given.
        """

        if 'F0' not in self.columns:
            return

        F0 = self.catalogue['F0']
        P0new = np.full(self.catalogue_len, np.nan)

        # find indices where P0 needs to be set from F0
        idx = np.isfinite(F0)
        P0new[idx] = 1./F0[idx]
        self.update(P0new, name='P0')

        # set the references
        if 'F0_REF' in self.columns:
            P0REFnew = np.full(self.catalogue_len, '', dtype='U32')
            F0REF = self.catalogue['F0_REF']
            P0REFnew[idx] = F0REF[idx]
            self.update(P0REFnew, 'P0_REF')

        # set the errors
        if 'F0_ERR' in self.columns:
            P0ERRnew = np.full(self.catalogue_len, np.nan)
            F0ERR = self.catalogue['F0_ERR']
            idx = idx & np.isfinite(F0ERR)
            P0ERRnew[idx] = F0ERR[idx]*P0new[idx]**2
            self.update(P0ERRnew, 'P0_ERR')

    def derived_f0(self):
        """
        Calculate the frequency from the period in cases where frequency is not
        given.
        """

        if 'P0' not in self.columns:
            return

        P0 = self.catalogue['P0']
        F0new = np.full(self.catalogue_len, np.nan)

        # find indices where F0 needs to be set from P0
        idx = np.isfinite(P0)
        F0new[idx] = 1./P0[idx]
        self.update(F0new, name='F0')

        # set the references
        if 'P0_REF' in self.columns:
            F0REFnew = np.full(self.catalogue_len, '', dtype='U32')
            P0REF = self.catalogue['P0_REF']
            F0REFnew[idx] = P0REF[idx]
            self.update(F0REFnew, name='F0_REF')

        # set the errors
        if 'P0_ERR' in self.columns:
            F0ERRnew = np.full(self.catalogue_len, np.nan)
            P0ERR = self.catalogue['P0_ERR']
            idx = idx & np.isfinite(P0ERR)
            F0ERRnew[idx] = P0ERR[idx]*F0new[idx]**2
            self.update(F0ERRnew, name='F0_ERR')

    def derived_p1(self):
        """
        Calculate the period derivative from the frequency derivative in cases
        where period derivative is not given.
        """

        reqpars = ['P0', 'F1']
        if not np.all([p in self.columns for p in reqpars]):
            return

        P0 = self.catalogue['P0']
        F1 = self.catalogue['F1']
        P1new = np.full(self.catalogue_len, np.nan)

        # find indices where P0 needs to be set from F0
        idx = np.isfinite(P0) & np.isfinite(F1)
        P1new[idx] = -(P0[idx]**2)*F1[idx]
        self.update(P1new, name='P1')

        # set the references
        if 'F1_REF' in self.columns:
            P1REFnew = np.full(self.catalogue_len, '', dtype='U32')
            F1REF = self.catalogue['F1_REF']
            P1REFnew[idx] = F1REF[idx]
            self.update(P1REFnew, name='P1_REF')

        # set the errors
        reqpars = ['F0_ERR', 'F1_ERR']
        if np.all([p in self.columns for p in reqpars]):
            P1ERRnew = np.full(self.catalogue_len, np.nan)
            F1ERR = self.catalogue['F1_ERR']
            F0ERR = self.catalogue['F0_ERR']
            idx = idx & (np.isfinite(F1ERR) & np.isfinite(F0ERR))
            P1ERRnew[idx] = np.sqrt(
                (P0[idx]**2*F1ERR[idx])**2
                + (2.0*P0[idx]**3*F1[idx]*F0ERR[idx])**2)
            self.update(P1ERRnew, name='P1_ERR')

    def derived_f1(self):
        """
        Calculate the frequency derivative from the period derivative in cases
        where frequency derivative is not given.
        """

        reqpars = ['F0', 'P1']
        if not np.all([p in self.columns for p in reqpars]):
            return

        F0 = self.catalogue['F0']
        P1 = self.catalogue['P1']
        F1new = np.full(self.catalogue_len, np.nan)

        # find indices where P0 needs to be set from F0
        idx = np.isfinite(P1) & np.isfinite(F0)
        F1new[idx] = -(F0[idx]**2)*P1[idx]
        self.update(F1new, name='F1')

        # set the references
        if 'P1_REF' in self.columns:
            F1REFnew = np.full(self.catalogue_len, '', dtype='U32')
            P1REF = self.catalogue['P1_REF']
            F1REFnew[idx] = P1REF[idx]
            self.update(F1REFnew, name='F1_REF')

        # set the errors
        reqpars = ['P0_ERR', 'P1_ERR']
        if np.all([p in self.columns for p in reqpars]):
            F1ERRnew = np.full(self.catalogue_len, np.nan)
            P1ERR = self.catalogue['P1_ERR']
            P0ERR = self.catalogue['P0_ERR']
            idx = idx & np.isfinite(P1ERR) & np.isfinite(P0ERR)
            F1ERRnew[idx] = np.sqrt(
                (F0[idx]**2*P1ERR[idx])**2
                + (2.0*F0[idx]**3*P1[idx]*P0ERR[idx])**2)
            self.update(F1ERRnew, name='F1_ERR')

    def derived_pb(self):
        """
        Calculate binary orbital period from orbital frequency.
        """

        if 'FB0' not in self.columns:
            return

        FB0 = self.catalogue['FB0']
        PBnew = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(FB0)
        PBnew[idx] = 1./(FB0[idx]*86400.)
        self.update(PBnew, name='PB')

        # set the references
        if 'FB0_REF' in self.columns:
            PBREFnew = np.full(self.catalogue_len, '', dtype='U32')
            FB0REF = self.catalogue['FB0_REF']
            PBREFnew[idx] = FB0REF[idx]
            self.update(PBREFnew, name='PB_REF')

        # set the errors
        if 'FB0_ERR' in self.columns:
            PBERRnew = np.full(self.catalogue_len, np.nan)
            FB0ERR = self.catalogue['FB0_ERR']
            idx = idx & np.isfinite(FB0ERR)
            PBERRnew[idx] = FB0ERR[idx]*PBnew[idx]**2*86400.
            self.update(PBERRnew, name='PB_ERR')

    def derived_pbdot(self):
        """
        Calculate binary orbital period derivative from orbital frequency
        derivative.
        """

        reqpars = ['FB1', 'PB']
        if not np.all([p in self.columns for p in reqpars]):
            return

        FB1 = self.catalogue['FB1']
        PB = self.catalogue['PB']
        PBDOTnew = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(PB) & np.isfinite(FB1)
        PBDOTnew[idx] = -(PB[idx]**2*FB1[idx])
        self.update(PBDOTnew, name='PBDOT')

        # set the references
        if 'FB1_REF' in self.columns:
            PBDOTREFnew = np.full(self.catalogue_len, '', dtype='U32')
            FB1REF = self.catalogue['FB1_REF']
            PBDOTREFnew[idx] = FB1REF[idx]
            self.update(PBDOTREFnew, name='PBDOT_REF')

        # set the errors
        reqpars = ['FB1_ERR', 'FB0_ERR']
        if np.all([p in self.columns for p in reqpars]):
            PBDOTERRnew = np.full(self.catalogue_len, np.nan)
            FB1ERR = self.catalogue['FB1_ERR']
            FB0ERR = self.catalogue['FB0_ERR']
            idx = idx & np.isfinite(FB1ERR) & np.isfinite(FB0ERR)
            PBDOTERRnew[idx] = np.sqrt((PB[idx]**2 * FB1ERR[idx])**2
                                       + (2.0 * PB[idx]**3 * FB1[idx]
                                          * FB0ERR[idx])**2)
            self.update(PBDOTERRnew, name='PBDOT_ERR')

    def derived_fb0(self):
        """
        Calculate orbital frequency from orbital period.
        """

        if 'PB' not in self.columns:
            return

        PB = self.catalogue['PB']
        FB0new = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(PB)
        FB0new[idx] = 1./(PB[idx]*86400.)
        self.update(FB0new, name='FB0')

        # set the references
        if 'PB_REF' in self.columns:
            FB0REFnew = np.full(self.catalogue_len, '', dtype='U32')
            PBREF = self.catalogue['PB_REF']
            FB0REFnew[idx] = PBREF[idx]
            self.update(FB0REFnew, name='FB0_REF')

        # set the errors
        if 'PB_ERR' in self.columns:
            FB0ERRnew = np.full(self.catalogue_len, np.nan)
            PBERR = self.catalogue['PB_ERR']
            idx = idx & np.isfinite(PBERR)
            FB0ERRnew[idx] = PBERR[idx]*(FB0new[idx]**2)*86400.
            self.update(FB0ERRnew, name='FB0_ERR')

    def derived_fb1(self):
        """
        Calculate the orbital frequency derivative from the binary orbital
        period derivative.
        """

        reqpars = ['PBDOT', 'FB0']
        if not np.all([p in self.columns for p in reqpars]):
            return

        PBDOT = self.catalogue['PBDOT']
        FB0 = self.catalogue['FB0']
        FB1new = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(FB0) & np.isfinite(PBDOT)
        FB1new[idx] = -(FB0[idx]**2*PBDOT[idx])
        self.update(FB1new, name='FB1')

        # set the references
        if 'PBDOT_REF' in self.columns:
            FB1REFnew = np.full(self.catalogue_len, '', dtype='U32')
            PBDOTREF = self.catalogue['PBDOT_REF']
            FB1REFnew[idx] = PBDOTREF[idx]
            self.update(FB1REFnew, name='FB1_REF')

        # set the errors
        reqpars = ['PBDOT_ERR', 'PB_ERR']
        if np.all([p in self.columns for p in reqpars]):
            FB1ERRnew = np.full(self.catalogue_len, np.nan)
            PBDOTERR = self.catalogue['PBDOT_ERR']
            PBERR = self.catalogue['PB_ERR']
            idx = idx & np.isfinite(PBERR) & np.isfinite(PBDOTERR)
            FB1ERRnew[idx] = np.sqrt(
                (FB0[idx]**2 * PBDOTERR[idx])**2
                + (2.0 * FB0[idx]**3 * PBDOT[idx] *
                   PBERR[idx] * 86400.)**2)
            self.update(FB1ERRnew, name='FB1_ERR')

    def derived_p1_i(self):
        """
        Calculate the intrinsic period derivative.
        """

        if 'VTRANS' not in self.columns:
            self.derived_vtrans()

        reqpars = ['VTRANS', 'P0', 'P1', 'DIST']
        if not np.all([p in self.columns for p in reqpars]):
            return

        # get required parameters
        VTRANS = self.catalogue['VTRANS']
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']
        DIST = self.catalogue['DIST']

        P1I = np.full(self.catalogue_len, np.nan)
        idx = (np.isfinite(P1) & np.isfinite(P0) & np.isfinite(VTRANS) &
               np.isfinite(DIST))
        P1I[idx] = ((P1[idx]/1.0e-15) -
                    VTRANS[idx]**2 * 1.0e10 * P0[idx] /
                    (DIST[idx] * 3.086e6)/2.9979e10) * 1.0e-15
        self.update(P1I, name='P1_I')

    def derived_age(self):
        """
        Calculate the characteristic age in years (see
        :func:`~psrqpy.utils.characteristic_age`, with an assumed braking index
        of n=3).
        """

        from .utils import characteristic_age

        if not np.all([p in self.columns for p in ['P0', 'P1']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']
        AGE = characteristic_age(P0, P1)
        self.update(AGE, name='AGE')

    def derived_age_i(self):
        """
        Calculate the characteristic age (in years), dervied from period and
        intrinsic period derivative.
        """

        from .utils import characteristic_age

        if 'P1_I' not in self.columns:
            self.derived_p1_i()

        if not np.all([p in self.columns for p in ['P0', 'P1_I']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1_I = self.catalogue['P1_I']
        AGEI = characteristic_age(P0, P1_I)
        self.update(AGEI, name='AGE_I')

    def derived_bsurf(self):
        """
        Calculate the surface magnetic field strength (see
        :func:`~psrqpy.utils.B_field`).
        """

        from .utils import B_field

        if not np.all([p in self.columns for p in ['P0', 'P1']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']
        BSURF = B_field(P0, P1)
        self.update(BSURF, name='BSURF')

    def derived_bsurf_i(self):
        """
        Calculate the surface magnetic field strength, dervied from period and
        intrinsic period derivative.
        """

        from .utils import B_field

        if 'P1_I' not in self.columns:
            self.derived_p1_i()

        if not np.all([p in self.columns for p in ['P0', 'P1_I']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1_I = self.catalogue['P1_I']
        BSURFI = B_field(P0, P1_I)
        self.update(BSURFI, name='BSURF_I')

    def derived_gw_h0_spindown_limit(self):
        """
        Calculate the limit on the gravitational-wave emission amplitude at
        Earth assuming all rotational kinetic energy is lost via
        gravitational-waves generated from an l=m=2 mass quadrupole (i.e., a
        braking index of n=5). This uses the intrinsic spin-down values rather
        than the observed spin-downs.
        """

        from .utils import gw_h0_spindown_limit, pdot_to_fdot

        if 'P1_I' not in self.columns:
            self.derived_p1_i()

        if not np.all([p in self.columns for p in ["F0", "P1_I", "P1", "DIST"]]):
            return

        F0 = self.catalogue["F0"]
        P1I = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(self.catalogue["P1_I"])
        P1I[idx] = self.catalogue["P1_I"][idx]

        # where P1_I is not present use P1
        idx = ~idx & np.isfinite(self.catalogue["P1"])
        P1I[idx] = self.catalogue["P1"][idx]

        F1_I = pdot_to_fdot(P1I, frequency=F0)
        DIST = self.catalogue["DIST"]

        H0UL = gw_h0_spindown_limit(frequency=F0, fdot=F1_I, distance=DIST)
        self.update(H0UL, name="H0_SD")

    def derived_b_lc(self):
        """
        Calculate the magnetic field strength at the light cylinder.
        """

        if not np.all([p in self.columns for p in ['P0', 'P1']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']

        BLC = np.full(self.catalogue_len, np.nan)
        idx = (P1 > 0.) & np.isfinite(P1) & np.isfinite(P0)
        BLC[idx] = 3.0e8*np.sqrt(P1[idx])*np.abs(P0[idx])**(-5./2.)
        self.update(BLC, name='B_LC')

    def derived_edot(self):
        """
        Calculate the spin-down luminosity.
        """

        if not np.all([p in self.columns for p in ['P0', 'P1']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']

        EDOT = np.full(self.catalogue_len, np.nan)
        idx = (P1 > 0.) & np.isfinite(P1) & np.isfinite(P0)
        EDOT[idx] = 4.0 * np.pi**2 * 1e45 * P1[idx] / P0[idx]**3
        self.update(EDOT, name='EDOT')

    def derived_edot_i(self):
        """
        Calculate the spin-down luminosity, dervied from period and intrinsic
        period derivative.
        """

        if 'P1_I' not in self.columns:
            self.derived_p1_i()

        if not np.all([p in self.columns for p in ['P0', 'P1_I']]):
            return

        # get period and period derivative
        P0 = self.catalogue['P0']
        P1_I = self.catalogue['P1_I']

        EDOT_I = np.full(self.catalogue_len, np.nan)
        idx = (P1_I > 0.) & np.isfinite(P1_I) & np.isfinite(P0)
        EDOT_I[idx] = 4.0 * np.pi**2 * 1e45 * P1_I[idx] / P0[idx]**3
        self.update(EDOT_I, name='EDOT_I')

    def derived_edotd2(self):
        """
        Calculate the spin-down luminosity flux at the Sun.
        """

        reqpars = ['P0', 'P1', 'DIST']
        if not np.all([p in self.columns for p in reqpars]):
            return

        # get period, period derivative and distance
        P0 = self.catalogue['P0']
        P1 = self.catalogue['P1']
        DIST = self.catalogue['DIST']

        EDOTD2 = np.full(self.catalogue_len, np.nan)
        idx = (P0 > 0.) & np.isfinite(P1) & np.isfinite(P0) & np.isfinite(DIST)
        EDOTD2[idx] = 4.0*np.pi**2*1e45*((P1[idx]/P0[idx]**3)/DIST[idx]**2)
        self.update(EDOTD2, name='EDOTD2')

    def derived_pmtot(self):
        """
        Calculate the total proper motion and error.
        """

        reqpars = ['PMRA', 'PMDEC', 'PMELONG', 'PMELAT']
        if not np.all([p in self.columns for p in reqpars]):
            return

        # get PMRA and PMDEC
        PMRA = self.catalogue['PMRA'].copy()
        PMDEC = self.catalogue['PMDEC'].copy()
        PMELONG = self.catalogue['PMELONG']
        PMELAT = self.catalogue['PMELAT']

        # use PM ELONG or ELAT if no RA and DEC
        useelong = ~np.isfinite(PMRA) & np.isfinite(PMELONG)
        useelat = ~np.isfinite(PMDEC) & np.isfinite(PMELAT)
        PMRA[useelong] = PMELONG[useelong]
        PMDEC[useelat] = PMELAT[useelat]

        PMTOT = np.sqrt(PMRA**2+PMDEC**2)
        self.update(PMTOT, name='PMTOT')

        # get the error
        reqpars1 = ['PMRA_ERR', 'PMDEC_ERR']
        reqpars2 = ['PMELONG_ERR', 'PMELAT_ERR']
        if (not np.all([p in self.columns for p in reqpars1]) and
                not np.all([p in self.columns for p in reqpars2])):
            return

        if 'PMRA_ERR' in self.columns:
            PMRA_ERR = self.catalogue['PMRA_ERR'].copy()
        else:
            PMRA_ERR = np.full(self.catalogue_len, np.nan)

        if 'PMDEC_ERR' in self.columns:
            PMDEC_ERR = self.catalogue['PMDEC_ERR'].copy()
        else:
            PMDEC_ERR = np.full(self.catalogue_len, np.nan)

        if 'PMELONG_ERR' in self.columns:
            PMELONG_ERR = self.catalogue['PMELONG_ERR'].copy()
        else:
            PMELONG_ERR = np.full(self.catalogue_len, np.nan)

        if 'PMELAT_ERR' in self.columns:
            PMELAT_ERR = self.catalogue['PMELAT_ERR'].copy()
        else:
            PMELAT_ERR = np.full(self.catalogue_len, np.nan)

        PMLAT = np.full(self.catalogue_len, np.nan)
        PMLONG = np.full(self.catalogue_len, np.nan)
        PMLATERR = np.full(self.catalogue_len, np.nan)
        PMLONGERR = np.full(self.catalogue_len, np.nan)

        idx = np.isfinite(PMELONG) & np.isfinite(PMELONG_ERR)
        PMLONG[idx] = PMELONG[idx]
        PMLONGERR[idx] = PMELONG_ERR[idx]
        idx = np.isfinite(PMRA) & np.isfinite(PMRA_ERR)
        PMLONG[idx] = PMRA[idx]
        PMLONGERR[idx] = PMRA_ERR[idx]

        idx = np.isfinite(PMELAT) & np.isfinite(PMELAT_ERR)
        PMLAT[idx] = PMELAT[idx]
        PMLATERR[idx] = PMELAT_ERR[idx]
        idx = np.isfinite(PMDEC) & np.isfinite(PMDEC_ERR)
        PMLAT[idx] = PMDEC[idx]
        PMLATERR[idx] = PMDEC_ERR[idx]

        PMTOTERR = np.sqrt(((PMLONG*PMLONGERR)**2+(PMLAT*PMLATERR)**2) /
                           (PMLONG**2 + PMLAT**2))
        self.update(PMTOTERR, name='PMTOT_ERR')

    def derived_vtrans(self):
        """
        Calculate the transverse velocity.
        """

        if 'PMTOT' not in self.columns:
            self.derived_pmtot()

        if not np.all([p in self.columns for p in ['PMTOT', 'DIST']]):
            return

        PMTOT = self.catalogue['PMTOT']
        DIST = self.catalogue['DIST']

        VTRANS = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(PMTOT) & np.isfinite(DIST)
        VTRANS[idx] = (PMTOT[idx] * np.pi / (1000.0*3600.0*180.0*365.25 *
                       86400.0))*3.086e16*DIST[idx]
        self.update(VTRANS, name='VTRANS')

    def derived_flux(self):
        """
        Calculate spectral index between 400 and 1400 MHz and radio
        flux at 400 and 1400 MHz.
        """

        if not np.all([p in self.columns for p in ['S1400', 'S400']]):
            return

        S1400 = self.catalogue['S1400']
        S400 = self.catalogue['S400']

        SI414 = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(S1400) & np.isfinite(S400) & (S1400 > 0.) & (S400 > 0.)
        fac = np.log10(400.0/1400.0)
        SI414[idx] = -(np.log10(S400[idx]/S1400[idx])/fac)
        self.update(SI414, name='SI414')

        # need distance for flux
        if 'DIST' not in self.columns:
            self.define_dist()

        if 'DIST' not in self.columns:
            return

        DIST = self.catalogue['DIST']

        R_LUM = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(S400) & np.isfinite(DIST)
        R_LUM[idx] = S400[idx] * DIST[idx]**2
        self.update(R_LUM, name='R_LUM')

        R_LUM14 = np.full(self.catalogue_len, np.nan)
        idx = np.isfinite(S1400) & np.isfinite(DIST)
        R_LUM14[idx] = S1400[idx] * DIST[idx]**2
        self.update(R_LUM14, name='R_LUM14')

    def get_pulsar(self, psr, selected=False):
        """
        Return the table row for a particular pulsar for all the catalogue
        parameters.

        Args:
            psr (str): The name of a pulsar to return.
            selected (bool): If True then output return a table row containing
                parameters specified by :meth:`~psrqpy.QueryATNF.query_params`,
                otherwise return all parameters. Defaults to False.

        Returns:
            :class:`astropy.table.Table`: a table row
        """

        namepars = ['PSRJ', 'PSRB', 'BNAME', 'JNAME', 'NAME']
        if not np.any([p in self.columns for p in namepars]):
            warnings.warn("No 'NAME' parameter in table!")
            return None

        # try searching for the name in each potential name-type
        for namepar in namepars:
            if namepar in self.columns:
                names = self.catalogue[namepar]
                if np.any(psr == names):
                    psrrow = self.catalogue_table[(psr == names).tolist()]
                    if selected:
                        return psrrow[self.query_params]
                    else:
                        return psrrow

        return None

    def get_ephemeris(self, psr, precision=15, selected=False):
        """
        Return the table row for a particular pulsar and output it as an
        ephemeris-style string output.

        Args:
            psr (str): The name of a pulsar to return.
            precision (int): The precision (number of decimal places) at which
                to output numbers. Defaults to 15.
            selected (bool): If True only output the parameters specified by
                :meth:`~psrqpy.QueryATNF.query_params`, otherwise output all
                parameters. Defaults to False.

        Returns:
            str: an ephemeris
        """

        psr = self.get_pulsar(psr, selected=selected)

        if psr is None:
            return None

        ephemstr = ''

        # get variables that are set
        variables = []
        values = []
        errors = []
        for par in PSR_ALL_PARS:
            if par in psr.columns:
                parval = psr[par]

                if type(parval) == MaskedColumn:
                    # required for astropy versions below v4
                    if not parval.mask[0]:
                        variables.append(par)
                        values.append(parval[0])

                        if par+'_ERR' in psr.columns:
                            errval = psr[par+'_ERR']

                            if type(errval) == MaskedColumn:
                                if not errval.mask[0]:
                                    errors.append(errval[0])
                                else:
                                    errors.append(None)
                            else:
                                errors.append(errval[0])
                        else:
                            errors.append(None)
                elif type(parval) == Column:
                    variables.append(par)
                    values.append(parval[0])

                    if par+'_ERR' in psr.columns:
                        errval = psr[par+'_ERR']

                        if type(errval) == MaskedColumn:
                            if not errval.mask[0]:
                                errors.append(errval[0])
                            else:
                                errors.append(None)
                        else:
                            errors.append(errval[0])
                    else:
                        errors.append(None)
                else:
                    raise TypeError("Pulsar contains a non-column type!")

        mkl = max([len(kn) for kn in variables])+2  # max key length for output alignment
        vlb = precision + 10  # allow extra space for minus sign/exponents
        outputstr = '{{name: <{0}}}{{value: <{1}}}\t{{error}}'.format(mkl, vlb)

        for varname, varval, varerr in zip(variables, values, errors):
            outputdic = {}
            outputdic['name'] = varname

            if isinstance(varval, float):
                if varval.is_integer():
                    precstr = '{0:.0f}'  # print out an integer
                else:
                    precstr = '{{0:.{}f}}'.format(precision)  # print out float

                if abs(varval) < 1e-6 or abs(varval) > 1e6:
                    # print out float in scientific notation
                    precstr = '{{0:.{}e}}'.format(precision)

                if varerr is not None:
                    if varerr.is_integer():
                        precstre = '{0:.0f}'  # print out an integer
                else:
                    precstre = '{{0:.{}f}}'.format(precision)  # print out float

                if varerr is not None:
                    if abs(varerr) < 1e-6 or abs(varerr) > 1e6:
                        # print out float in scientific notation
                        precstre = '{{0:.{}e}}'.format(precision)

                outputdic['value'] = precstr.format(varval)
                outputdic['error'] = precstre.format(varerr) if varerr is not None else ''
            else:
                outputdic['value'] = varval

                if isinstance(varerr, float):
                    precstre = '{{0:.{}f}}'.format(precision)  # print out float
                    if abs(varerr) < 1e-6 or abs(varerr) > 1e6:
                        # print out float in scientific notation
                        precstre = '{{0:.{}e}}'.format(precision)

                    outputdic['error'] = precstre.format(varerr)
                else:
                    outputdic['error'] = ''

            ephemstr += outputstr.format(**outputdic).strip()+'\n'

        return ephemstr

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
            qparams = deepcopy(self.query_params)
            if isinstance(qparams, list):
                if 'JNAME' not in qparams:
                    self.query_params = self.query_params + ['JNAME']
            else:
                self.query_params = ['JNAME']

            psrtable = self.table
            for row in psrtable:
                P = Pulsar(row['JNAME'], query=self)
                self._pulsars.add_pulsar(P)

            self.query_params = qparams  # revert to previous query parameters

        return self._pulsars

    @property
    def get_version(self):
        """
        Return a string with the ATNF version number, or None if not found.

        Returns:
            str: the ATNF version number.
        """

        return self.catalogue.version

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
                    if not isinstance(p, str):
                        raise Exception("Non-string value '{}' found in pulsar type list"
                                        .format(p))
                self._query_psr_types = psrtype
            else:
                if isinstance(psrtype, str):
                    self._query_psr_types = [psrtype]
                else:
                    raise Exception("'psrtype' must be a list or string")

            for p in list(self._query_psr_types):
                if p.upper() not in PSR_TYPE:
                    warnings.warn("Pulsar type '{}' is not recognised, no type will be required"
                                  .format(p))
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
                    if not isinstance(p, str):
                        raise Exception("Non-string value '{}' found in associations list"
                                        .format(p))
                self._query_assocs = assoc
            else:
                if isinstance(assoc, str):
                    self._query_assocs = [assoc]
                else:
                    raise Exception("'assoc' must be a list or string")

            for p in list(self._query_assocs):
                if p.upper() not in PSR_ASSOC_TYPE:
                    warnings.warn("Pulsar association '{}' is not recognised, "
                                  "no type will be required".format(p))
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
                    if not isinstance(p, str):
                        raise Exception("Non-string value '{}' found in binary "
                                        "companions list".format(p))
                self._query_bincomps = bincomp
            else:
                if isinstance(bincomp, str):
                    self._query_bincomps = [bincomp]
                else:
                    raise Exception("'bincomp' must be a list or string")

            for p in list(self._query_bincomps):
                if p.upper() not in PSR_BINARY_TYPE:
                    warnings.warn("Pulsar binary companion '{}' is not recognised, "
                                  "no type will be required".format(p))
                    self._query_bincomps.remove(p)
                else:
                    if len(conditionparse) == 0:
                        conditionparse = 'bincomp({})'.format(p.upper())
                    else:
                        conditionparse += ' && bincomp({})'.format(p.upper())

        return conditionparse

    @property
    def catalogue_shape(self):
        """
        The shape of the entire catalogue table as a tuple containing the
        number of rows and the number of columns.
        """

        return self.catalogue.shape

    @property
    def catalogue_nrows(self):
        """
        The number of rows in the entire catalogue, i.e. the number of pulsars
        it contains.
        """

        return self.catalogue.shape[0]

    @property
    def catalogue_ncols(self):
        """
        The number of columns in the entire catalogue, i.e. the number of
        parameters it contains.
        """

        return self.catalogue.shape[1]

    @property
    def catalogue_len(self):
        """
        The length of the entire catalogue, i.e., the number of pulsars it
        contains. This should be the same as ``catalogue_nrows``.
        """

        return len(self.catalogue)

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
              periodlims=None, usecondition=True, usepsrs=True, rcparams={}):
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
                pulsars that conform to the original query condition values.
                Defaults to True.
            usepsrs (bool): if True create the P-Pdot diagram only with pulsars
                specified in the original query. Defaults to True.
            rcparams (dict): a dictionary of :py:obj:`matplotlib.rcParams`
                setup parameters for the plot.

        Returns:
            :class:`matplotlib.figure.Figure`: the figure object
        """

        try:
            import matplotlib as mpl
            from matplotlib import pyplot as pl
        except ImportError:
            raise ImportError('Cannot produce P-Pdot plot as Matplotlib is '
                              'not available')

        from .utils import death_line, label_line

        # get table containing all required parameters
        table = self.query_table(usecondition=usecondition,
                                 usepsrs=usepsrs,
                                 query_params=['P0', 'P1', 'P1_I', 'ASSOC',
                                               'BINARY', 'TYPE'])

        if len(table) == 0:
            print("No pulsars found, so no P-Pdot plot has been produced")
            return None

        if isinstance(showtypes, str):
            nshowtypes = [showtypes]
        else:
            nshowtypes = showtypes

        for stype in list(nshowtypes):
            if 'ALL' == stype.upper():
                nshowtypes = list(PSR_TYPE)
                # remove radio as none are returned as this
                del nshowtypes[nshowtypes.index('RADIO')]
                break
            elif stype.upper() not in list(PSR_TYPE):
                warnings.warn('"TYPE" {} is not recognised, so will not be '
                              'included'.format(stype))
                del nshowtypes[nshowtypes.index(stype)]
            if 'SGR' == stype.upper():  # synonym for AXP
                nshowtypes[nshowtypes.index(stype)] = 'AXP'

        # set plot parameters
        rcparams['figure.figsize'] = rcparams['figure.figsize'] if \
            'figure.figsize' in rcparams else (9, 9.5)
        rcparams['figure.dpi'] = rcparams['figure.dpi'] if \
            'figure.dpi' in rcparams else 250
        rcparams['text.usetex'] = rcparams['text.usetex'] if \
            'text.usetex' in rcparams else True
        rcparams['axes.linewidth'] = rcparams['axes.linewidth'] if \
            'axes.linewidth' in rcparams else 0.5
        rcparams['axes.grid'] = rcparams['axes.grid'] if \
            'axes.grid' in rcparams else False
        rcparams['font.family'] = rcparams['font.family'] if \
            'font.family' in rcparams else 'sans-serif'
        rcparams['font.sans-serif'] = rcparams['font.sans-serif'] if \
            'font.sans-serif' in rcparams else \
            'Avant Garde, Helvetica, Computer Modern Sans serif'
        rcparams['font.size'] = rcparams['font.size'] if \
            'font.size' in rcparams else 20
        rcparams['legend.fontsize'] = rcparams['legend.fontsize'] if \
            'legend.fontsize' in rcparams else 16
        rcparams['legend.frameon'] = rcparams['legend.frameon'] if \
            'legend.frameon' in rcparams else False

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

        if 'ASSOC' in table.columns:
            assocs = table['ASSOC'][pidx]     # associations
        if 'TYPE' in table.columns:
            types = table['TYPE'][pidx]       # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = table['BINARY'][pidx]  # binary pulsars

        # now get only positive pdot values
        pidx = pdots > 0.
        periods = periods[pidx]
        pdots = pdots[pidx]
        if 'ASSOC' in table.columns:
            assocs = assocs[pidx]      # associations
        if 'TYPE' in table.columns:
            types = types[pidx]        # pulsar types
        if 'BINARY' in nshowtypes:
            binaries = binaries[pidx]  # binary pulsars

        # check whether to exclude globular cluster pulsars that could have
        # contaminated spin-down value
        if excludeGCs:
            # use '!=' to find GC indexes
            nongcidxs = np.flatnonzero(
                np.char.find(np.array(assocs.tolist(),
                                      dtype=str), 'GC:') == -1)
            periods = periods[nongcidxs]
            pdots = pdots[nongcidxs]
            if 'ASSOC' in table.columns:
                assocs = assocs[nongcidxs]
            if 'TYPE' in table.columns:
                types = types[nongcidxs]
            if 'BINARY' in nshowtypes:
                binaries = binaries[nongcidxs]

        # plot pulsars
        ax.loglog(periods, pdots, marker='.', color='dimgrey',
                  linestyle='none')
        ax.set_xlabel(r'Period (s)')
        ax.set_ylabel(r'Period Derivative')

        # get limits
        if periodlims is None:
            periodlims = [10**np.floor(np.min(np.log10(periods))),
                          10.*int(np.ceil(np.max(pdots)/10.))]
        if pdotlims is None:
            pdotlims = [10**np.floor(np.min(np.log10(pdots))),
                        10**np.ceil(np.max(np.log10(pdots)))]
        ax.set_xlim(periodlims)
        ax.set_ylim(pdotlims)

        if deathline:
            deathpdots = 10**death_line(np.log10(periodlims),
                                        linemodel=deathmodel)
            ax.loglog(periodlims, deathpdots, 'k--', linewidth=0.5)

            if filldeath:
                if not filldeathtype:
                    filldeathtype = {}

                filldeathtype['linestyle'] = filldeathtype['linestyle'] if \
                    'linestyle' in filldeathtype else '-'
                filldeathtype['alpha'] = filldeathtype['alpha'] if \
                    'alpha' in filldeathtype else 0.15
                filldeathtype['facecolor'] = filldeathtype['facecolor'] if \
                    'facecolor' in filldeathtype else 'darkorange'
                filldeathtype['hatch'] = filldeathtype['hatch'] if \
                    'hatch' in filldeathtype else ''
                ax.fill_between(periodlims, deathpdots, pdotlims[0],
                                **filldeathtype)

        # add markers for each pulsar type
        if not markertypes:
            markertypes = {}

        # check if markers have been defined by the user or not
        markertypes['AXP'] = {'marker': 's', 'markeredgecolor': 'red'} if \
            'AXP' not in markertypes else markertypes['AXP']
        markertypes['BINARY'] = {'marker': 'o', 'markeredgecolor': 'grey'} if \
            'BINARY' not in markertypes else markertypes['BINARY']
        markertypes['HE'] = {'marker': 'D', 'markeredgecolor': 'orange'} if \
            'HE' not in markertypes else markertypes['HE']
        markertypes['RRAT'] = {'marker': 'h', 'markeredgecolor': 'green'} if \
            'RRAT' not in markertypes else markertypes['RRAT']
        markertypes['NRAD'] = {'marker': 'v', 'markeredgecolor': 'blue'} if \
            'NRAD' not in markertypes else markertypes['NRAD']
        markertypes['XINS'] = {'marker': '^', 'markeredgecolor': 'magenta'} if \
            'XINS' not in markertypes else markertypes['XINS']
        markertypes['GC'] = {'marker': '8', 'markeredgecolor': 'cyan'} if \
            'GC' not in markertypes else markertypes['GC']
        markertypes['SNR'] = {'marker': '*', 'markeredgecolor': 'darkorchid'} if \
            'SNR' not in markertypes else markertypes['SNR']

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

        handles = dict()

        for stype in nshowtypes:
            if stype.upper() in PSR_TYPE + ['GC', 'SNR']:
                thistype = stype.upper()
                if thistype == 'BINARY':
                    # for binaries used the 'BINARY' column in the table
                    if type(binaries) == MaskedColumn:
                        typeidx = np.flatnonzero(~binaries.mask)
                    else:
                        typeidx = np.ones(len(binaries))
                elif thistype in ['GC', 'SNR']:
                    typeidx = np.flatnonzero(
                        np.char.find(np.array(assocs.tolist(),
                                              dtype=str), thistype) != -1)
                else:
                    typeidx = np.flatnonzero(
                        np.char.find(np.array(types.tolist(),
                                              dtype=str), thistype) != -1)

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

                ax.legend(handles.values(), handles.keys(), loc='upper left',
                          numpoints=1)

        # add characteristic age lines
        tlines = dict()
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
                    tlines[r'${{{0:.1f}}}!\times\!10^{{{1:d}}}\,{{\rm yr}}$'
                           .format(numv, taupow)] = tline

        # add magnetic field lines
        Blines = dict()
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
                    Blines[r'${{{0:.1f}}}!\times\!10^{{{1:d}}}\,{{\rm G}}$'
                           .format(numv, Bpow)] = bline

        fig.tight_layout()

        # add text for characteristic age lines and magnetic field strength lines
        for l in tlines:
            _ = label_line(ax, tlines[l], l, color='k', fs=18, frachoffset=0.05)

        for l in Blines:
            _ = label_line(ax, Blines[l], l, color='k', fs=18, frachoffset=0.90)

        # return the figure
        return fig

    def gw_mass_quadrupole(self):
        """
        Return the :math:`l=m=2` mass quadrupoles based on the spin-down limits
        for the pulsars in the catalogue using :func:`~psrqpy.utils.h0_to_q22`.

        Returns:
            :class:`astropy.table.Table`: a table containing ``PSRJ`` names and
                :math:`Q_{22}` values as a column called ``Q22``.
        """

        # parameters required to calculate the mass quadrupole
        requiredpars = ["PSRJ", "H0_SD", "F0", "DIST"]

        table = self.query_table(query_params=requiredpars)

        q22 = h0_to_q22(
            table["H0_SD"],
            table["F0"],
            table["DIST"]
        )

        idx = np.isfinite(q22)
        return Table(
            data=[table["PSRJ"][idx], q22[idx]],
            names=["PSRJ", "Q22"],
            units=[None, aunits.kg * aunits.m ** 2],
        )

    def gw_ellipticity(self):
        """
        Return the :math:`l=m=2` mass quadrupoles based on the spin-down limits
        for the pulsars in the catalogue using
        :func:`~psrqpy.utils.h0_to_ellipticity`.

        Returns:
            :class:`astropy.table.Table`: a table containing ``PSRJ`` names and
                :math:`\\varepsilon` values as a column called ``ELL``.
        """

        q22 = self.gw_mass_quadrupole()
        ell = q22_to_ellipticity(q22["Q22"])

        return Table(
            data=[q22["PSRJ"], ell],
            names=["PSRJ", "ELL"],
        )
