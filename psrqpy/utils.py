"""
A selection of useful functions used by the module.
"""

from __future__ import print_function, division

import warnings
import re
import os
import numpy as np
import requests
import tarfile
from bs4 import BeautifulSoup

from six import string_types

from collections import OrderedDict

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as aunits
from astropy.utils.data import download_file, clear_download_cache
from pandas import DataFrame

from .config import (ATNF_BASE_URL, ADS_URL, ATNF_TARBALL,
                     PSR_ALL, PSR_ALL_PARS, GLITCH_URL)


# problematic references that are hard to parse
PROB_REFS = ['bwck08', 'crf+18']


def get_catalogue(path_to_db=None, cache=True, update=False, pandas=False):
    """
    This function will attempt to download and cache the entire ATNF Pulsar
    Catalogue database `tarball
    <http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz>`_,
    or read in database file from a provided path. The database will be
    converted into an :class:`astropy.table.Table` or
    :class:`pandas.DataFrame`. This was originally based on the method in the
    `ATNF.ipynb
    <https://github.com/astrophysically/ATNF-Pulsar-Cat/blob/master/ATNF.ipynb>`_
    notebook by Joshua Tan (`@astrophysically
    <https://github.com/astrophysically/>`_).

    Args:
        path_to_db (str): if the path to a local version of the database file
            is given then that will be read in rather than attempting to
            download the file (defaults to None).
        cache (bool): cache the downloaded ATNF Pulsar Catalogue file. Defaults
            to True. This is ignored if `path_to_db` is given.
        update (bool): if True the ATNF Pulsar Catalogue will be
            re-downloaded and cached if there has been a change compared to the
            currently cached version. This is ignored if `path_to_db` is given.
        pandas (bool): if True the catalogue will be returned as a
            :class:`pandas.DataFrame` rather than the default of an
            :class:`~astropy.table.Table`.

    Returns:
        :class:`~astropy.table.Table` or :class:`~pandas.DataFrame`: a table
        containing the entire catalogue.

    """

    if path_to_db is None:
        # remove any cached file if requested
        if update:
            if check_update():
                clear_download_cache(ATNF_TARBALL)

        # get the tarball
        try:
            dbtarfile = download_file(ATNF_TARBALL, cache=cache)
        except IOError:
            raise IOError('Problem accessing ATNF catalogue tarball')

        try:
            # open tarball
            pulsargz = tarfile.open(dbtarfile, mode='r:gz')

            # extract the database file
            dbfile = pulsargz.extractfile('psrcat_tar/psrcat.db')
        except IOError:
            raise IOError('Problem extracting the database file')
    else:
        try:
            dbfile = open(path_to_db, 'r')
        except IOError:
            raise IOError('Error loading given database file')

    breakstring = '@'    # break between each pulsar
    commentstring = '#'  # specifies line is a comment

    # create list of dictionaries - one for each pulsar
    psrlist = [{}]

    version = None  # catalogue version

    # loop through lines in dbfile
    for line in dbfile.readlines():
        if isinstance(line, string_types):
            dataline = line.split()
        else:
            dataline = line.decode().split()   # Splits on whitespace

        if dataline[0][0] == commentstring:
            # get catalogue version (should be in first comment string)
            if dataline[0] == '#CATALOGUE' and len(dataline) == 2:
                version = dataline[1]
            continue

        if dataline[0][0] == breakstring:
            # First break comes at the end of the first object and so forth
            psrlist.append({})  # New object!
            continue

        try:
            psrlist[-1][dataline[0]] = float(dataline[1])
        except ValueError:
            psrlist[-1][dataline[0]] = dataline[1]

        if len(dataline) > 2:
            # check whether 3rd value is a float (so its an error value) or not
            try:
                float(dataline[2])
                isfloat = True
            except ValueError:
                isfloat = False

            if isfloat:
                # error values are last digit errors, so convert to actual
                # errors by finding the number of decimal places after the
                # '.' in the value string
                val = dataline[1].split(':')[-1]  # account for RA and DEC strings

                try:
                    float(val)
                except ValueError:
                    raise ValueError("Value with error is not convertable to a float")

                if dataline[2][0] == '-' or '.' in dataline[2]:
                    # negative errors or those with decimal points are absolute values
                    scalefac = 1.
                else:
                    # split on exponent
                    valsplit = re.split('e|E|d|D', val)
                    scalefac = 1.
                    if len(valsplit) == 2:
                        scalefac = 10**(-int(valsplit[1]))

                    dpidx = valsplit[0].find('.')  # find position of decimal point
                    if dpidx != -1:  # a point is found
                        scalefac *= 10**(len(valsplit[0])-dpidx-1)

                # add error column if required
                psrlist[-1][dataline[0]+'_ERR'] = float(dataline[2])/scalefac  # error entry
            else:
                # add reference column if required
                psrlist[-1][dataline[0]+'_REF'] = dataline[2]  # reference entry

            if len(dataline) > 3:
                # last entry must(!) be a reference
                psrlist[-1][dataline[0]+'_REF'] = dataline[3]  # reference entry

    dbfile.close()   # close tar file
    if not path_to_db:
        pulsargz.close()

    del psrlist[-1]  # Final breakstring comes at the end of the file

    # add RA and DEC in degs and JNAME/BNAME
    for i, psr in enumerate(list(psrlist)):
        if 'RAJ' in psr.keys() and 'DECJ' in psr.keys():
            # check if the string can be converted to a float (there are a few
            # cases where the position is just a decimal value)
            try:
                rad = float(psr['RAJ'])
                ras = Angle(rad*aunits.hourangle)
                psr['RAJ'] = ras.to_string(sep=':', pad=True)
            except ValueError:
                pass

            try:
                decd = float(psr['DECJ'])
                decs = Angle(decd*aunits.deg)
                psr['DECJ'] = decs.to_string(sep=':', pad=True, alwayssign=True)
            except ValueError:
                pass

            coord = SkyCoord(psr['RAJ'], psr['DECJ'],
                             unit=(aunits.hourangle, aunits.deg))
            psrlist[i]['RAJD'] = coord.ra.deg    # right ascension in degrees
            psrlist[i]['DECJD'] = coord.dec.deg  # declination in degrees

        # add 'JNAME', 'BNAME' and 'NAME'
        if 'PSRJ' in psr.keys():
            psrlist[i]['JNAME'] = psr['PSRJ']
            psrlist[i]['NAME'] = psr['PSRJ']
            if 'PSRJ_REF' in psr.keys():
                psrlist[i]['JNAME_REF'] = psr['PSRJ_REF']
                psrlist[i]['NAME_REF'] = psr['PSRJ_REF']

        if 'PSRB' in psr.keys():
            psrlist[i]['BNAME'] = psr['PSRB']
            if 'PSRB_REF' in psr.keys():
                psrlist[i]['BNAME_REF'] = psr['PSRB_REF']

            if 'NAME' not in psrlist[i].keys():
                psrlist[i]['NAME'] = psr['PSRB']
                if 'PSRB_REF' in psr.keys():
                    psrlist[i]['NAME_REF'] = psr['PSRB_REF']

    # convert to a pandas DataFrame - this will fill in empty spaces
    dftable = DataFrame(psrlist)

    if pandas:
        # return pandas DataFrame
        dftable.version = version

        return dftable

    # convert into an astropy table
    psrtable = Table.from_pandas(dftable)

    # add units if known
    for key in PSR_ALL_PARS:
        if key in psrtable.colnames:
            if PSR_ALL[key]['units']:
                psrtable.columns[key].unit = PSR_ALL[key]['units']

                if PSR_ALL[key]['err'] and key+'_ERR' in psrtable.colnames:
                    psrtable.columns[key+'_ERR'].unit = PSR_ALL[key]['units']

    # add metadata
    if not path_to_db:
        if version is not None:
            psrtable.meta['version'] = version
        else:
            psrtable.meta['version'] = None
            warnings.warn('No version number found in the database file',
                          UserWarning)
        psrtable.meta['ATNF Pulsar Catalogue'] = ATNF_BASE_URL

    if path_to_db:
        psrtable.meta['Database file'] = path_to_db

    return psrtable


def check_update():
    """
    Check if the ATNF Pulsar Catalogue has been updated compared to the version
    in the cache.

    Returns:
       bool: True if the cache can be updated.

    """

    from astropy.utils.data import download_file, get_cached_urls, compute_hash

    if ATNF_TARBALL not in get_cached_urls():
        # can update cache as file is not cached yet
        return True

    # get the cached file name
    cachefile = download_file(ATNF_TARBALL, cache=True)

    # download a new version of the file and check the hash
    tmpcache = download_file(ATNF_TARBALL, cache=False)

    curhash = compute_hash(cachefile)
    tmphash = compute_hash(tmpcache)

    if curhash == tmphash:
        # no update needed
        return False
    else:
        # an update can be obtained
        return True


def get_glitch_catalogue(psr=None):
    """
    Return a :class:`~astropy.table.Table` containing the `Jodrell Bank pulsar
    glitch catalogue <http://www.jb.man.ac.uk/pulsar/glitches/gTable.html>`_.
    If using data from the glitch catalogue then please cite `Espinoza et al.
    (2011) <http://adsabs.harvard.edu/abs/2011MNRAS.414.1679E>`_ and the URL
    `<http://www.jb.man.ac.uk/pulsar/glitches.html>`_.

    The output table will contain the following columns:

     * `NAME`: the pulsars common name
     * `JNAME`: the pulsar name based on J2000 coordinates
     * `Glitch number`: the number of the glitch for a particular pulsar in chronological order
     * `MJD`: the time of the glitch in Modified Julian Days
     * `MJD_ERR`: the uncertainty on the glitch time in days
     * `DeltaF/F`: the fractional frequency change
     * `DeltaF/F_ERR`: the uncertainty on the fractional frequency change
     * `DeltaF1/F1`: the fractional frequency derivative change
     * `DeltaF1/F1_ERR`: the uncertainty on the fractional frequency derivative change
     * `Reference`: the glitch publication reference

    Args:
        psr (str): if a pulsar name is given then only the glitches for that
            pulsar are returned, otherwise all glitches are returned.

    Returns:
        :class:`~astropy.table.Table`: a table containing the entire glitch
        catalogue.

    Example:
        An example of using this to extract the glitches for the Crab Pulsar
        would be:

        >>> import psrqpy
        >>> gtable = psrqpy.get_glitch_catalogue(psr='J0534+2200')
        >>> print("There have been {} glitches observed from the Crab pulsar".format(len(gtable)))
        27
    """

    try:
        from astropy.table import Table
        from astropy.units import Unit
    except ImportError:
        raise ImportError('Problem importing astropy')

    # get webpage
    try:
        gt = requests.get(GLITCH_URL)
    except Exception as e:
        raise RuntimeError("Error downloading glitch catalogue: {}".format(str(e)))

    if gt.status_code != 200:
        warnings.warn("Count not query the glitch catalogue.", UserWarning)
        return None

    # parse HTML
    try:
        soup = BeautifulSoup(gt.content, 'html.parser')
    except Exception as e:
        warnings.warn("Count not parse the glitch catalogue: {}".format(str(e)), UserWarning)
        return None

    # get table rows
    rows = soup.table.find_all('tr')

    # set the table headings
    tabledict = OrderedDict()
    tabledict['NAME'] = []
    tabledict['JNAME'] = []
    tabledict['Glitch number'] = []
    tabledict['MJD'] = []
    tabledict['MJD_ERR'] = []
    tabledict['DeltaF/F'] = []
    tabledict['DeltaF/F_ERR'] = []
    tabledict['DeltaF1/F1'] = []
    tabledict['DeltaF1/F1_ERR'] = []
    tabledict['Reference'] = []

    # loop through rows: rows with glitches have their first column as an index
    for row in rows:
        tds = row.find_all('td')

        if tds[0].contents[0].string is None:
            continue

        try:
            tabledict['NAME'].append(tds[1].contents[0].string)
            jname = 'J'+tds[2].contents[0].string if \
                    'J' != tds[2].contents[0].string[0] else tds[2].contents[0].string
            tabledict['JNAME'].append(jname)
            tabledict['Glitch number'].append(int(tds[3].contents[0].string))

            for j, pname in enumerate(['MJD', 'MJD_ERR', 'DeltaF/F',
                                       'DeltaF/F_ERR', 'DeltaF1/F1',
                                       'DeltaF1/F1_ERR']):
                try:
                    val = float(tds[4+j].contents[0].string)
                except ValueError:
                    val = np.nan

                tabledict[pname].append(val)

            # get reference link if present
            try:
                ref = tds[10].contents[0].a.attrs['href']
            except AttributeError:
                ref = tds[10].contents[0].string
            tabledict['Reference'].append(ref)
        except RuntimeError:
            warnings.warn("Problem parsing glitch table", UserWarning)
            return None

    # convert to an astropy table
    table = Table(tabledict)
    table.columns['MJD'].unit = Unit('d')     # add units of days to glitch time
    table.columns['MJD_ERR'].unit = Unit('d')

    # correct scaling of parameters
    table['DeltaF/F'] *= 1e-9
    table['DeltaF/F_ERR'] *= 1e-9
    table['DeltaF1/F1'] *= 1e-3
    table['DeltaF1/F1_ERR'] *= 1e-3

    if psr is None:
        return table
    else:
        if psr not in table['NAME'] and psr not in table['JNAME']:
            warnings.warn("Pulsar '{}' not found in glitch catalogue".format(psr), UserWarning)
            return None
        else:
            if psr in table['NAME']:
                return table[table['NAME'] == psr]
            else:
                return table[table['JNAME'] == psr]


def get_references(useads=False, cache=True):
    """
    Return a dictionary of paper
    `reference <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_ref.html>`_
    in the ATNF catalogue. The keys are the ref strings given in the ATNF
    catalogue.

    Args:
        useads (bool): boolean to set whether to use the python mod:`ads`
            module to get the NASA ADS URL for the references.
        cache (bool): use cached, or cache, the reference bundled with the
            catalogue tarball.

    Returns:
        dict: a dictionary of references.
    """

    import tempfile
    import json

    # get the tarball
    try:
        dbtarfile = download_file(ATNF_TARBALL, cache=cache)
    except IOError:
        raise IOError('Problem accessing ATNF catalogue tarball')

    try:
        # open tarball
        pulsargz = tarfile.open(dbtarfile, mode='r:gz')

        # extract the references
        reffile = pulsargz.extractfile('psrcat_tar/psrcat_ref')
    except IOError:
        raise IOError('Problem extracting the database file')

    refdic = {}
    refidx = 0
    thisref = ''

    for line in reffile.readlines():
        if isinstance(line, string_types):
            thisline = line
        else:
            thisline = line.decode()

        if thisline[0:3] == '***':
            if refidx > 0:
                # return reference making sure to only have single spaces
                refdic[thisname] = re.sub(r'\s+', ' ', thisref)
            thisref = ''
            refidx += 1
            thisname = thisline.split()[0].strip('***')
            thisref += thisline[thisline.find(':')+1:]
        else:
            # make sure there is a space so words don't get concatenated
            thisref += ' '
            thisref += thisline.strip()

    reffile.close()
    pulsargz.close()  # close tar file

    # if not requiring ADS references just return the current dictionary
    if not useads:
        return refdic
    else:
        try:
            import ads
            from ads.exceptions import APIResponseError
        except ImportError:
            warnings.warn('Could not import ADS module, so no ADS information '
                          'will be included', UserWarning)
            return refdic, None

    # try getting cached references
    if not cache:
        adsrefs = {}
    else:
        from astropy.utils.data import is_url_in_cache

        tmpdir = tempfile.gettempdir()  # get system "temporary" directory
        dummyurl = 'file://{}/ads_cache'.format(tmpdir)
        dummyfile = os.path.join('{}'.format(tmpdir), 'ads_cache')

        # check if cached ADS refs list exists (using dummy URL)
        if is_url_in_cache(dummyurl):
            adsfile = download_file(dummyurl, cache=True, show_progress=False)

            try:
                fp = open(adsfile, 'r')
            except IOError:
                warnings.warn('Could not load ADS URL cache for references',
                              UserWarning)
                return refdic, None

            adsrefs = json.load(fp)
            fp.close()

            return refdic, adsrefs
        else:
            adsrefs = {}

    # loop over references
    j = 0
    for reftag in refdic:
        j = j + 1

        if reftag in PROB_REFS:
            continue

        refstring = refdic[reftag]

        # try getting the year from the string and split on this (allows years
        # between 1000-2999 and followed by a lowercase letter, e.g. 2009 or
        # 2009a)
        match = re.match(r'.*([1-2][0-9]{3}[az]{1}|[1-2][0-9]{3})', refstring)
        if match is None:
            continue

        # do splitting
        spl = re.split(r'([1-2][0-9]{3}[az]{1}|[1-2][0-9]{3})', refstring)

        if len(spl) != 3:
            # more than 1 "year", so ignore!
            continue

        year = spl[1] if len(spl[1]) == 4 else spl[1][:4]

        try:
            int(year)
        except ValueError:
            # "year" is not an integer
            continue

        # get the authors (remove line breaks/extra spaces and final full-stop)
        authors = spl[0].strip().strip('.')

        # remove " Jr." from any author names (as it causes issues!)
        authors = authors.replace(' Jr.', '')

        # separate out authors
        sepauthors = authors.split('.,')[:-1]

        if len(sepauthors) == 0:
            # no authors were parsed
            continue

        # remove any "'s for umlauts in author names
        sepauthors = [a.replace(r'"', '') for a in sepauthors]

        # split any authors that are seperated by an ampersand
        if '&' in sepauthors[-1] or 'and' in sepauthors[-1]:
            lastauthors = [a.strip() for a in re.split(r'& | and ', sepauthors.pop(-1))]
            sepauthors = sepauthors + lastauthors
            for i in range(len(sepauthors)-2):
                sepauthors[i] += '.'  # re-add final full stops where needed
            sepauthors[-1] += '.'
        else:
            sepauthors = [a+'.' for a in sepauthors]  # re-add final full stops

        # get the title
        try:
            # remove preceding or trailing full stops
            title = spl[2].strip('.').split('.')[0].strip()
        except RuntimeError:
            # could not get title so ignore this entry
            continue

        # try getting ADS references
        try:
            article = ads.SearchQuery(year=year, first_author=sepauthors[0],
                                      title=title)
        except APIResponseError:
            warnings.warn('Could not get reference information, so no ADS '
                          'information will be included', UserWarning)
            continue

        try:
            adsrefs[reftag] = ADS_URL.format(list(article)[0].bibcode)
        except (IndexError, APIResponseError):
            pass

    if cache:
        # output adsrefs to cache file
        try:
            # output to dummy temporary file and then "download" to cache
            fp = open(dummyfile, 'w')
            json.dump(adsrefs, fp, indent=2)
            fp.close()
        except IOError:
            raise IOError("Could not output the ADS references to a file")

        # cache the file
        _ = download_file(dummyurl, cache=True, show_progress=False)

        # remove the temporary file
        os.remove(dummyfile)

    return refdic, adsrefs


# string of logical expressions for use in regex parser
LOGEXPRS = (r'(\bAND\b'        # logical AND
            r'|\band\b'        # logical AND
            r'|\&\&'           # logical AND
            r'|\bOR\b'         # logical OR
            r'|\bor\b'         # logical OR
            r'|\|\|'           # logical OR
            r'|!='             # not equal to
            r'|=='             # equal to
            r'|<='             # less than or equal to
            r'|>='             # greater than or equal to
            r'|<'              # less than
            r'|>'              # greater than
            r'|\('             # left opening bracket
            r'|\)'             # right closing bracket
            r'|\bNOT\b'        # logical NOT
            r'|\bnot\b'        # logical NOT
            r'|!'              # logical NOT
            r'|~'              # logical NOT
            r'|\bASSOC\b'      # pulsar association
            r'|\bassoc\b'      # pulsar association
            r'|\bTYPE\b'       # pulsar type
            r'|\btype\b'       # pulsar type
            r'|\bBINCOMP\b'    # pulsar binary companion type
            r'|\bbincomp\b'    # pulsar binary companion type
            r'|\bEXIST\b'      # pulsar parameter exists in the catalogue
            r'|\bexist\b'      # pulsar parameter exists in the catalogue
            r'|\bERROR\b'      # condition on parameter error
            r'|\berror\b)')    # condition on parameter error


def condition(table, expression, exactMatch=False):
    """
    Apply a logical expression to a table of values.

    Args:
        table (:class:`astropy.table.Table` or :class:`pandas.DataFrame`): a
            table of pulsar data
        expression (str, :class:`~numpy.ndarray`): a string containing a set of
            logical conditions with respect to pulsar parameter names (also
            containing `conditions
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#condition>`_
            allowed when accessing the ATNF Pulsar Catalogue), or a boolean
            array of the same length as the table.
        exactMatch (bool): set to true to exactly match some string comparison
            expressions, e.g., if asking for `'ASSOC(SNR)'` and `exactMatch` is
            True then only pulsar with an association that is just `'SNR'` will
            be returned, whereas if it is False then there could be multiple
            associations including `'SNR'`.

    Returns:
        :class:`astropy.table.Table` or :class:`pandas.DataFrame`: the table of
        values conforming to the input condition. Depending on the type of
        input table the returned table will either be a
        :class:`astropy.table.Table` or :class:`pandas.DataFrame`.

    Example:
        Some examples of this might are:

        1. finding all pulsars with frequencies greater than 100 Hz

        >>> newtable = condition(psrtable, 'F0 > 100')

        2. finding all pulsars with frequencies greater than 50 Hz and
        period derivatives less than 1e-15 s/s.

        >>> newtable = condition(psrtable, '(F0 > 50) & (P1 < 1e-15)')

        3. finding all pulsars in binary systems

        >>> newtable = condition(psrtable, 'TYPE(BINARY)')

        4. parsing a boolean array equivalent to the first example

        >>> newtable = condition(psrtable, psrtable['F0'] > 100)

    """

    from astropy.table import Table
    from pandas import DataFrame

    if expression is None:
        return table

    # check if expression is just a boolean array
    if isinstance(expression, np.ndarray):
        if expression.dtype != np.bool:
            raise TypeError("Numpy array must be a boolean array")
        elif len(expression) != len(table):
            raise Exception("Boolean array and table must be the same length")
        else:
            return table[expression]
    else:
        if not isinstance(expression, string_types):
            raise TypeError("Expression must be a boolean array or a string")
        else:
            if len(expression) == 0:
                return table

    # parse the expression string and split into tokens
    reg = re.compile(LOGEXPRS)
    tokens = reg.split(expression)
    tokens = [t.strip() for t in tokens if t.strip() != '']

    if isinstance(table, Table):
        # convert astropy table to pandas DataFrame
        tab = table.to_pandas()
    elif not isinstance(table, DataFrame):
        raise TypeError("Table must be a pandas DataFrame or astropy Table")
    else:
        tab = table

    matchTypes = ['ASSOC', 'TYPE', 'BINCOMP', 'EXIST', 'ERROR']

    # parse through tokens and replace as required
    ntokens = len(tokens)
    newtokens = []
    i = 0
    while i < ntokens:
        if tokens[i] in [r'&&', r'AND', r'and']:
            # replace synonyms for '&' or 'and'
            newtokens.append(r'&')
        elif tokens[i] in [r'||', r'OR', r'or']:
            # replace synonyms for '|' or 'or'
            newtokens.append(r'|')
        elif tokens[i] in [r'!', r'NOT', r'not']:
            # replace synonyms for '~'
            newtokens.append(r'~')
        elif tokens[i].upper() in matchTypes:
            if ntokens < i+3:
                warnings.warn("A '{}' must be followed by a '(NAME)': "
                              "ignoring in query".format(tokens[i].upper()),
                              UserWarning)
            elif tokens[i+1] != '(' or tokens[i+3] != ')':
                warnings.warn("A '{}' must be followed by a '(NAME)': "
                              "ignoring in query".format(tokens[i].upper()),
                              UserWarning)
            else:
                if tokens[i].upper() == 'ASSOC':
                    if 'ASSOC' not in tab.columns:
                        warnings.warn("'ASSOC' parameter not in table: "
                                      "ignoring in query", UserWarning)
                    elif exactMatch:
                        newtokens.append(r'(ASSOC == "{}")'.format(tokens[i+2].upper()))
                    else:
                        assoc = np.array([tokens[i+2] in str(a) for a in table['ASSOC']])
                        newtokens.append(r'(@assoc)')
                        i += 1
                elif tokens[i].upper() == 'TYPE':
                    if tokens[i+2].upper() == 'BINARY':
                        if 'BINARY' not in tab.keys():
                            warnings.warn("'BINARY' parameter not in table: "
                                          "ignoring in query", UserWarning)
                        else:
                            binary = ~tab['BINARY'].isna()
                            newtokens.append(r'(@binary)')
                            i += 1
                    else:
                        if 'TYPE' not in tab.keys():
                            warnings.warn("'TYPE' parameter not in table: ignoring in query",
                                          UserWarning)
                        elif exactMatch:
                            newtokens.append(r'(TYPE == "{}")'.format(tokens[i+2].upper()))
                        else:
                            ttype = np.array([tokens[i+2] in str(a) for a in table['TYPE']])
                            newtokens.append(r'(@ttype)')
                            i += 1
                elif tokens[i].upper() == 'BINCOMP':
                    if 'BINCOMP' not in tab.columns:
                        warnings.warn("'BINCOMP' parameter not in table: ignoring in query",
                                      UserWarning)
                    elif exactMatch:
                        newtokens.append(r'(BINCOMP == "{}")'.format(tokens[i+2].upper()))
                    else:
                        bincomp = np.array([tokens[i+2] in str(a) for a in table['BINCOMP']])
                        newtokens.append(r'(@bincomp)')
                        i += 1
                elif tokens[i].upper() == 'EXIST':
                    if tokens[i+2].upper() not in tab.columns:
                        warnings.warn("'{}' does not exist for any pulsar".format(tokens[i+2]),
                                      UserWarning)
                        # create an empty DataFrame
                        tab = DataFrame(columns=tab.columns)
                        break
                    else:
                        exists = ~tab[tokens[i+2].upper()].isna()
                        newtokens.append(r'(@exists)')
                        i += 1
                elif tokens[i].upper() == 'ERROR':
                    if tokens[i+2].upper()+'_ERR' not in tab.columns:
                        warnings.warn("Error value for '{}' not present: "
                                      "ignoring in query".format(tokens[i+2]),
                                      UserWarning)
                    else:
                        newtokens.append(r'{}_ERR'.format(tokens[i+2].upper()))
            i += 2
        else:
            newtokens.append(tokens[i].upper())

        i += 1

    # evaluate the expression
    try:
        newtab = tab.query(''.join(newtokens))
    except RuntimeError:
        raise RuntimeError("Could not parse the query")

    if isinstance(table, Table):
        # convert back to an astropy table
        newtab = Table.from_pandas(newtab)

        # re-add any units/types
        for key in table.colnames:
            newtab.columns[key].unit = table.columns[key].unit
            newtab[key] = newtab[key].astype(table[key].dtype)

    return newtab


def characteristic_age(period, pdot, braking_idx=3.):
    """
    Function defining the characteristic age of a pulsar. Returns the
    characteristic age in years using

    .. math::

       \\tau = \\frac{P}{\\dot{P}(n-1)}

    NaNs are returned for any negative period derivates, or NaN imput values.

    Args:
        period (float, array_like): the pulsar period in seconds
        pdot (float, array_like): the pulsar period derivative
        braking_idx (float): the pulsar braking index (defaults to :math:`n=3`)

    Returns:
        float: the characteristic age in years
    """

    # try converting period and pdot to numpy arrays
    try:
        periodarr = np.array(period).flatten()
        pdotarr = np.array(pdot).flatten()
    except Exception as e:
        raise ValueError("Could not convert period/pdot to "
                         "array: {}".format(str(e)))

    assert periodarr.dtype == np.float, "Periods must be floats"
    assert pdotarr.dtype == np.float, "Period derivatives must be floats"
    assert len(periodarr) == len(pdotarr), "Period and derivative arrays must be equal lengths"
    assert braking_idx > 1., "Braking index must be greater than 1"

    # check everything is positive, otherwise return NaN
    age = np.full(len(periodarr), np.nan)
    with np.errstate(invalid='ignore'):
        idx = np.isfinite(pdotarr) & np.isfinite(periodarr) & (pdotarr > 0.)

    age[idx] = (periodarr[idx]/(pdotarr[idx] * (braking_idx - 1.))) / (365.25*86400.)

    # if period and period derivates were just floats return values rather
    # than arrays
    if isinstance(pdot, float) and isinstance(period, float):
        return age[0]

    return age


def age_pdot(period, tau=1e6, braking_idx=3.):
    """
    Function returning the period derivative for a pulsar with a given period
    and characteristic age, using

    .. math::

       \\dot{P} = \\frac{P}{\\tau(n - 1)}

    Args:
        period (list, :class:`numpy.ndarray`): the pulsar period in seconds
        tau (float): the characteristic age in years
        braking_idx (float): the pulsar braking index (defaults to :math:`n=3`)

    Returns:
        :class:`numpy.ndarray`: an array of period derivatives.
    """

    periods = period
    if not isinstance(periods, np.ndarray):
        periods = np.array(periods)

    taus = tau*365.25*86400.  # tau in seconds

    pdots = (periods/(taus * (braking_idx - 1.)))
    pdots[pdots < 0] = np.nan  # set any non zero values to NaN

    return pdots


def B_field(period, pdot):
    """
    Function defining the polar magnetic field strength at the surface of the
    pulsar in gauss (Equation 5.12 of Lyne & Graham-Smith, Pulsar Astronmy, 2nd
    edition) with

    .. math::

       B = 3.2\\!\\times\\!10^{19} \\sqrt{P\\dot{P}}

    NaNs are returned for any negative period derivates, or NaN imput values.

    Args:
        period (float): a pulsar period (s)
        pdot (float): a period derivative

    Returns:
        float: the magnetic field strength in gauss.
    """

    # try converting period and pdot to numpy arrays
    try:
        periodarr = np.array(period).flatten()
        pdotarr = np.array(pdot).flatten()
    except Exception as e:
        raise ValueError("Could not convert period/pdot to "
                         "array: {}".format(str(e)))

    assert periodarr.dtype == np.float, "Periods must be floats"
    assert pdotarr.dtype == np.float, "Period derivatives must be floats"
    assert len(periodarr) == len(pdotarr), "Period and derivative arrays must be equal lengths"

    # check pdot is positive, otherwise return NaN
    bfield = np.full(len(periodarr), np.nan)
    with np.errstate(invalid='ignore'):
        idx = np.isfinite(pdotarr) & np.isfinite(periodarr) & (pdotarr > 0.)
    bfield[idx] = 3.2e19 * np.sqrt(periodarr[idx] * pdotarr[idx])

    # if period and period derivates were just floats return values rather
    # than arrays
    if isinstance(pdot, float) and isinstance(period, float):
        return bfield[0]

    return bfield


def B_field_pdot(period, Bfield=1e10):
    """
    Function to get the period derivative from a given pulsar period and
    magnetic field strength using

    .. math::

       \\dot{P} = \\frac{1}{P}\\left( \\frac{B}{3.2\\!\\times\\!10^{19}} \\right)^2

    Args:
        period (list, :class:`~numpy.ndarray`): a list of period values
        Bfield (float): the polar magnetic field strength (Defaults to
            :math:`10^{10}` G)

    Returns:
        :class:`numpy.ndarray`: an array of period derivatives
    """

    periods = period
    if not isinstance(periods, np.ndarray):
        periods = np.array(periods)

    pdots = (Bfield/3.2e19)**2/periods
    pdots[pdots < 0] = np.nan  # set any non zero values to NaN

    return pdots


def death_line(logP, linemodel='Ip', rho6=1.):
    """
    The pulsar death line. Returns the base-10 logarithm of the period
    derivative for the given values of the period.

    Args:
        logP (list, :class:`~numpy.ndarray`): the base-10 log values of period.
        linemodel (str): a string with one of the above model names. Defaults
            to ``'Ip'``.
        rho6 (float): the value of the :math:`\\rho_6` parameter from [ZHM]_ .
            Defaults to 1 is, which is equivalent to :math:`10^6` cm.

    Returns:
        :class:`numpy.ndarray`: a vector of period derivative values

    .. note::

        The death line models can be:

        * 'I' - Equation 3 of [ZHM]
        * 'Ip' - Equation 4 of [ZHM]
        * 'II' - Equation 5 of [ZHM]
        * 'IIp' - Equation 6 of [ZHM]
        * 'III' - Equation 8 of [ZHM]
        * 'IIIp' - Equation 9 of [ZHM]
        * 'IV' - Equation 10 of [ZHM]
        * 'IVp' - Equation 11 of [ZHM]

    .. [ZHM] Zhang, Harding & Muslimov, *ApJ*, **531**, L135-L138 (2000),
        `arXiv:astro-ph/0001341 <https://arxiv.org/abs/astro-ph/0001341>`_

    """

    gradvals = {'I': (11./4),
                'Ip': (9./4.),
                'II': (2./11.),
                'IIp': -(2./11.),
                'III': (5./2.),
                'IIIp': 2.,
                'IV': -(3./11.),
                'IVp': -(7./11.)}
    intercept = {'I': 14.62,
                 'Ip': 16.58,
                 'II': 13.07,
                 'IIp': 14.50,
                 'III': 14.56,
                 'IIIp': 16.52,
                 'IV': 15.36,
                 'IVp': 16.79}
    rho = {'I': 0.,
           'Ip': 1.,
           'II': 0.,
           'IIp': (8./11.),
           'III': 0.,
           'IIIp': 1.,
           'IV': 0.,
           'IVp': (8./11.)}

    lp = logP
    if not isinstance(lp, np.ndarray):
        lp = np.array(lp)

    return lp*gradvals[linemodel] - intercept[linemodel] + rho[linemodel]*np.log10(rho6)


def label_line(ax, line, label, color='k', fs=14, frachoffset=0.1):
    """
    Add an annotation to the given line with appropriate placement and
    rotation.

    Based on code from `"How to rotate matplotlib annotation to match a line?"
    <http://stackoverflow.com/a/18800233/230468>`_ and `this
    <https://stackoverflow.com/a/38414616/1862861>`_ answer.

    Args:
        ax (:class:`matplotlib.axes.Axes`): Axes on which the label should be
            added.
        line (:class:`matplotlib.lines.Line2D`): Line which is being labeled.
        label (str): Text which should be drawn as the label.
        color (str): a color string for the label text. Defaults to ``'k'``
        fs (int): the font size for the label text. Defaults to 14.
        frachoffset (float): a number between 0 and 1 giving the fractional
            offset of the label text along the x-axis. Defaults to 0.1, i.e.,
            10%.

    Returns:
        :class:`matplotlib.text.Text`: an object containing the label
        information

    """
    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    # use fractional horizontal offset frachoffset to set the x position of the label by default
    # otherwise use the halign value
    if frachoffset >= 0 and frachoffset <= 1:
        if ax.get_xscale() == 'log':
            xx = np.log10(x1) + frachoffset*(np.log10(x2) - np.log10(x1))
        else:
            xx = x1 + frachoffset*(x2 - x1)
    else:
        raise ValueError("frachoffset must be between 0 and 1")

    if ax.get_xscale() == 'log' and ax.get_yscale() == 'log':
        yy = np.interp(xx, np.log10(xdata), np.log10(ydata))
        xx = 10**xx
        yy = 10**yy
    elif ax.get_xscale() == 'log' and ax.get_yscale() != 'log':
        yy = np.interp(xx, np.log10(xdata), ydata)
        xx = 10**xx
    else:
        yy = np.interp(xx, xdata, ydata)

    ylim = ax.get_ylim()
    xytext = (0, 5)
    text = ax.annotate(label, xy=(xx, yy), xytext=xytext, textcoords='offset points',
                       size=fs, color=color, zorder=1,
                       horizontalalignment='left', verticalalignment='center_baseline')

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation_mode('anchor')
    text.set_rotation(slope_degrees)
    ax.set_ylim(ylim)
    return text
