"""
A selection of useful functions used by the module.
"""

import warnings
import re
import os
import functools
import numpy as np
import requests
import tarfile
from bs4 import BeautifulSoup

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as aunits
from astropy.utils.data import download_file, clear_download_cache
from pandas import DataFrame

from .config import (
    ATNF_BASE_URL,
    ADS_URL,
    ATNF_TARBALL,
    ATNF_VERSION_TARBALL,
    MSP_URL,
    PSR_ALL,
    PSR_ALL_PARS,
    GLITCH_URL,
    GC_URL,
)


def get_catalogue(path_to_db=None, cache=True, update=False, pandas=False, version="latest"):
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
        version (str): the version string (without the leading "v") of the ATNF
            catalogue version to download. This defaults to "latest" to get the
            most up-to-date version.

    Returns:
        :class:`~astropy.table.Table` or :class:`~pandas.DataFrame`: a table
        containing the entire catalogue.

    """

    if path_to_db is None:
        # remove any cached file if requested
        if update:
            if check_update():
                clear_download_cache(ATNF_TARBALL)

        if version == "latest":
            atnftarball = ATNF_TARBALL
        else:
            atnftarball = ATNF_VERSION_TARBALL.format(version)

        # get the tarball
        try:
            dbtarfile = download_file(atnftarball, cache=cache)
        except IOError:
            raise IOError(
                "Problem accessing ATNF catalogue tarball: {}".format(atnftarball)
            )

        try:
            # open tarball
            pulsargz = tarfile.open(dbtarfile, mode="r:gz")

            # extract the database file
            dbfile = pulsargz.extractfile("psrcat_tar/psrcat.db")
        except IOError:
            raise IOError("Problem extracting the database file")
    else:
        try:
            dbfile = open(path_to_db, "r")
        except IOError:
            raise IOError("Error loading given database file")

    breakstring = "@"  # break between each pulsar
    commentstring = "#"  # specifies line is a comment

    # create list of dictionaries - one for each pulsar
    psrlist = [{}]

    version = None  # catalogue version

    # loop through lines in dbfile
    for line in dbfile.readlines():
        if isinstance(line, str):
            dataline = line.split()
        else:
            dataline = line.decode().split()  # Splits on whitespace

        if dataline[0][0] == commentstring:
            # get catalogue version (should be in first comment string)
            if dataline[0] == "#CATALOGUE" and len(dataline) == 2:
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
                val = dataline[1].split(":")[-1]  # account for RA and DEC strings

                try:
                    float(val)
                except ValueError:
                    raise ValueError("Value with error is not convertable to a float")

                if dataline[2][0] == "-" or "." in dataline[2]:
                    # negative errors or those with decimal points are absolute values
                    scalefac = 1.0
                else:
                    # split on exponent
                    valsplit = re.split("e|E|d|D", val)
                    scalefac = 1.0
                    if len(valsplit) == 2:
                        scalefac = 10 ** (-int(valsplit[1]))

                    dpidx = valsplit[0].find(".")  # find position of decimal point
                    if dpidx != -1:  # a point is found
                        scalefac *= 10 ** (len(valsplit[0]) - dpidx - 1)

                # add error column if required
                psrlist[-1][dataline[0] + "_ERR"] = (
                    float(dataline[2]) / scalefac
                )  # error entry
            else:
                # add reference column if required
                psrlist[-1][dataline[0] + "_REF"] = dataline[2]  # reference entry

            if len(dataline) > 3:
                # last entry must(!) be a reference
                psrlist[-1][dataline[0] + "_REF"] = dataline[3]  # reference entry

    dbfile.close()  # close tar file
    if not path_to_db:
        pulsargz.close()

    del psrlist[-1]  # Final breakstring comes at the end of the file

    # add RA and DEC in degs and JNAME/BNAME
    for i, psr in enumerate(list(psrlist)):
        # add 'JNAME', 'BNAME' and 'NAME'
        if "PSRJ" in psr.keys():
            psrlist[i]["JNAME"] = psr["PSRJ"]
            psrlist[i]["NAME"] = psr["PSRJ"]
            if "PSRJ_REF" in psr.keys():
                psrlist[i]["JNAME_REF"] = psr["PSRJ_REF"]
                psrlist[i]["NAME_REF"] = psr["PSRJ_REF"]

        if "PSRB" in psr.keys():
            psrlist[i]["BNAME"] = psr["PSRB"]
            if "PSRB_REF" in psr.keys():
                psrlist[i]["BNAME_REF"] = psr["PSRB_REF"]

            if "NAME" not in psrlist[i].keys():
                psrlist[i]["NAME"] = psr["PSRB"]
                if "PSRB_REF" in psr.keys():
                    psrlist[i]["NAME_REF"] = psr["PSRB_REF"]

        if "RAJ" in psr.keys() and "DECJ" in psr.keys():
            # check if the string can be converted to a float (there are a few
            # cases where the position is just a decimal value)
            try:
                rad = float(psr["RAJ"])
                ras = Angle(rad * aunits.hourangle)
                psr["RAJ"] = ras.to_string(sep=":", pad=True)
            except ValueError:
                pass

            try:
                decd = float(psr["DECJ"])
                decs = Angle(decd * aunits.deg)
                psr["DECJ"] = decs.to_string(sep=":", pad=True, alwayssign=True)
            except ValueError:
                pass

            try:
                coord = SkyCoord(
                    psr["RAJ"], psr["DECJ"], unit=(aunits.hourangle, aunits.deg)
                )
                psrlist[i]["RAJD"] = coord.ra.deg  # right ascension in degrees
                psrlist[i]["DECJD"] = coord.dec.deg  # declination in degrees
            except Exception as e:
                warnings.warn(
                    "Error converting RAJ/DECJ strings to degrees for {}: {}".format(
                        psrlist[i]["NAME"], e
                    )
                )
                psrlist[i]["RAJD"] = np.nan
                psrlist[i]["DECJD"] = np.nan

            # set errors on positions in degrees
            if "RAJ_ERR" in psr.keys():
                # get the units for the error
                nvals = len(psr["RAJ"].split(":"))
                rajscale = 1.0 / (60.0 ** (len(psr["RAJ"].split(":")) - 1))

                psrlist[i]["RAJD_ERR"] = (
                    (psr["RAJ_ERR"] * rajscale * aunits.hourangle).to("deg").value
                )

            if "DECJ_ERR" in psr.keys():
                # get the units for the error
                nvals = len(psr["DECJ"].split(":"))
                if nvals == 1:
                    decunit = aunits.deg
                elif nvals == 2:
                    decunit = aunits.arcmin
                else:
                    decunit = aunits.arcsec

                psrlist[i]["DECJD_ERR"] = (psr["DECJ_ERR"] * decunit).to("deg").value

        # always include a value for POSEPOCH if PEPOCH is present
        # this matches the behaviour of the psrcat v1.66 CLI (definePosEpoch.c:202-213)
        if "POSEPOCH" not in psr.keys() and "PEPOCH" in psr.keys():
            psrlist[i]["POSEPOCH"] = psrlist[i]["PEPOCH"]
            if "PEPOCH_REF" in psr.keys():
                psrlist[i]["POSEPOCH_REF"] = psrlist[i]["PEPOCH_REF"]

        # derive the 'DATE' parameter from the 'POSEPOCH' reference
        # implementation based upon the psrcat v1.66 CLI (definePosEpoch.c:214-242)
        if "POSEPOCH_REF" in psr.keys():
            try:
                refyear = int(re.search(r"(\d{2})(?!.*\d)", psrlist[i]["POSEPOCH_REF"]).group(0))
                if refyear > 65:
                    refyear += 1900
                else:
                    refyear += 2000
                psrlist[i]['DATE'] = refyear
            except (AttributeError, TypeError):
                pass

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
            if PSR_ALL[key]["units"]:
                psrtable.columns[key].unit = PSR_ALL[key]["units"]

                if PSR_ALL[key]["err"] and key + "_ERR" in psrtable.colnames:
                    psrtable.columns[key + "_ERR"].unit = PSR_ALL[key]["units"]

    # add metadata
    if not path_to_db:
        if version is not None:
            psrtable.meta["version"] = version
        else:
            psrtable.meta["version"] = None
            warnings.warn("No version number found in the database file", UserWarning)
        psrtable.meta["ATNF Pulsar Catalogue"] = ATNF_BASE_URL

    if path_to_db:
        psrtable.meta["Database file"] = path_to_db

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
        soup = BeautifulSoup(gt.content, "html.parser")
    except Exception as e:
        warnings.warn(
            "Count not parse the glitch catalogue: {}".format(str(e)), UserWarning
        )
        return None

    # get table rows
    rows = soup.table.find_all("tr")

    # set the table headings
    tabledict = dict()
    tabledict["NAME"] = []
    tabledict["JNAME"] = []
    tabledict["Glitch number"] = []
    tabledict["MJD"] = []
    tabledict["MJD_ERR"] = []
    tabledict["DeltaF/F"] = []
    tabledict["DeltaF/F_ERR"] = []
    tabledict["DeltaF1/F1"] = []
    tabledict["DeltaF1/F1_ERR"] = []
    tabledict["Reference"] = []

    # loop through rows: rows with glitches have their first column as an index
    for row in rows:
        tds = row.find_all("td")

        if tds[0].contents[0].string is None:
            continue

        try:
            tabledict["NAME"].append(tds[1].contents[0].string)
            jname = (
                "J" + tds[2].contents[0].string
                if "J" != tds[2].contents[0].string[0]
                else tds[2].contents[0].string
            )
            tabledict["JNAME"].append(jname)
            tabledict["Glitch number"].append(int(tds[3].contents[0].string))

            for j, pname in enumerate(
                [
                    "MJD",
                    "MJD_ERR",
                    "DeltaF/F",
                    "DeltaF/F_ERR",
                    "DeltaF1/F1",
                    "DeltaF1/F1_ERR",
                ]
            ):
                try:
                    val = float(tds[4 + j].contents[0].string)
                except ValueError:
                    val = np.nan

                tabledict[pname].append(val)

            # get reference link if present
            try:
                ref = tds[10].contents[0].a.attrs["href"]
            except AttributeError:
                ref = tds[10].contents[0].string
            tabledict["Reference"].append(ref)
        except RuntimeError:
            warnings.warn("Problem parsing glitch table", UserWarning)
            return None

    # convert to an astropy table
    table = Table(tabledict)
    table.columns["MJD"].unit = aunits.Unit("d")  # add units of days to glitch time
    table.columns["MJD_ERR"].unit = aunits.Unit("d")

    # correct scaling of parameters
    table["DeltaF/F"] *= 1e-9
    table["DeltaF/F_ERR"] *= 1e-9
    table["DeltaF1/F1"] *= 1e-3
    table["DeltaF1/F1_ERR"] *= 1e-3

    if psr is None:
        return table
    else:
        if psr not in table["NAME"] and psr not in table["JNAME"]:
            warnings.warn(
                "Pulsar '{}' not found in glitch catalogue".format(psr), UserWarning
            )
            return None
        else:
            if psr in table["NAME"]:
                return table[table["NAME"] == psr]
            else:
                return table[table["JNAME"] == psr]


def get_gc_catalogue():
    """
    Download and parse Paolo Freire's table of pulsars in globular clusters
    from http://www.naic.edu/~pfreire/GCpsr.txt. This will be returned as a
    :class:`astropy.table.Table`, but with some additional methods to extract
    information for individual clusters. If the webpage returned an error
    ``None`` is returned.

    The additional methods of the returned :class:`~astropy.table.Table` are:

    * ``clusters``: by default this returns a list if the common names of the
      globular clusters in the table, or using the ``ngc=True`` argument will
      return the NGC names of all clusters.
    * ``cluster_pulsars``: returns a sub-table only containing pulsars within
      a supplied globular cluster name.
    """

    # get the webpage
    try:
        gct = requests.get(GC_URL)
    except Exception as e:
        raise RuntimeError(
            "Error downloading the globular cluster pulsar table: {}".format(str(e))
        )

    if gct.status_code != 200:
        warnings.warn("Count not query the globular cluster pulsar table.", UserWarning)
        return None

    # convert from btyes to utf-8
    tabledata = gct.content.decode("utf-8")

    names = [
        "Cluster",
        "NGC",
        "Pulsar",
        "Offset",
        "Period",
        "dP/dt",
        "dP/dt error",
        "DM",
        "DM error",
        "Pb",
        "x",
        "x error",
        "e",
        "e error",
        "m2",
    ]

    # create astropy table
    gctable = Table(
        names=names,
        dtype=[
            str,
            str,
            str,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
        ],
        units=[
            None,
            None,
            None,
            aunits.arcmin,
            1e-3 * aunits.s,
            1e-20 * aunits.s / aunits.s,
            1e-20 * aunits.s / aunits.s,
            aunits.cm ** -3 * aunits.pc,
            aunits.cm ** -3 * aunits.pc,
            aunits.d,
            aunits.s,
            aunits.s,
            None,
            None,
            aunits.Msun,
        ],
        masked=True,
    )

    # loop over data and add in rows
    for line in tabledata.split("\n"):
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == "#":  # comment line
            continue

        # make sure any en-dashes "−" are replaced with -
        linevalues = re.sub("−", "-", line).split()

        # check if line is a pulsar
        ispulsar = False
        if linevalues[0][0] in ["J", "B"]:
            try:
                if linevalues[0][5] in ["+", "-"]:
                    ispulsar = True
            except IndexError:
                pass

        # get globular cluster name
        if not ispulsar:
            # check if there's an NGC name
            ngc = re.search(r"NGC\s\d+", line)

            if ngc is not None:
                ngc = ngc.group()

            # common name
            gcname = line.split("(")[0].strip()
            continue

        psrentry = {key: None for key in names}
        psrentry["Cluster"] = gcname
        psrentry["NGC"] = str(ngc)
        psrentry["Pulsar"] = linevalues[0]

        psrentry["Offset"] = float(linevalues[1]) if linevalues[1] != "*" else np.nan
        psrentry["Period"] = float(linevalues[2]) if linevalues[2] != "*" else np.nan

        def parse_value_error(value):
            """
            Parse a string value and extract the error if given.

            Args:
                value (str): A string containing a value with associated last
                    digit error in parentheses, e.g., "-1.2324(5)"

            Returns:
                tuple: A tuple containing the value and the error as floats, or
                    NaN for the error is not present
            """
            if "(" in value:
                error = float(value[value.find("(") + 1 : value.find(")")])
                exponent = value.find("(") - value.find(".") - 1

                if "*10-15" in value:
                    # scale to units of 10-20
                    scale = 1e-15 / 1e-20
                else:
                    scale = 1.0

                val = float(value[: value.find("(")]) * scale
                error *= 10 ** -exponent * scale
            else:
                if "*10-15" in value:
                    # scale to units of 10-20
                    scale = 1e-15 / 1e-20
                    value = value.rstrip("*10-15")
                else:
                    scale = 1.0

                # remove any inequalities
                val = (
                    float(value.strip(">").strip("<"))
                    if value != "*" and value != "i"
                    else np.nan
                )
                error = np.nan

            return val, error

        # dP/dt
        val, error = parse_value_error(linevalues[3])
        psrentry["dP/dt"] = val
        psrentry["dP/dt error"] = error

        # DM
        val, error = parse_value_error(linevalues[4])
        psrentry["DM"] = val
        psrentry["DM error"] = error

        # binary parameters
        val, _ = parse_value_error(linevalues[5])
        psrentry["Pb"] = val

        val, error = parse_value_error(linevalues[6])
        psrentry["x"] = val
        psrentry["x error"] = error

        val, error = parse_value_error(linevalues[7])
        psrentry["e"] = val
        psrentry["e error"] = error

        # due to an entry without proper spacing (J1953+1846A) one line gives an index error for m2
        try:
            val, _ = parse_value_error(linevalues[8])
            psrentry["m2"] = val
        except IndexError:
            psrentry["m2"] = np.nan

        mask = {}
        for key, value in psrentry.items():
            mask[key] = (
                True if value is None or value == "None" or value is np.nan else False
            )

        gctable.add_row(psrentry, mask=mask)

    class GCTable(Table):
        def __init__(self, data=None, **kwargs):
            super().__init__(data=data, **kwargs)

        def clusters(self, ngc=False):
            """
            Return a list of the globular cluster names.

            Args:
                ngc (bool): Set to True to return the NGC names rather than
                    common names (some of which may be NGC names anyway!)
            """

            if ngc:
                return [val[0] for val in self.group_by("NGC").groups.keys.as_array()]
            else:
                return [
                    val[0] for val in self.group_by("Cluster").groups.keys.as_array()
                ]

        def cluster_pulsars(self, cluster):
            """
            Return a table just containing the pulsars in a requested globular
            cluster. If the cluster is not found then return None.

            Args:
                cluster (str): The name of the globular cluster for which to
                    return the pulsars.

            Returns:
                :class:`~astropy.table.Table`: a table the pulsars in the given
                    cluster.
            """

            clustercolumns = ["Cluster", "NGC"]

            # check if name is just an integer value, which could be an NGC id
            try:
                gcnum = int(cluster)
            except ValueError:
                gcnum = None

            if gcnum is not None:
                cluster = "NGC {}".format(gcnum)

            # check for cluster name in the two columns containing name info
            gctable = None
            for cc in clustercolumns:
                # remove spaces from cluster names, so, e.g., 47Tuc will work
                # as well as 47 Tuc
                cluster = cluster.replace(" ", "")
                clusters = [
                    gc.replace(" ", "") for gc in self[cc] if isinstance(gc, str)
                ]

                if cluster in clusters:
                    mask = clusters == cluster
                    gctable = self[mask]
                    break

            return gctable

    # return the table
    return GCTable(gctable)


def get_msp_catalogue():
    """
    Download and parse Dunc Lorimer's `catalogue
    <http://astro.phys.wvu.edu/GalacticMSPs/GalacticMSPs.txt>`_ of galactic
    millisecond pulsars. This makes use of code from the now deleted `galmsps
    <https://github.com/astrogewgaw/galmsps>`_ respository that scraped the
    database into JSON every month, contributed by `@astrogewgaw
    <https://github.com/astrogewgaw>`_ before deletion via `PR #88
    <https://github.com/mattpitkin/psrqpy/pull/88>`_.

    Returns:
        :class:`~astropy.table.Table`: a table the MSPs.
    """

    import re

    # get the webpage
    try:
        mspt = requests.get(MSP_URL)
    except Exception as e:
        raise RuntimeError("Error downloading the MSP table: {}".format(str(e)))

    if mspt.status_code != 200:
        warnings.warn("Count not query the MSP table.", UserWarning)
        return None

    numeric = (
        lambda _: float(_)
        if re.search(re.compile(r"^[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$"), _)
        else None
    )

    # parse the text file into a dict
    tabledata = {
        str(i + 1): {
            key: conv(value)  # type: ignore
            for (key, conv), value in zip(
                [
                    ("NAME", str),
                    ("P0", numeric),
                    ("DM", numeric),
                    ("GL", numeric),
                    ("GB", numeric),
                    ("PB", numeric),
                    ("A1", numeric),
                    ("DISCOVERY YEAR", int),
                    ("NOTES", str),
                ],
                values,
            )
        }
        for i, values in enumerate(
            [
                [_ for _ in re.split(r"\s+", _) if _]
                for _ in re.split(
                    r"\n+",
                    re.sub(
                        re.compile(r"^([^#]*)[#](.*)$", re.MULTILINE),
                        "",
                        mspt.content.decode(),
                    ).strip(),
                )
            ]
        )
    }

    # column names
    colnames = [
        "NAME",
        "P0",  # period
        "DM",  # dispersion measure
        "GL",  # galactic longitude
        "GB",  # galactic latitude,
        "PB",  # binary period
        "A1",  # projected semi-major axis
        "DISCOVERY YEAR",  # discovery year
        "NOTES",
    ]

    # units for the columns
    units = {
        "NAME": None,
        "P0": aunits.s * 1e-3,
        "DM": aunits.cm ** -3 * aunits.pc,
        "GL": aunits.deg,
        "GB": aunits.deg,
        "PB": aunits.d,
        "A1": aunits.s,
        "DISCOVERY YEAR": aunits.year,
        "NOTES": None,
    }

    # create dictionary in form usable by astropy Table
    tabledict = {
        colnames[i]: [item[colnames[i]] for item in tabledata.values()]
        for i in range(len(colnames))
    }

    # convert into astropy Table
    msptable = Table(
        data=tabledict,
        units=units,
        dtype=[str, float, float, float, float, float, float, int, str],
        masked=True,
    )

    return msptable


def check_old_references(func):
    @functools.wraps(func)
    def wrapper_check_old_references(*args, **kwargs):
        # check version
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if kwargs.get("useads", False):
                from ads.exceptions import APIResponseError

                if type(e) is APIResponseError:
                    raise (e)

            warnings.warn(
                "The way references are stored in the ATNF catalogue has "
                "changed and 'get_references' will no longer parse the old "
                "style references. Please update your cached version of the "
                "ATNF pulsar catalogue with:\n\n"
                ">>> from psrqpy import QueryATNF\n"
                ">>> QueryATNF(checkupdate=True)\n\n"
                "and update any cached references with:\n\n"
                ">>> from psrqpy import get_references\n"
                ">>> refs = get_references(updaterefcache=True)"
            )
            # get the correct number of expected outputs
            output = [None]
            if kwargs.get("useads", False):
                output.append(None)
            elif len(args) > 0:
                if args[0]:
                    output.append(None)

            if kwargs.get("bibtex", False):
                output.append(None)
            elif len(args) == 4:
                if args[3]:
                    output.append(None)

            if kwargs.get("showfails", False):
                output.append(None)
            elif len(args) == 5:
                if args[4]:
                    output.append(None)

            return tuple(output)

    return wrapper_check_old_references


@check_old_references
def get_references(
    useads=False, cache=True, updaterefcache=False, bibtex=False, showfails=False
):
    """
    Return a dictionary of paper
    `reference <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_ref.html>`_
    in the ATNF catalogue. The keys are the ref strings given in the ATNF
    catalogue.

    Note: The way that the ATNF references are stored has changed, so if you
    downloaded the catalogue with a version of psrqpy before v1.0.8 you may
    need to run this function with ``updaterefcache=True`` to allow references
    to work. You may also want to update the ATNF catalogue tarball with:

    >>> import psrqpy
    >>> psrqpy.QueryATNF(checkupdate=True)

    Args:
        useads (bool): boolean to set whether to use the python mod:`ads`
            module to get the NASA ADS URL for the references.
        cache (bool): use cached, or cache, the reference bundled with the
            catalogue tarball.
        updaterefcache (bool): update the cached references.
        bibtex (bool): if using ADS return the bibtex for the reference along
            with the ADS URL.
        showfails (bool): if outputting NASA ADS references set this flag to
            True to output the reference tags of references that fail to be
            found (mainly for debugging purposes).

    Returns:
        dict: a dictionary of references.
    """

    import tempfile
    import json

    # get the tarball
    try:
        dbtarfile = download_file(ATNF_TARBALL, cache=not updaterefcache)
    except IOError:
        raise IOError("Problem accessing ATNF catalogue tarball")

    try:
        # open tarball
        pulsargz = tarfile.open(dbtarfile, mode="r:gz")

        # extract the references
        reffile = pulsargz.extractfile("psrcat_tar/psrcat_ref")
    except IOError:
        raise IOError("Problem extracting the database file")

    refdic = {
        line.split()[0]: " ".join(line.split()[2:])
        for line in reffile.read().decode("utf-8").strip().split("***")
        if len(line) > 0
    }

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
            warnings.warn(
                "Could not import ADS module, so no ADS information "
                "will be included",
                UserWarning,
            )
            return refdic, None

    # try getting cached references
    if not cache:
        adsrefs = {}
    else:
        from astropy.utils.data import is_url_in_cache

        tmpdir = tempfile.gettempdir()  # get system "temporary" directory
        dummyurl = "file://{}/ads_cache".format(tmpdir)
        dummyfile = os.path.join("{}".format(tmpdir), "ads_cache")

        # check if cached ADS refs list exists (using dummy URL)
        if is_url_in_cache(dummyurl) and not updaterefcache:
            adsfile = download_file(dummyurl, cache=True, show_progress=False)

            try:
                fp = open(adsfile, "r")
            except IOError:
                warnings.warn(
                    "Could not load ADS URL cache for references", UserWarning
                )
                return refdic, None

            cachedrefs = json.load(fp)
            fp.close()

            adsrefs = None
            adsbibtex = None
            failures = None
            if "urls" in cachedrefs:
                adsrefs = cachedrefs["urls"]
            if bibtex and "bibtex" in cachedrefs:
                adsbibtex = cachedrefs["bibtex"]
            if showfails and "failures" in cacherefs:
                failures = cachedrefs["failures"]

            if bibtex:
                if failures is None:
                    return refdic, adsrefs, adsbibtex
                else:
                    return refdic, adsrefs, adsbibtex, failures
            else:
                if failures is None:
                    return refdic, adsrefs
                else:
                    return refdic, adsrefs, failures
        else:
            adsrefs = {}

    # loop over references
    j = 0
    bibcodes = {}
    failures = []
    for reftag in refdic:
        j = j + 1

        refstring = refdic[reftag]

        # check if IAU Circular or PhD thesis
        iaucirc = True if "IAU Circ" in refstring else False
        thesis = True if "PhD thesis" in refstring else False

        sepauthors = ""

        # check for arXiv identifier
        arxivid = None
        if "arXiv:" in refstring or "ArXiv:" in refstring:
            for searchterm in [
                r"[Aa]rXiv:[0-9]{4}.[0-9]*",
                r"[Aa]rXiv:astro-ph/[0-9]{7}",
            ]:
                match = re.search(searchterm, refstring)

                if match is not None:
                    arxivid = match.group().lower()
                    break
        else:
            if iaucirc:
                # get circular number (value after IAU Circ. No.)
                spl = re.split(r"([0-9]{4})", refstring)
                noidx = 1
                for val in spl:
                    if "IAU Circ" in val:
                        break
                    noidx += 1
                volume = spl[noidx]
            else:
                # do splitting on the year (allows between 1000-2999)
                spl = re.split(r"([1-2][0-9]{3})", refstring)

                if len(spl) < 2:
                    # no authors + year, so ignore!
                    failures.append(reftag)
                    continue

                year = spl[1] if len(spl[1]) == 4 else None

                try:
                    int(year)
                except (ValueError, TypeError):
                    # "year" is not an integer
                    failures.append(reftag)
                    continue

                # get the authors (remove line breaks/extra spaces and final full-stop)
                authors = spl[0].strip().strip(".")

                # remove " Jr." from any author names (as it causes issues!)
                authors = authors.replace(" Jr.", "")

                # replace ampersands/and with ".," for separation
                authors = authors.replace(" &", ".,").replace(" and", ".,")

                # separate out authors
                sepauthors = [
                    auth.lstrip()
                    for auth in authors.split(".,")
                    if len(auth.strip()) > 0 and "et al" not in auth
                ]

                # remove any "'s for umlauts in author names
                sepauthors = [a.replace(r'"', "") for a in sepauthors]

                if len(sepauthors) == 0:
                    # no authors were parsed
                    failures.append(reftag)
                    continue

            if not thesis and not iaucirc:
                volume = None
                page = None
                if len(spl) > 2:
                    # join the remaining values and split on ","
                    extrainfo = [
                        info
                        for info in ("".join(spl[2:])).lstrip(".").split(",")
                        if len(info.strip()) > 0
                    ]

                    # get the journal volume (assumed to be second from last)
                    try:
                        # in case volume contains issue number in brackets perform split
                        volume = int(extrainfo[-2].strip().split("(")[0])
                    except (IndexError, TypeError, ValueError):
                        # could not get the volume
                        pass

                    # get the page if given (assumed to be th last value)
                    try:
                        testpage = re.sub(
                            r"[\+\-\.]", "", extrainfo[-1].strip().split("-")[0]
                        )
                        if not testpage.startswith(
                            "eaao"
                        ):  # Science Advances page string
                            if (
                                testpage[0].upper() in ["L", "A", "E"]
                                or testpage[0:4] == ""
                            ):  # e.g. for ApJL, A&A, PASA
                                _ = int(testpage[1:])
                            elif testpage[-1].upper() == "P":  # e.g., for early MNRAS
                                _ = int(testpage[:-1])
                            else:
                                _ = int(testpage)
                        page = testpage
                    except (IndexError, TypeError, ValueError):
                        # could not get the page
                        pass

                if volume is None or page is None:
                    failures.append(reftag)
                    continue

        # generate the query string
        if arxivid is None:
            if not thesis:
                if iaucirc:
                    myquery = 'bibstem:"IAUC" volume:"{}"'.format(volume)
                else:
                    # default query without authors
                    myquery = "year:{} AND volume:{} AND page:{}".format(
                        year, volume, page
                    )

                    # add author if given
                    if len(sepauthors) > 0:
                        # check if authors have spaces in last names (a few cases due to formating of some accented names),
                        # if so try next author...
                        for k, thisauthor in enumerate(sepauthors):
                            if len(thisauthor.split(",")[0].split()) == 1:
                                myquery += ' AND author:"{}{}"'.format(
                                    "^" if k == 0 else "", thisauthor
                                )
                                break
            else:
                myquery = 'year: {} AND author:"^{}" AND bibstem:"PhDT"'.format(
                    year, sepauthors[0]
                )
        else:
            myquery = arxivid

        try:
            article = ads.SearchQuery(q=myquery)
        except APIResponseError:
            failures.append(reftag)
            warnings.warn(
                "Could not get reference information, so no ADS "
                "information for {} will be included".format(reftag),
                UserWarning,
            )
            continue

        for paper in article:
            bibcodes[reftag] = paper.bibcode
            adsrefs[reftag] = ADS_URL.format(bibcodes[reftag])

        # check if paper bibcode was found
        if reftag not in bibcodes:
            failures.append(reftag)

    if bibtex:
        # use ExportQuery to get bibtex
        expquery = ads.ExportQuery(list(bibcodes.values())).execute().split("\n\n")

        adsbibtex = {}
        for reftag in bibcodes:
            for equery in expquery:
                if bibcodes[reftag] in equery:
                    adsbibtex[reftag] = equery
                    break

    if cache:
        # output adsrefs to cache file
        try:
            # output to dummy temporary file and then "download" to cache
            fp = open(dummyfile, "w")

            cachedic = {}
            cachedic["urls"] = adsrefs

            if bibtex:
                cachedic["bibtex"] = adsbibtex

            if showfails:
                cachedic["failures"] = failures

            json.dump(cachedic, fp, indent=2)
            fp.close()
        except IOError:
            raise IOError("Could not output the ADS references to a file")

        # cache the file
        _ = download_file(dummyurl, cache=True, show_progress=False)

        # remove the temporary file
        os.remove(dummyfile)

    if bibtex:
        if showfails:
            return refdic, adsrefs, adsbibtex, failures
        else:
            return refdic, adsrefs, adsbibtex
    else:
        if showfails:
            return refdic, adsrefs, failures
        else:
            return refdic, adsrefs


# string of logical expressions for use in regex parser
LOGEXPRS = (
    r"(\bAND\b"  # logical AND
    r"|\band\b"  # logical AND
    r"|\&\&"  # logical AND
    r"|\bOR\b"  # logical OR
    r"|\bor\b"  # logical OR
    r"|\|\|"  # logical OR
    r"|!="  # not equal to
    r"|=="  # equal to
    r"|<="  # less than or equal to
    r"|>="  # greater than or equal to
    r"|<"  # less than
    r"|>"  # greater than
    r"|\("  # left opening bracket
    r"|\)"  # right closing bracket
    r"|\bNOT\b"  # logical NOT
    r"|\bnot\b"  # logical NOT
    r"|!"  # logical NOT
    r"|~"  # logical NOT
    r"|\bASSOC\b"  # pulsar association
    r"|\bassoc\b"  # pulsar association
    r"|\bTYPE\b"  # pulsar type
    r"|\btype\b"  # pulsar type
    r"|\bBINCOMP\b"  # pulsar binary companion type
    r"|\bbincomp\b"  # pulsar binary companion type
    r"|\bSURVEY\b"  # pulsar observation survey
    r"|\bsurvey\b"  # pulsar observation survey
    r"|\bDISCOVERY\b"  # pulsar discovery survey
    r"|\bdiscovery\b"  # pulsar discovery survey
    r"|\bEXIST\b"  # pulsar parameter exists in the catalogue
    r"|\bexist\b"  # pulsar parameter exists in the catalogue
    r"|\bERROR\b"  # condition on parameter error
    r"|\berror\b)"
)  # condition on parameter error


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
        if expression.dtype != bool:
            raise TypeError("Numpy array must be a boolean array")
        elif len(expression) != len(table):
            raise Exception("Boolean array and table must be the same length")
        else:
            return table[expression]
    else:
        if not isinstance(expression, str):
            raise TypeError("Expression must be a boolean array or a string")
        else:
            if len(expression) == 0:
                return table

    # parse the expression string and split into tokens
    reg = re.compile(LOGEXPRS)
    tokens = reg.split(expression)
    tokens = [t.strip() for t in tokens if t.strip() != ""]

    if isinstance(table, Table):
        # convert astropy table to pandas DataFrame
        tab = table.to_pandas()
    elif not isinstance(table, DataFrame):
        raise TypeError("Table must be a pandas DataFrame or astropy Table")
    else:
        tab = table

    matchTypes = ["ASSOC", "TYPE", "BINCOMP", "EXIST", "ERROR", "SURVEY", "DISCOVERY"]

    # parse through tokens and replace as required
    ntokens = len(tokens)
    newtokens = []
    i = 0
    while i < ntokens:
        if tokens[i] in [r"&&", r"AND", r"and"]:
            # replace synonyms for '&' or 'and'
            newtokens.append(r"&")
        elif tokens[i] in [r"||", r"OR", r"or"]:
            # replace synonyms for '|' or 'or'
            newtokens.append(r"|")
        elif tokens[i] in [r"!", r"NOT", r"not"]:
            # replace synonyms for '~'
            newtokens.append(r"~")
        elif tokens[i].upper() in matchTypes:
            if ntokens < i + 3:
                warnings.warn(
                    "A '{}' must be followed by a '(NAME)': "
                    "ignoring in query".format(tokens[i].upper()),
                    UserWarning,
                )
            elif tokens[i + 1] != "(" or tokens[i + 3] != ")":
                warnings.warn(
                    "A '{}' must be followed by a '(NAME)': "
                    "ignoring in query".format(tokens[i].upper()),
                    UserWarning,
                )
            else:
                if tokens[i].upper() == "ASSOC":
                    if "ASSOC" not in tab.columns:
                        warnings.warn(
                            "'ASSOC' parameter not in table: " "ignoring in query",
                            UserWarning,
                        )
                    elif exactMatch:
                        newtokens.append(
                            r'(ASSOC == "{}")'.format(tokens[i + 2].upper())
                        )
                    else:
                        assoc = np.array(
                            [tokens[i + 2] in str(a) for a in table["ASSOC"]]
                        )
                        newtokens.append(r"(@assoc)")
                        i += 1
                elif tokens[i].upper() == "TYPE":
                    if tokens[i + 2].upper() == "BINARY":
                        if "BINARY" not in tab.keys():
                            warnings.warn(
                                "'BINARY' parameter not in table: " "ignoring in query",
                                UserWarning,
                            )
                        else:
                            binary = ~tab["BINARY"].isna()
                            newtokens.append(r"(@binary)")
                            i += 1
                    else:
                        if "TYPE" not in tab.keys():
                            warnings.warn(
                                "'TYPE' parameter not in table: ignoring in query",
                                UserWarning,
                            )
                        elif exactMatch:
                            newtokens.append(
                                r'(TYPE == "{}")'.format(tokens[i + 2].upper())
                            )
                        else:
                            ttype = np.array(
                                [tokens[i + 2] in str(a) for a in table["TYPE"]]
                            )
                            newtokens.append(r"(@ttype)")
                            i += 1
                elif tokens[i].upper() == "BINCOMP":
                    if "BINCOMP" not in tab.columns:
                        warnings.warn(
                            "'BINCOMP' parameter not in table: ignoring in query",
                            UserWarning,
                        )
                    elif exactMatch:
                        newtokens.append(
                            r'(BINCOMP == "{}")'.format(tokens[i + 2].upper())
                        )
                    else:
                        bincomp = np.array(
                            [tokens[i + 2] in str(a) for a in table["BINCOMP"]]
                        )
                        newtokens.append(r"(@bincomp)")
                        i += 1
                elif tokens[i].upper() == "SURVEY":
                    if "SURVEY" not in tab.columns:
                        warnings.warn(
                            "'SURVEY' parameter not in table: ignoring in query",
                            UserWarning,
                        )
                    elif exactMatch:
                        newtokens.append(
                            r'(SURVEY == "{}")'.format(tokens[i + 2].upper())
                        )
                    else:
                        survey = np.array(
                            [tokens[i + 2] in str(a) for a in table["SURVEY"]]
                        )
                        newtokens.append(r"(@survey)")
                        i += 1
                elif tokens[i].upper() == "DISCOVERY":
                    if "SURVEY" not in tab.columns:
                        warnings.warn(
                            "'SURVEY' parameter not in table: ignoring in query",
                            UserWarning,
                        )
                    else:
                        if exactMatch:
                            discovery = np.array(
                                [
                                    tokens[i + 2] == str(a).split(",")[0]
                                    for a in table["SURVEY"]
                                ]
                            )
                        else:
                            discovery = np.array(
                                [
                                    tokens[i + 2] in str(a).split(",")[0]
                                    for a in table["SURVEY"]
                                ]
                            )
                        newtokens.append(r"(@discovery)")
                        i += 1
                elif tokens[i].upper() == "EXIST":
                    if tokens[i + 2].upper() not in tab.columns:
                        warnings.warn(
                            "'{}' does not exist for any pulsar".format(tokens[i + 2]),
                            UserWarning,
                        )
                        # create an empty DataFrame
                        tab = DataFrame(columns=tab.columns)
                        break
                    else:
                        exists = ~tab[tokens[i + 2].upper()].isna()
                        newtokens.append(r"(@exists)")
                        i += 1
                elif tokens[i].upper() == "ERROR":
                    if tokens[i + 2].upper() + "_ERR" not in tab.columns:
                        warnings.warn(
                            "Error value for '{}' not present: "
                            "ignoring in query".format(tokens[i + 2]),
                            UserWarning,
                        )
                    else:
                        newtokens.append(r"{}_ERR".format(tokens[i + 2].upper()))
            i += 2
        else:
            newtokens.append(tokens[i].upper())

        i += 1

    # evaluate the expression
    try:
        newtab = tab.query("".join(newtokens))
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


def characteristic_age(period, pdot, braking_idx=3.0):
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

    # convert period and pdot to numpy arrays
    periodarr = np.array(period).flatten()
    pdotarr = np.array(pdot).flatten()

    assert periodarr.dtype == float, "Periods must be floats"
    assert pdotarr.dtype == float, "Period derivatives must be floats"
    assert len(periodarr) == len(
        pdotarr
    ), "Period and derivative arrays must be equal lengths"
    assert braking_idx > 1.0, "Braking index must be greater than 1"

    # check everything is positive, otherwise return NaN
    age = np.full(len(periodarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(pdotarr) & np.isfinite(periodarr) & (pdotarr > 0.0)

    age[idx] = (periodarr[idx] / (pdotarr[idx] * (braking_idx - 1.0))) / (
        365.25 * 86400.0
    )

    # if period and period derivates were just floats return values rather
    # than arrays
    if isinstance(pdot, float) and isinstance(period, float):
        return age[0]

    return age


def age_pdot(period, tau=1e6, braking_idx=3.0):
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

    taus = tau * 365.25 * 86400.0  # tau in seconds

    pdots = periods / (taus * (braking_idx - 1.0))
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

    # convert period and pdot to numpy arrays
    periodarr = np.array(period).flatten()
    pdotarr = np.array(pdot).flatten()

    assert periodarr.dtype == float, "Periods must be floats"
    assert pdotarr.dtype == float, "Period derivatives must be floats"
    assert len(periodarr) == len(
        pdotarr
    ), "Period and derivative arrays must be equal lengths"

    # check pdot is positive, otherwise return NaN
    bfield = np.full(len(periodarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(pdotarr) & np.isfinite(periodarr) & (pdotarr > 0.0)
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

    pdots = (Bfield / 3.2e19) ** 2 / periods
    pdots[pdots < 0] = np.nan  # set any non zero values to NaN

    return pdots


def pdot_to_fdot(pdot, period=None, frequency=None):
    """
    Convert period derivatives to frequency derivatives. Either periods or
    frequencies must be supplied.

    .. math::

       \\dot{f} = -\\frac{\\dot{P}}{P^2}

    Args:
        pdot (list, :class:`~numpy.ndarray`): a list of period derivative
            values.
        period (list, :class:`~numpy.ndarray`): a list of periods (seconds)
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz)

    Returns:
        :class:`numpy.ndarray`: an array of frequency derivatives.
    """

    if period is None and frequency is None:
        raise ValueError("Either periods or frequencies must be provided")

    # convert to numpy arrays
    if period is not None:
        periodarr = np.array(period).flatten()
        pdotarr = np.array(pdot).flatten()
    else:
        periodarr = 1.0 / np.array(frequency).flatten()
        pdotarr = np.array(pdot).flatten()

    assert periodarr.dtype == float, "Periods must be floats"
    assert pdotarr.dtype == float, "Period derivatives must be floats"
    assert len(periodarr) == len(
        pdotarr
    ), "Period and period derivative arrays must be equal lengths"

    fdot = np.full(len(periodarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(pdotarr) & np.isfinite(periodarr)
    fdot[idx] = -pdotarr[idx] / periodarr[idx] ** 2

    # period derivates were just a float return values rather than arrays
    if isinstance(pdot, float):
        return fdot[0]

    return fdot


def fdot_to_pdot(fdot, period=None, frequency=None):
    """
    Convert frequency derivatives to period derivatives. Either periods or
    frequencies must be supplied.

    .. math::

       \\dot{P} = -\\frac{\\dot{f}}{f^2}

    Args:
        fdot (list, :class:`~numpy.ndarray`): a list of frequency derivative
            values.
        period (list, :class:`~numpy.ndarray`): a list of periods (seconds)
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz)

    Returns:
        :class:`numpy.ndarray`: an array of period derivatives.
    """

    if period is None and frequency is None:
        raise ValueError("Either periods or frequencies must be provided")

    # convert to numpy arrays
    if period is not None:
        frequencyarr = 1.0 / np.array(period).flatten()
        fdotarr = np.array(fdot).flatten()
    else:
        frequencyarr = np.array(frequency).flatten()
        fdotarr = np.array(fdot).flatten()

    assert frequencyarr.dtype == float, "Frequencies must be floats"
    assert fdotarr.dtype == float, "Frequency derivatives must be floats"
    assert len(frequencyarr) == len(
        fdotarr
    ), "Frequency and frequency derivative arrays must be equal lengths"

    pdot = np.full(len(frequencyarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(fdotarr) & np.isfinite(frequencyarr)
    pdot[idx] = -fdotarr[idx] / frequencyarr[idx] ** 2

    # period derivates were just a float return values rather than arrays
    if isinstance(fdot, float):
        return pdot[0]

    return pdot


def gw_h0_spindown_limit(frequency, fdot, distance, Izz=1e38):
    """
    The gravitational-wave amplutitude at Earth assuming all rotational energy
    loss is radiated via emission from an :math:`l=m=2` mass quadrupole mode,
    using (see Equation A7 of [LVC]_)

    .. math::

        h_0^{\\text{sd}} = \\frac{1}{d}\\left(\\frac{5}{2} \\frac{G I_{zz}}{c^3} \\frac{|\\dot{f}|}{f}\\right)^{1/2}

    Args:
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz).
        fdot (list, :class:`~numpy.ndarray`): a list of frequency derivatives
            (Hz/s).
        distance (list, :class:`~numpy.ndarray`): a list of distances (kpc).
        Izz (float): the moment of inertia along the principle axis (defaults
            to 1e38 kg m\\ :sup:`2`).

    Returns:
        :class:`numpy.ndarray`: an array of spin-down limits.

    .. [LVC] Abbott et al, *ApJ*, **879**, 28 (2019),
        `arXiv:1902.08507 <https://arxiv.org/abs/1902.08507>`_
    """

    from astropy.constants import G, c

    # convert values to numpy arrays
    frequencyarr = np.array(frequency).flatten()
    fdotarr = np.array(fdot).flatten()
    distancearr = (np.array(distance).flatten() * aunits.pc * 1e3).to("m").value

    assert frequencyarr.dtype == float, "Frequencies must be floats"
    assert fdotarr.dtype == float, "Frequency derivatives must be floats"
    assert len(frequencyarr) == len(fdotarr) and len(frequencyarr) == len(
        distancearr
    ), "Input arrays must be equal lengths"

    h0ul = np.full(len(frequencyarr), np.nan)
    # only use negative fdots
    with np.errstate(invalid="ignore"):
        idx = (
            np.isfinite(fdotarr)
            & np.isfinite(frequencyarr)
            & np.isfinite(distancearr)
            & (fdotarr < 0.0)
        )
    h0ul[idx] = (
        np.sqrt(
            (5.0 * G.value * Izz * np.fabs(fdotarr[idx]))
            / (2.0 * c.value ** 3 * frequencyarr[idx])
        )
        / distancearr[idx]
    )

    if isinstance(frequency, float):
        return h0ul[0]

    return h0ul


def gw_luminosity(h0, frequency, distance):
    """
    The gravitational-wave luminosity as derived from a given gravitatonal-wave
    amplitude measured at Earth, assuming all rotational energy loss is
    radiated via emission from an :math:`l=m=2` mass quadrupole mode, using
    (see Equation A5 of [LVC]_)

    .. math::

        \\dot{E}_{\\text{gw}} = \\frac{8\\pi^2c^3}{5G} f^2 h_0^2 d^2

    Args:
        h0 (list, :class:`~numpy.ndarray`): a list of gravitational-wave
            amplitudes.
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz).
        distance (list, :class:`~numpy.ndarray`): a list of distances (kpc).

    Returns:
        :class:`numpy.ndarray`: an array of luminosities.
    """

    from astropy.constants import G, c

    # convert values to numpy arrays
    frequencyarr = np.array(frequency).flatten()
    h0arr = np.array(h0).flatten()
    distancearr = (np.array(distance).flatten() * aunits.pc * 1e3).to("m").value

    assert frequencyarr.dtype == float, "Frequencies must be floats"
    assert h0arr.dtype == float, "h0 must be floats"
    assert len(frequencyarr) == len(h0arr) and len(frequencyarr) == len(
        distancearr
    ), "Input arrays must be equal lengths"

    edot = np.full(len(frequencyarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(h0arr) & np.isfinite(frequencyarr) & np.isfinite(distancearr)
    edot[idx] = (
        8.0
        * np.pi ** 2
        * c.value ** 3
        * (frequencyarr[idx] * h0arr[idx] * distancearr[idx]) ** 2
        / (5.0 * G.value)
    )

    if isinstance(frequency, float):
        return edot[0]

    return edot


def h0_to_q22(h0, frequency, distance):
    """
    Convert a gravitational-wave amplitude :math:`h_0` at Earth to an
    equivalent :math:`l=m=2` mass quadrupole :math:`Q_{22}` of the pulsar
    using (see Equations A3 and A4 of [LVC]_):

    .. math::

        Q_{22} = h_0 \\left(\\frac{c^4 d}{16\\pi^2 G f^2} \\right) \\sqrt{\\frac{15}{8\\pi}}

    Args:
        h0 (list, :class:`~numpy.ndarray`): a list of gravitational-wave
            amplitudes.
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz).
        distance (list, :class:`~numpy.ndarray`): a list of distances (kpc).

    Returns:
        :class:`numpy.ndarray`: an array of :math:`Q_{22}` values.
    """

    from astropy.constants import G, c

    # convert values to numpy arrays
    frequencyarr = np.array(frequency).flatten()
    h0arr = np.array(h0).flatten()
    distancearr = (np.array(distance).flatten() * aunits.pc * 1e3).to("m").value

    assert frequencyarr.dtype == float, "Frequencies must be floats"
    assert h0arr.dtype == float, "h0 must be floats"
    assert len(frequencyarr) == len(h0arr) and len(frequencyarr) == len(
        distancearr
    ), "Input arrays must be equal lengths"

    q22 = np.full(len(frequencyarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(h0arr) & np.isfinite(frequencyarr) & np.isfinite(distancearr)
    q22[idx] = (
        h0arr[idx]
        * np.sqrt(15.0 / (8.0 * np.pi))
        * (
            c.value ** 4
            * distancearr[idx]
            / (16.0 * np.pi ** 2 * G.value * frequencyarr[idx] ** 2)
        )
    )

    if isinstance(frequency, float):
        return q22[0]

    return q22


def q22_to_ellipticity(q22, Izz=1e38):
    """
    Convert a :math:`l=m=2` mass quadrupole :math:`Q_{22}` to an
    equivalent fiducial ellipticity :math:`\\varepsilon` of the pulsar
    using (see Equation A3 of [LVC]_):

    .. math::

        \\varepsilon = \\frac{Q_{22}}{I_{zz}}\\sqrt{\\frac{8\\pi}{15}}

    Args:
        q22 (list, :class:`~numpy.ndarray`): a list of mass quadrupole values
            (kg m\\ :sup:`2`).
        Izz (float): the moment of inertia of the principle axis (defaults to
            1e38 kg m\\ :sup:`2`).

    Returns:
        :class:`numpy.ndarray`: an array of :math:`\\varepsilon` values.
    """

    q22arr = np.array(q22).flatten()

    ellipticity = np.full(len(q22arr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(q22arr)

    ellipticity[idx] = q22arr[idx] * np.sqrt(8.0 * np.pi / 15.0) / Izz

    if isinstance(q22, float):
        return ellipticity[0]

    return ellipticity


def ellipticity_to_q22(ellipticity, Izz=1e38):
    """
    Convert a fiducial ellipticity :math:`\\varepsilon` of the pulsar to the
    equivalent :math:`l=m=2` mass quadrupole :math:`Q_{22}` using (see
    Equation A3 of [LVC]_):

    .. math::

        Q_{22} = \\varepsilon I_{zz} \\sqrt{\\frac{15}{8\\pi}}

    Args:
        ellipticity (list, :class:`~numpy.ndarray`): a list of ellipticity
            values.
        Izz (float): the moment of inertia of the principle axis (defaults to
            1e38 kg m\\ :sup:`2`).

    Returns:
        :class:`numpy.ndarray`: an array of :math:`Q_{22}` values.
    """

    ellipticityarr = np.array(ellipticity).flatten()

    q22 = np.full(len(ellipticityarr), np.nan)
    with np.errstate(invalid="ignore"):
        idx = np.isfinite(ellipticityarr)

    q22[idx] = ellipticityarr[idx] * np.sqrt(15.0 / (8.0 * np.pi)) * Izz

    if isinstance(ellipticity, float):
        return q22[0]

    return q22


def h0_to_ellipticity(h0, frequency, distance, Izz=1e38):
    """
    Convert a gravitational-wave amplitude :math:`h_0` at Earth to an
    equivalent fiducial ellipticity :math:`\\varepsilon` of the pulsar
    using (see Equation A4 of [LVC]_):

    .. math::

        \\varepsilon = h_0 \\left(\\frac{c^4 d}{16\\pi^2 G f^2} \\right) \\sqrt{\\frac{15}{8\\pi}}

    Args:
        h0 (list, :class:`~numpy.ndarray`): a list of gravitational-wave
            amplitudes.
        frequency (list, :class:`~numpy.ndarray`): a list of frequencies (Hz).
        distance (list, :class:`~numpy.ndarray`): a list of distances (kpc).
        Izz (float): the moment of inertia of the principle axis (defaults to
            1e38 kg m\\ :sup:`2`).

    Returns:
        :class:`numpy.ndarray`: an array of :math:`Q_{22}` values.
    """

    # use q22
    q22 = h0_to_q22(h0, frequency, distance)
    ellipticity = q22_to_ellipticity(q22)

    return ellipticity


def death_line(logP, linemodel="Ip", rho6=1.0):
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

    gradvals = {
        "I": (11.0 / 4),
        "Ip": (9.0 / 4.0),
        "II": (2.0 / 11.0),
        "IIp": -(2.0 / 11.0),
        "III": (5.0 / 2.0),
        "IIIp": 2.0,
        "IV": -(3.0 / 11.0),
        "IVp": -(7.0 / 11.0),
    }
    intercept = {
        "I": 14.62,
        "Ip": 16.58,
        "II": 13.07,
        "IIp": 14.50,
        "III": 14.56,
        "IIIp": 16.52,
        "IV": 15.36,
        "IVp": 16.79,
    }
    rho = {
        "I": 0.0,
        "Ip": 1.0,
        "II": 0.0,
        "IIp": (8.0 / 11.0),
        "III": 0.0,
        "IIIp": 1.0,
        "IV": 0.0,
        "IVp": (8.0 / 11.0),
    }

    lp = logP
    if not isinstance(lp, np.ndarray):
        lp = np.array(lp)

    return (
        lp * gradvals[linemodel]
        - intercept[linemodel]
        + rho[linemodel] * np.log10(rho6)
    )


def label_line(ax, line, label, color="k", fs=14, frachoffset=0.1):
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
        if ax.get_xscale() == "log":
            xx = np.log10(x1) + frachoffset * (np.log10(x2) - np.log10(x1))
        else:
            xx = x1 + frachoffset * (x2 - x1)
    else:
        raise ValueError("frachoffset must be between 0 and 1")

    if ax.get_xscale() == "log" and ax.get_yscale() == "log":
        yy = np.interp(xx, np.log10(xdata), np.log10(ydata))
        xx = 10 ** xx
        yy = 10 ** yy
    elif ax.get_xscale() == "log" and ax.get_yscale() != "log":
        yy = np.interp(xx, np.log10(xdata), ydata)
        xx = 10 ** xx
    else:
        yy = np.interp(xx, xdata, ydata)

    ylim = ax.get_ylim()
    xytext = (0, 10)

    text = ax.annotate(
        label,
        xy=(xx, yy),
        xytext=xytext,
        textcoords="offset points",
        size=fs,
        color=color,
        zorder=1,
        horizontalalignment="left",
        verticalalignment="center_baseline",
    )

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = sp2[1] - sp1[1]
    run = sp2[0] - sp1[0]

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation_mode("anchor")
    text.set_rotation(slope_degrees)
    ax.set_ylim(ylim)
    return text

