"""
A selection of useful functions used by the module.
"""

from __future__ import print_function, division

import warnings
import re
import datetime
import numpy as np
import requests

from bs4 import BeautifulSoup

from .config import ATNF_BASE_URL, ATNF_VERSION, ADS_URL

# problematic references that are hard to parse
PROB_REFS = ['bwck08']

def get_version():
    """
    Return a string with the ATNF version number, or default to that defined in ATNF_VERSION
    """

    site = requests.get(ATNF_BASE_URL)

    if site.status_code != 200:
        warnings.warn("Could not get ATNF version number, defaulting to {}".format(ATNF_VERSION), UserWarning)
        atnfversion = ATNF_VERSION
    else:
        # parse the site content with BeautifulSoup
        vsoup = BeautifulSoup(site.content, 'html.parser')

        try:
            vsoup = BeautifulSoup(site.content, 'html.parser')

            version = vsoup.find(attrs={'name': 'version'})
            atnfversion = version['value']
        except IOError:
            warnings.warn("Could not get ATNF version number, defaulting to {}".format(ATNF_VERSION), UserWarning)
            atnfversion = ATNF_VERSION

    return atnfversion

def get_references(useads=False):
    """
    Return a dictionary of paper
    `reference <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_ref.html>`_ in the ATNF
    catalogue. The keys are the ref strings given in the ATNF catalogue.

    Args:
        useads (bool): boolean to set whether to use the python mod:`ads` module to get
            the NASA ADS URL for the references

    Returns:
        dict: a dictionary of references.
    """

    refs = {}

    queryrefs = requests.get(ATNF_BASE_URL + 'psrcat_ref.html')

    if queryrefs.status_code != 200:
        warnings.warn("Could query the ATNF references. No references returned", UserWarning)
    else:
        try:
            refsoup = BeautifulSoup(queryrefs.content, 'html.parser')

            # get table containing the references
            pattern = re.compile('References') # References are in a h2 tag containing 'References'
            table = refsoup.find('h2', text=pattern).parent.find('table') # get the table in the same parent element as the 'References' header

            trows = table.find_all('tr')
        except IOError:
            warnings.warn("Could not get ATNF reference list", UserWarning)
            return refs

        # loop over rows
        j = 0
        for tr in trows:
            j = j + 1
            reftag = tr.b.text # the reference string is contained in a <b> tag

            if reftag in PROB_REFS:
                continue

            refs[reftag] = {}
            tds = tr.find_all('td') # get the two <td> tags - reference info is in the second

            # check if publication is 'awkward', i.e. if has a year surrounded by '.'s, e.g, '.1969.' or '.1969a.'
            utext = re.sub(r'\s+', ' ', tds[1].text)
            dotyeardot = re.compile(r'\.(\d+\D?)\.')
            dotyeardotlist = dotyeardot.split(utext)
            if len(dotyeardotlist) != 3:
                utext = None

            refdata = list(tds[1].contents) # copy list so contents of table aren't changed in the journal name substitution step below

            # check that the tag contains a string (the paper/book title) within <i> (paper) or <b> (book) - there are a few exceptions to this rule
            titlestr = None
            booktitlestr = None
            if tds[1].find('i') is not None:
                titlestr = tds[1].i.text
            if tds[1].find('b') is not None:
                booktitlestr = tds[1].b.text

            # change some journal refs that contain '.' (this causes issues when splitting authors based on '.,')
            for ridx, rdf in enumerate(list(refdata)):
                # subtitute some journal names to abbreviated versions
                journalsubs = {'Chin. J. Astron. Astrophys.': 'ChJAA',
                               'Astrophys. Lett.': 'ApJL',
                               'Res. Astron. Astrophys.': 'RAA',
                               'J. Astrophys. Astr.': 'JApA',
                               'Curr. Sci.': 'Current Science',
                               'Astrophys. Space Sci.': 'Ap&SS',
                               'Nature Phys. Sci.': 'NPhS',
                               'Sov. Astron. Lett.': 'SvAL',
                               'ATel.': 'ATel'}

                if isinstance(rdf, basestring): # only run on string values
                    rdfs = re.sub(r'\s+', ' ', rdf) # make sure only single spaces are present
                    for js in journalsubs:
                        if js in rdfs:
                            refdata[ridx] = re.sub(js, journalsubs[js], rdfs)

            if (titlestr is not None or booktitlestr is not None) and utext is None:
                authors = re.sub(r'\s+', ' ', refdata[0]).strip().strip('.') # remove line breaks and extra spaces (and final full-stop)
                sepauthors = authors.split('.,')
            elif utext is not None:
                year = int(re.sub('\D', '', dotyeardotlist[1])) # remove any non-number values
                authors = dotyeardotlist[0]
                sepauthors = authors.split('.,')
            else:
                sepauthors = re.sub(r'\s+', ' ', refdata[0]).split('.,')[:-1]

            if (titlestr is not None or booktitlestr is not None) and utext is None:
                try:
                    year = int(''.join(filter(lambda x: x.isdigit(), sepauthors.pop(-1).strip('.')))) # strip any non-digit characters (e.g. from '1976a')
                except ValueError:
                    # get year from reftag
                    year = int(''.join(filter(lambda x: x.isdigit(), reftag)))
                    thisyear = int(str(datetime.datetime.now().year)[-2:])
                    if year > thisyear:
                        year += 1900
                    else:
                        year += 2000
            elif utext is None:
                rd = re.sub(r'\s+', ' ', refdata[0]).split('.,')[-1].split()
                try:
                    year = int(''.join(filter(lambda x: x.isdigit(), rd[0].strip('.'))))
                except ValueError:
                    # get year from reftag
                    year = int(''.join(filter(lambda x: x.isdigit(), reftag)))
                    thisyear = int(str(datetime.datetime.now().year)[-2:])
                    if year > thisyear:
                        year += 1900
                    else:
                        year += 2000

            if '&' in sepauthors[-1] or 'and' in sepauthors[-1]: # split any authors that are seperated by an ampersand
                lastauthors = [a.strip() for a in re.split(r'& | and ', sepauthors.pop(-1))]
                sepauthors = sepauthors + lastauthors
                for i in xrange(len(sepauthors)-2):
                    sepauthors[i] += '.' # re-add final full stops where needed
                sepauthors[-1] += '.'
            else:
                sepauthors = [a+'.' for a in sepauthors] # re-add final full stops

            refs[reftag]['authorlist'] = ', '.join(sepauthors)
            refs[reftag]['authors'] = sepauthors
            refs[reftag]['year'] = year
            refs[reftag]['journal'] = ''
            refs[reftag]['volume'] = ''
            refs[reftag]['pages'] = ''

            if titlestr is not None:
                title = (re.sub(r'\s+', ' ', titlestr)).lstrip() # remove any leading spaces
            else:
                title = ''
            refs[reftag]['title'] = title

            if booktitlestr is not None:
                booktitle = (re.sub(r'\s+', ' ', booktitlestr)).lstrip()
                refs[reftag]['booktitle'] = booktitle

            if titlestr is not None:
                # seperate journal name, volume and pages
                journalref = [a.strip() for a in refdata[-1].strip('.').split(',')]
                if len(journalref) == 3:
                    refs[reftag]['journal'] = journalref[0]
                    refs[reftag]['volume'] = journalref[1]
                    refs[reftag]['pages'] = journalref[2]
                else:
                    if 'arxiv' in refdata[-1].strip('.').lower():
                        axvparts = refdata[-1].strip('.').split(':')
                        if len(axvparts) == 2: # if an arXiv number of found
                            axv = 'arXiv:{}'.format(re.split(', |. ', axvparts[1])[0])
                        else:
                            axv = 'arXiv' # no arXiv number can be set
                        refs[reftag]['journal'] = axv
            elif booktitlestr is not None:
                # seperate book volume and other editorial/publisher info
                bookref = [a.strip() for a in refdata[-1].strip('.').split('eds')]
                refs[reftag]['volume'] = re.sub(r', |. |\s+', '', bookref[0])
                refs[reftag]['eds'] = bookref[1]
            else:
                refs[reftag]['year'] = year

                # split on year
                if utext is None:
                    rd = re.sub(r'\s+', ' ', refdata[0]).split('{}'.format(year))[1].split(',')
                else:
                    rd = re.sub(r'\s+', ' ', dotyeardotlist[-1]).split(',')

                if 'PhD thesis' in rd[0]:
                    refs[reftag]['journal'] = 'PhD thesis'
                    refs[reftag]['thesis pub. info.'] = ' '.join(rd[1:]).lstrip()
                else:
                    if len(rd) >= 1:
                        refs[reftag]['journal'] = rd[0].strip()
                    if len(rd) >= 2:
                        refs[reftag]['volume'] = rd[1].strip()
                    if len(rd) >= 3:
                        refs[reftag]['pages'] = rd[2].strip().strip('.')

            # get ADS entry
            if useads:
                try:
                    import ads
                except ImportError:
                    warnings.warn('Could not import ADS module, so no ADS information will be included', UserWarning)
                    continue

                refs[reftag]['ADS'] = None
                refs[reftag]['ADS URL'] = ''

                try:
                    article = ads.SearchQuery(year=refs[reftag]['year'], first_author=refs[reftag]['authors'][0], title=refs[reftag]['title'])
                except IOError:
                    warnings.warn('Could not get reference information, so no ADS information will be included', UserWarning)
                    continue

                article = list(article)

                if len(article) > 0:
                    refs[reftag]['ADS'] = list(article)[0]
                    refs[reftag]['ADS URL'] = ADS_URL.format(list(article)[0].bibcode)

    return refs


def characteristic_age(period, pdot, braking_idx=3.):
    """
    Function defining the characteristic age of a pulsar. Returns the characteristic
    age in using

    .. math::
    
       \\tau = \\frac{P}{\dot{P}(n-1)}

    Args:
        period (float): the pulsar period in seconds
        pdot (float): the pulsar period derivative
        braking_idx (float): the pulsar braking index (defaults to :math:`n=3`)

    Returns:
        float: the characteristic age in years
    """

    # check everything is positive, otherwise return NaN
    if period < 0.:
        warnings.warn("The period must be positive to define a characteristic age", UserWarning)
        return np.nan

    if pdot < 0.:
        warnings.warn("The period derivative must be positive to define a characteristic age", UserWarning)
        return np.nan

    if braking_idx < 0.:
        warnings.warn("The braking index must be positive to define a characteristic age", UserWarning)
        return np.nan

    return (period/(pdot * (braking_idx - 1.)))/(365.25*86400.)


def age_pdot(period, tau=1e6, braking_idx=3.):
    """
    Function returning the period derivative for a pulsar with a given period
    and characteristic age, using

    .. math::

       \dot{P} = \\frac{\dot{P}}{\\tau(n - 1)}

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

    taus = tau*365.25*86400. # tau in seconds

    pdots = (periods/(taus * (braking_idx - 1.)))
    pdots[pdots < 0] = np.nan # set any non zero values to NaN

    return pdots


def B_field(period, pdot):
    """
    Function defining the polar magnetic field strength at the surface of the pulsar
    in gauss (Equation 5.12 of Lyne & Graham-Smith, Pulsar Astronmy, 2nd edition) with
    
    .. math::

       B = 3.2\!\\times\!10^{19} \\sqrt{P\dot{P}}

    Args:
        period (float): a pulsar period (s)
        pdot (float): a period derivative

    Returns:
        float: the magnetic field strength in gauss.
    """

    assert isinstance(period, float) or isinstance(period, int), "Period '{}' must be a number".format(period)
    assert isinstance(pdot, float) or isinstance(pdot, int), "Period derivtaive '{}' must be a number".format(pdot)

    # check everything is positive, otherwise return 0
    if period < 0.:
        warnings.warn("The period must be positive to define a magnetic field strength", UserWarning)
        return 0.

    if pdot < 0.:
        warnings.warn("The period derivative must be positive to define a magnetic field streng", UserWarning)
        return 0.

    return 3.2e19 * np.sqrt(period * pdot)


def B_field_pdot(period, Bfield=1e10):
    """
    Function to get the period derivative from a given pulsar period and magnetic
    field strength using
    
    .. math::
    
       \dot{P} = \\frac{1}{P}\left( \\frac{B}{3.2\!\\times\!10^{19}} \\right)^2

    Args:
        period (list, :class:`~numpy.ndarray`): a list of period values
        Bfield (float): the polar magnetic field strength (Defaults to :math:`10^{10}` G)
        
    Returns:
        :class:`numpy.ndarray`: an array of period derivatives
    """

    periods = period
    if not isinstance(periods, np.ndarray):
        periods = np.array(periods)

    pdots = (Bfield/3.2e19)**2/periods
    pdots[pdots < 0] = np.nan # set any non zero values to NaN

    return pdots


def death_line(logP, linemodel='Ip', rho6=1.):
    """
    The pulsar death line. Returns the base-10 logarithm of the period derivative for the given
    values of the period.

    Args:
        logP (list, :class:`~numpy.ndarray`): the base-10 log values of period.
        linemodel (str): a string with one of the above model names. Defaults to ``'Ip'``.
        rho6 (float): the value of the :math:`\\rho_6` parameter from [ZHM]_ . Defaults to 1 is,
            which is equivalent to :math:`10^6` cm.

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

    gradvals = {'I': (11./4), 'Ip': (9./4.), 'II': (2./11.), 'IIp': -(2./11.), 'III': (5./2.), 'IIIp': 2., 'IV': -(3./11.), 'IVp': -(7./11.)}
    intercept = {'I': 14.62, 'Ip': 16.58, 'II': 13.07, 'IIp': 14.50, 'III': 14.56, 'IIIp': 16.52, 'IV': 15.36, 'IVp': 16.79}
    rho = {'I': 0., 'Ip': 1., 'II': 0., 'IIp': (8./11.), 'III': 0., 'IIIp': 1., 'IV': 0., 'IVp': (8./11.)}

    lp = logP
    if not isinstance(lp, np.ndarray):
        lp = np.array(lp)

    return lp*gradvals[linemodel] - intercept[linemodel] + rho[linemodel]*np.log10(rho6)


def label_line(ax, line, label, color='k', fs=14, frachoffset=0.1):
    """
    Add an annotation to the given line with appropriate placement and rotation.

    Based on code from `"How to rotate matplotlib annotation to match a line?"
    <http://stackoverflow.com/a/18800233/230468>`_ and `this
    <https://stackoverflow.com/a/38414616/1862861>`_ answer.

    Args:
        ax (:class:`matplotlib.axes.Axes`): Axes on which the label should be added.
        line (:class:`matplotlib.lines.Line2D`): Line which is being labeled.
        label (str): Text which should be drawn as the label.
        color (str): a color string for the label text. Defaults to ``'k'``
        fs (int): the font size for the label text. Defaults to 14.
        frachoffset (float): a number between 0 and 1 giving the fractional offset of the label
            text along the x-axis. Defaults to 0.1, i.e. 10%.


    Returns:
        :class:`matplotlib.text.Text`: an object containing the label information

    """
    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    # use fractional horizontal offset frachoffset to set the x position of the label by default
    # other wise use the halign value
    if frachoffset >= 0 and frachoffset <= 1:
        if ax.get_xscale() == 'log':
            xx = np.log10(x1) + frachoffset*(np.log10(x2) - np.log10(x1))
        else:
            xx = x1 + frachoffset*(x2 - x1)
    else:
        raise ValueError("frachoffset must be between 0 and 1".format(halign))

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
