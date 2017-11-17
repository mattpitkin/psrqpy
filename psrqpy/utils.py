"""
Various useful functions
"""

import warnings
import requests
import re
import datetime
from bs4 import BeautifulSoup

from .config import *

# problematic references that are hard to parse
prob_refs = ['bwck08']

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
        except:
            warnings.warn("Could not get ATNF version number, defaulting to {}".format(ATNF_VERSION), UserWarning)
            atnfversion = ATNF_VERSION

    return atnfversion

def get_references(useads=False):
    """
    Return a dictionary of paper reference in the ATNF catalogue. The keys are the ref strings given
    in the ATNF catalogue.

    :param useads: boolean to set whether to use the python 'ads' module to get the NASA ADS URL for the references
    """

    refs = {}

    queryrefs = requests.get(ATNF_BASE_URL + 'psrcat_ref.html')

    if queryrefs.status_code != 200:
        warnings.warn("Could query the ATNF references. No references returned".format(ATNF_VERSION), UserWarning)
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

            if reftag in prob_refs:
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

                try:
                    article = list(ads.SearchQuery(year=refs[reftag]['year'], first_author=refs[reftag]['authors'][0], title=refs[reftag]['title']))[0]
                    refs[reftag]['ADS'] = article
                    refs[reftag]['ADS URL'] = ADS_URL.format(article.bibcode)
                except IOError:
                    warnings.warn('Could not import ADS module, so no ADS information will be included', UserWarning)

    return refs
