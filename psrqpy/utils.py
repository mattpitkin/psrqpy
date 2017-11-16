"""
Various useful functions
"""

import warnings
import requests
import re
from bs4 import BeautifulSoup

from .config import *

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
        except:
            warnings.warn("Could not get ATNF reference list", UserWarning)
            return refs

        # loop over rows
        for tr in trows:
            reftag = tr.b.text # the reference string is contained in a <b> tag

            refs[reftag] = {}
            tds = tr.find_all('td') # get the two <td> tags - reference info is in the second

            refdata = tds[1].contents

            # check that the tag contains a string (the paper title) within <i> - there are a few exceptions to this rule
            if tds[1].find('i') is not None:
                authors = re.sub(r'\s+', ' ', refdata[0]).strip().strip('.') # remove line breaks and extra spaces (and final full-stop)
            else:
                authors = re.sub(r'\s+', ' ', refdata[0]).split('.,')[:-2]
                # authors, year, and journal info are all together
                rd = re.sub(r'\s+', ' ', refdata[0]).split('.,')[-1].split() # split on whitespace
                
            sepauthors = [a for a in authors.split('.,')]
                
            if tds[1].find('i') is not None:
                year = int(''.join(filter(lambda x: x.isdigit(), sepauthors.pop(-1).strip('.')))) # strip any non-digit characters (e.g. from '1976a')
            else:
                rd = re.sub(r'\s+', ' ', refdata[0]).split('.,')[-1].split() # split on whitespace
                year = int(''.join(filter(lambda x: x.isdigit(), rd[0].strip('.'))))

            if '&' in sepauthors[-1]: # split any authors that are seperated by an ampersand
                lastauthors = [a.strip() for a in sepauthors.pop(-1).split('&')]
                sepauthors = sepauthors + lastauthors
                for i in xrange(len(sepauthors)-2):
                    sepauthors[i] += '.' # re-add final full stops where needed
            else:
                sepauthors = [a+'.' for a in sepauthors] # re-add final full stops

            refs[reftag]['authorlist'] = ', '.join(sepauthors)
            refs[reftag]['authors'] = sepauthors
            refs[reftag]['year'] = year

            if tds[1].find('i') is not None:
                title = re.sub(r'\s+', ' ', tds[1].i.text)
            else:
                title = ''
            refs[reftag]['title'] = title

            if tds[1].find('i') is not None:
                # seperate journal name, volume and pages
                journalref = [a.strip() for a in refdata[-1].strip('.').split(',')]
                if len(journalref) == 3:
                    refs[reftag]['journal'] = journalref[0]
                    refs[reftag]['volume'] = journalref[1]
                    refs[reftag]['pages'] = journalref[2]
                else:
                    if 'arxiv' in refdata[-1].strip('.').lower():
                        axvparts = refdata[-1].strip('.').split(':')
                        axv = 'arXiv:{}:{}'.format(axvparts[0].split()[-1], axvparts[1].split()[0])
                        refs[reftag]['journal'] = axv
                    else:
                        refs[reftag]['journal'] = ''
                    refs[reftag]['volume'] = ''
                    refs[reftag]['pages'] = ''
            else:
                refs[reftag]['year'] = year
                refs[reftag]['journal'] = rd[1].strip(',')
                refs[reftag]['volume'] = rd[2].strip(',')
                refs[reftag]['pages'] = rd[3].strip('.')

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
                except:
                    warnings.warn('Could not import ADS module, so no ADS information will be included', UserWarning)

    return refs
