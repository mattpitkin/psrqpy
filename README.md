# PSRQpy

This module aims to provide a python interface for querying the [ATNF pulsar catalogue](http://www.atnf.csiro.au/people/pulsar/psrcat/).
It is inspired by, and has some minor similarities to, the [`ads`](https://ads.readthedocs.io) module for interfacing with the
[NASA ADS](https://ui.adsabs.harvard.edu/) [API](https://github.com/adsabs/adsabs-dev-api). It is an unofficial
package and is not endorsed by or affiliated with the ATNF.

Full documentation of the module can be found [here](http://psrqpy.readthedocs.io/).

Any comments or suggestions are welcome.

## Installation

To install the code from source, clone the git repository and run either:

```
python setup.py install --user
```

to install as a user, or

```
sudo python setup.py install
```

to install as root.

The module can also be installed using `pip` with:

```
pip install psrqpy
```

### Requirements

The [requirements](requirements.txt) for installing the code are:

 * `six`
 * `requests`
 * `beautifulsoup4`
 * `numpy`
 * `astropy` (for Python 2 astropy versions before [3.0](http://docs.astropy.org/en/latest/whatsnew/3.0.html#whatsnew-3-0-python3) must be used)
 * `datetime`

The `ads` module is an optional requirement that is needed to get ADS URLs for references.

## Examples

A simple query of the catalogue, e.g., to just return all pulsar frequencies, would be:

```python
import psrqpy

q = psrqpy.QueryATNF(params='F0')

# get frequencies as an astropy table
t = q.table()

print t['F0']
```

You can query multiple parameters, e.g.:

```python
import psrqpy

q = psrqpy.QueryATNF(params=['F0', 'F1', 'RAJ', 'DecJ'])

# get values as an astropy table
t = q.table()

print t['F0']
```

You can query specific pulsars, e.g.:

```
import psrqpy

q = psrqpy.QueryATNF(params=['JName', 'F0', 'F1', 'RAJ', 'DecJ'], psrs=['J0534+2200', 'J0537-6910'])

# get values as an astropy table
t = q.table()

# print the table
print(t)
  JNAME          F0       F0_ERR       F1      F1_ERR     RAJ      RAJ_ERR     DECJ     DECJ_ERR
                 Hz         Hz       1 / s2    1 / s2                                           
---------- ------------- ------- ------------- ------ ------------ ------- ------------ --------
J0534+2200     29.946923   1e-06  -3.77535e-10  2e-15 05:34:31.973   0.005 +22:00:52.06     0.06
J0537-6910 62.0261895958 1.3e-09 -1.992272e-10  4e-17 05:37:47.416    0.11 -69:10:19.88      0.6
```

You can set [conditions](http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#condition) for the searches,
e.g.:

```python
import psrqpy
q = psrqpy.QueryATNF(params=['Jname', 'f0'], condition='f0 > 100 && f0 < 200', assoc='GC')
```

where `assoc=GC` looks for all pulsars in globular clusters.

If you really want to query the catalogue many times in quick succession it is probably preferable to [download
the catalogue](http://www.atnf.csiro.au/research/pulsar/psrcat/download.html) and query it with the software
provided. The enitre catalogue can be downloaded using:

```python
import psrqpy
catalogue = psrqpy.get_catalogue()
```

You can also [generate](http://psrqpy.readthedocs.io/en/latest/query.html#psrqpy.search.QueryATNF.ppdot) a
_lovely_ period vs. period derivative plot based on the latest catalogue information, using
just three lines of code, e.g.:

```python
from psrqpy import QueryATNF
query = QueryATNF(params=['P0', 'P1', 'ASSOC', 'BINARY', 'TYPE', 'P1_I'])
query.ppdot(showSNRs=True, showtypes='all')
```

gives

![PPdot](../master/docs/source/images/ppdot.png)

## Development and Support

Code development is done via the package's [GitHib repository](https://github.com/mattpitkin/psrqpy).
Any contributions can be made via a [fork and pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) model
from that repository, and must adhere to the [MIT license](#License). Any problems with the code
or support requests can be submitted via the repository's [Issue tracker](https://github.com/mattpitkin/psrqpy/issues).

## Copyright and referencing for the catalogue

Regarding the use of the catalogue and software behind it, the [following statements](http://www.atnf.csiro.au/research/pulsar/psrcat/download.html) apply:

> PSRCAT is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. PSRCAT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
>
> PSRCAT makes use of "evaluateExpression: A Simple Expression Evaluator". Copyright &copy; 1996 - 1999 Parsifal Software, All Rights Reserved.
>
> The programs and databases remain the property of the Australia Telescope National Facility, CSIRO, and are covered by the [CSIRO Legal Notice and Disclaimer](http://www.csiro.au/en/About/Footer/Legal-notice).
>
> If you make use of information from the ATNF Pulsar Catalogue in a publication, we would appreciate acknowledgement by reference to the publication "[The ATNF Pulsar Catalogue](http://adsabs.harvard.edu/abs/2005AJ....129.1993M)", R. N. Manchester, G. B. Hobbs, A. Teoh & M. Hobbs, Astronomical Journal, 129, 1993-2006 (2005) and by quoting the web address http://www.atnf.csiro.au/research/pulsar/psrcat for updated versions.

If making use of this code to access the catalogue, or produce plots, I would be grateful if (as well as citing the ATNF pulsar catalogue [paper](http://adsabs.harvard.edu/abs/2005AJ....129.1993M) and [URL](http://www.atnf.csiro.au/research/pulsar/psrcat) given above) you consider citing the [JOSS](http://joss.theoj.org/) [paper](https://doi.org/10.21105/joss.00538) for this software:

```tex
@article{psrqpy,
  author = {{Pitkin}, M.},
   title = "{psrqpy: a python interface for querying the ATNF pulsar catalogue}",
  volume = 3,
  number = 22,
   pages = 538,
   month = feb,
    year = 2018,
 journal = "{Journal of Open Source Software}",
     doi = {10.21105/joss.00538},
     url = {https://doi.org/10.21105/joss.00538}
}
```

## License

This code is licensed under the [MIT License](http://opensource.org/licenses/MIT).

&copy; Matt Pitkin, 2017

[![PyPI version](https://badge.fury.io/py/psrqpy.svg)](https://badge.fury.io/py/psrqpy)
[![Build Status](https://travis-ci.org/mattpitkin/psrqpy.svg?branch=master)](https://travis-ci.org/mattpitkin/psrqpy)
[![Documentation Status](https://readthedocs.org/projects/psrqpy/badge/?version=latest)](http://psrqpy.readthedocs.io/en/latest/?badge=latest)
[![status](http://joss.theoj.org/papers/711dc5566159f6e9f8ea5d07dbfaf5d2/status.svg)](http://joss.theoj.org/papers/711dc5566159f6e9f8ea5d07dbfaf5d2)

