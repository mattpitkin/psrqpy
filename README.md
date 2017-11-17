# PSRQpy

This module aims to provide a python interface for querying the [ATNF pulsar catalogue](http://www.atnf.csiro.au/people/pulsar/psrcat/).
It is inspired by, and has some minor similarities to, the [`ads`](https://ads.readthedocs.io) module for interfacing with the
[NASA ADS](https://ui.adsabs.harvard.edu/) [API](https://github.com/adsabs/adsabs-dev-api). It is an unofficial
package and is not endorsed or affiliated with the ATNF.

Any comments or suggestions are welcome.

## Installation

To install the code from source clone the git repository and run either:

```
python setup.py install --user
```

to install as a user, or

```
sudo python setup.py install --user
```

to install as root.

The module can also be installed using `pip` with:

```
pip install psrqpy
```

### Requirements

The requirements for installing the code are:

 * `six`
 * `requests`
 * `beautifulsoup4`
 * `numpy`
 * `astropy`
 * `datetime`

The `ads` module is an optional requirement that is needed to get ADS URLs for references.

## Example

A simple query of the catalogue, e.g., to just return all pulsar frequencies, would be:

```
import psrqpy

q = psrqpy.QueryATNF(params='F0')

# get frequencies as an astropy table
t = q.table()

print t['F0']
```

You can query multiply parameters, e.g.:

```
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
```

If you really want to query the catalogue many times in quick succession it is advisable not to use this module, as
it could result in too much load on the ATNF catalogue's server. Instead it is probably preferable to [download
the catalogue](http://www.atnf.csiro.au/research/pulsar/psrcat/download.html) and query it with the software
provided.

## License

This code is licensed under the [MIT License](http://opensource.org/licenses/MIT).

&copy; Matt Pitkin, 2017

[![PyPI version](https://badge.fury.io/py/psrqpy.svg)](https://badge.fury.io/py/psrqpy)
