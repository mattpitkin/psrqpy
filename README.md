# PSRQpy

This module aims to provide a python interface for querying the [ATNF pulsar catalogue](http://www.atnf.csiro.au/people/pulsar/psrcat/).
It is inspired by, and has some minor similarities to, the [`ads`](https://ads.readthedocs.io) module for interfacing with the
[NASA ADS](https://ui.adsabs.harvard.edu/) [API](https://github.com/adsabs/adsabs-dev-api).

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

### Requirements

The requirements for installing the code are:

 * `six`
 * `requests`
 * `beautifulsoup4`
 * `numpy`
 * `astropy`
 * `datetime`

## License

This code is licensed under the [MIT License](http://opensource.org/licenses/MIT).

&copy; Matt Pitkin, 2017