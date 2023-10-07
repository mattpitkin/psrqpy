# coding: utf-8

""" A Python tool for interacting with the ATNF pulsar catalogue """

import warnings
from .search import QueryATNF
from .pulsar import Pulsar, Pulsars
from .utils import *

try:
    from ._version import version as __version__
except ModuleNotFoundError:
    __version__ = ""


__citation__ = """@article{psrqpy,
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
"""


# set formatting of warnings to not include line number and code (see
# e.g. https://pymotw.com/3/warnings/#formatting)
def warning_format(message, category, filename, lineno, file=None, line=None):
    return '{}: {}\n'.format(category.__name__, message)


warnings.formatwarning = warning_format
