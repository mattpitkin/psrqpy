# coding: utf-8

""" A Python tool for interacting with the ATNF pulsar catalogue """

__version__ = "0.5.1"

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

from .search import QueryATNF
from .pulsar import Pulsar, Pulsars
from .utils import *
