# coding: utf-8

"""
A Python module for accessing the ATNF pulsar catalogue
"""

import os
import re
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

major, minor1, minor2, release, serial =  sys.version_info

readfile_kwargs = {"encoding": "utf-8"} if major >= 3 else {}

def readfile(filename):
    with open(filename, **readfile_kwargs) as fp:
        contents = fp.read()
    return contents

version_regex = re.compile("__version__ = \"(.*?)\"")
contents = readfile(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "psrqpy", "__init__.py"))

version = version_regex.findall(contents)[0]

setup(name="psrqpy",
      version=version,
      author="Matthew Pitkin",
      author_email="matthew.pitkin@glasgow.ac.uk",
      packages=["psrqpy"],
      url="http://www.github.com/mattpitkin/psrqpy/",
      license="MIT",
      description="A Python module for querying the ATNF pulsar catalogue",
      long_description=\
          readfile(os.path.join(os.path.dirname(__file__), "README.md")),
      install_requires=\
          readfile(os.path.join(os.path.dirname(__file__), "requirements.txt"))
     )
