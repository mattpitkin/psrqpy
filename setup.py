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

MAJOR, MINOR1, MINOR2, RELEASE, SERIAL = sys.version_info

READFILE_KWARGS = {"encoding": "utf-8"} if MAJOR >= 3 else {}

def readfile(filename):
    with open(filename, **READFILE_KWARGS) as fp:
        filecontents = fp.read()
    return filecontents

VERSION_REGEX = re.compile("__version__ = \"(.*?)\"")
CONTENTS = readfile(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "psrqpy", "__init__.py"))

VERSION = VERSION_REGEX.findall(CONTENTS)[0]

# 'setup.py publish' shortcut for publishing (e.g. setup.py from requests https://github.com/requests/requests/blob/master/setup.py)
# 'setup.py publish' shortcut.
if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist bdist_wheel --universal')
    os.system('twine upload dist/*')
    sys.exit()

if sys.argv[-1] == 'testpublish':
    os.system('python setup.py sdist bdist_wheel --universal')
    os.system('twine upload --repository-url https://test.pypi.org/legacy/ dist/*')
    sys.exit()

setup(name="psrqpy",
      version=VERSION,
      author="Matthew Pitkin",
      author_email="matthew.pitkin@glasgow.ac.uk",
      packages=["psrqpy"],
      url="http://www.github.com/mattpitkin/psrqpy/",
      license="MIT",
      description="A Python module for querying the ATNF pulsar catalogue",
      long_description=\
          readfile(os.path.join(os.path.dirname(__file__), "README.md")),
      install_requires=\
          readfile(os.path.join(os.path.dirname(__file__), "requirements.txt")),
      classifiers=[
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7"])
