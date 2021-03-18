# coding: utf-8

"""
A Python module for accessing the ATNF pulsar catalogue
"""

import os
import re

from setuptools import setup


def readfile(filename):
    with open(filename, encoding="utf-8") as fp:
        filecontents = fp.read()
    return filecontents


VERSION_REGEX = re.compile('__version__ = "(.*?)"')
CONTENTS = readfile(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "psrqpy", "__init__.py")
)

VERSION = VERSION_REGEX.findall(CONTENTS)[0]

setup(
    name="psrqpy",
    version=VERSION,
    author="Matthew Pitkin",
    author_email="matthew.pitkin@glasgow.ac.uk",
    packages=["psrqpy"],
    url="http://www.github.com/mattpitkin/psrqpy/",
    license="MIT",
    description="A Python module for querying the ATNF pulsar catalogue",
    long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
    long_description_content_type="text/markdown",
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")
    ),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)
