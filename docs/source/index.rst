.. psrqpy documentation master file, created by
   sphinx-quickstart on Thu Nov 23 21:34:43 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _reference:

The psrqpy package
==================

.. automodule:: psrqpy

This package provides a way to directly query the `ATNF Pulsar Catalogue <http://www.atnf.csiro.au/people/pulsar/psrcat/>`_ [1]_ using python.

Installation
============

This package can be installed using ``pip`` via ``pip install psrqpy``. Alternatively
the source code can be obtained from `here <https://github.com/mattpitkin/psrqpy>`_, and installed using::

    python setup.py install

with ``sudo`` if wanted to install system wide, and with the ``--user`` flag
if just installing for an individual user.


Contents:

.. toctree::
   :maxdepth: 1

   query
   test

.. include:: references.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

