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

Examples
========

You can query the ATNF catalogue for any combination of the pulsar parameters listed
`here <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#par_list>`_.

A simple example of a query is to get the 'JName' (pulsar name based on the J2000 coordinates) and frequency 'F0' for all pulsars in the catalogue. This could be done with

   >>> from psrqpy import QueryATNF
   >>> query = QueryATNF(params=['JName', 'F0'])

where the parameter names are case insensitive. The returned :class:`psrqpy.QueryATNF` query will
contain a dictionary keyed on the parameter names (converted
to upper case), with the
dictionary values being :class:`numpy.ndarray` arrays containing the parameter values. This can
be accessed with

   >>> qdict = query.get_dict()

or, perhaps more conveniently, the output can be viewed as an :class:`astropy.table.Table` via

   >>> qtable = query.table()

The number of pulsars can easily be accessed, e.g.,

   >>> numstring = 'Version {} of the ATNF catalogue contains {} pulsars'
   >>> print(numstring.format(query.get_version, query.num_pulsars))
   Version 1.57 of the ATNF catalogue contains 2627 pulsars

The code will automatically attempt to query the current version of the
ATNF catalogue, the value of which is printed in the above example.

More complex queries
--------------------



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

