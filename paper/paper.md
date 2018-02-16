---
title: 'psrqpy: a python interface for querying the ATNF pulsar catalogue'
tags:
  - Python
  - pulsars
authors:
  - name: Matthew Pitkin
    orcid: 0000-0003-4548-526X
    affiliation: 1
affiliations:
  - name: Institute for Gravitational Research, SUPA, University of Glasgow, University Avenue, Glasgow, UK, G12 8QQ
    index: 1
date: 16 January 2018
bibliography: paper.bib
---

# Summary

This Python module provides an interface for querying the [Australia Telescope
National Facility (ATNF) pulsar catalogue](http://www.atnf.csiro.au/people/pulsar/psrcat/) [@ATNF].
The intended users are astronomers wanting to extract data from the catalogue through a
script rather than having to download and parse text tables output using the standard web interface. It allows users to access
information, such as pulsar frequencies and sky locations, on all pulsars in
the catalogue. Querying of the catalogue can easily be incorporated into Python scripts.

The module can also be used to create plots of pulsar period against period
derivative ($P$ vs.\ $\dot{P}$ plots) using `matplotlib` [@matplotlib] as shown
below.

![A plot of pulsar period vs.\ period derivative as produced using *psrqpy*](ppdot.png)

If requested the module can also return references for parameter values for
pulsars using the `ads` Python module [@ADS].

Development of *psrqpy* happens on Github [@psrqpy_github] and the documentation
is provided [here](http://psrqpy.readthedocs.io).
