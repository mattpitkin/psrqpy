---
title: 'psrqpy: a python interface for querying the ATNF pulsar catalogue'
tags:
  - Python
  - pulsars
authors:
  - name: Matthew Pitkin
    orcid: 0000-0003-4548-526X
    affiliation: Institute for Gravitational Research, SUPA, University of Glasgow, University Avenue, UK, G12 8QQ
date: 16 January 2018
bibliography: paper.bib
---

# Summary

This Python module provides an interface for querying the Australia Telescope
National Facility (ATNF) pulsar catalogue [@ATNF]. It allows users to access
information (such as pulsar frequencies and sky locations) on all pulsars in
the catalogue without having to use the current web interface. As such
querying of the catalogue can easily be incorporated into scripts.

The module can also be used to create plots of pulsar period against period
derivative ($P$ vs.\ $\dot{P} plots) using `matplotlib` [@matplotlib]. It can
also return references for parameter values for pulsars using the `ads` Python
module [@ADS].

Development of *psrqpy* happens on Github [@psrqpy_github].

