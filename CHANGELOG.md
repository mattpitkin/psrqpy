# Notable changes between versions

## [1.0.5] 2019-09-04

Changes for this release:

- The version 1.61 of the ATNF Pulsar Catalogue containing the declinations with the unicode "−" character in place of an ascii minus sign "-" has been fixed upstream, so the code fix to account for this has been removed (it was also causing issues with the Python 2 version of the code).

## [1.0.4] 2019-09-04

Changes for this release:

- Update the P-Pdot diagram, so that the user can specify which pulsars to include, either via the pulsars passed to the original query request (which is the default option) or by directly passing names to the P-Pdot function. If no particular pulsars were requested then all pulsars will be included.
- Version 1.61 of the ATNF Pulsar Catalogue contains some declinations with the unicode "−" character in place of an ascii minus sign "-". This release switched the unicode "−" to the "-", so that the values can be parsed to float numbers. This version of the ATNF Pulsar Catalogue also contains a couple of negative parallax values, so this release does not allow these to be used for calculating distances.

## [1.0.3] 2019-07-29

Changes for this release:

- Updates to the [`BarycentricTrueEcliptic`](http://docs.astropy.org/en/stable/api/astropy.coordinates.BarycentricTrueEcliptic.html) in astropy v3.2 caused test to fail (due to ["correcting"](http://docs.astropy.org/en/v3.2.1/changelog.html#id12) bugs and including nutation). This version updates to use [`BarycentricMeanEcliptic`](http://docs.astropy.org/en/stable/api/astropy.coordinates.BarycentricMeanEcliptic.html), for which the tests no-longer fail. For earlier versions of astropy the former function is still used.
- The `version` attribute of the catalogue table is now set if loading a query from a previously read table (see [#53](https://github.com/mattpitkin/psrqpy/pull/53)).

## [1.0.0] 2018-11-15

This release involves major changes to the API.

- PSRQpy now downloads the full ATNF Pulsar Catalogue database and stores it internally.
- PSRQpy will no longer generate queries of the ATNF Pulsar Catalogue via the web interface.
- Derived parameters are calculated internally rather than being values returned be the catalogue web interface.
- Full string paper references are no longer included in the output table, but can be returned using the `parse_ref`
method of the `QueryATNF` class.
