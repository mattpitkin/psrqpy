# Notable changes between versions

## [1.1.7] 2022-02-3

- Set position epoch to period epoch if the latter is not given. See [#96](https://github.com/mattpitkin/psrqpy/pull/96). 
- Allow discovery data tp be derived from the catalogue. See [#97](https://github.com/mattpitkin/psrqpy/pull/97). 

## [1.1.6] 2021-12-13

- Allow specific versions of the ATNF pulsar catalogue to be downloaded. See [#95](https://github.com/mattpitkin/psrqpy/pull/95).

## [1.1.5] 2021-11-18

- Extend `__getitem__` for `QueryATNF` class so that individual pulsar can be extracted using keys from a class instance. See [#92](https://github.com/mattpitkin/psrqpy/pull/92).

## [1.1.4] 2021-11-17

- Propagate errors on ecliptical coordinates to equatorial. See [#90](https://github.com/mattpitkin/psrqpy/pull/90).
- Update build system. See [#91](https://github.com/mattpitkin/psrqpy/pull/91).

## [1.1.3] 2021-10-11

- Move scraping of Galactic MSP catalogue into PSRQPy. See [#88](https://github.com/mattpitkin/psrqpy/pull/88).

## [1.1.2] 2021-06-22

- Fix parsing of the galactic MSP table after changes to the upstream format. See [#85](https://github.com/mattpitkin/psrqpy/pull/85).
- Remove explicit use of certain NumPy built-in types following their [deprecation](https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations). See [#83](https://github.com/mattpitkin/psrqpy/pull/83).

## [1.1.1] 2021-03-19

- Add function to download and parse the galactic MSP table. See [#80](https://github.com/mattpitkin/psrqpy/pull/80).
- Add functions to calculate derived gravitational wave parameters. See [#79](https://github.com/mattpitkin/psrqpy/pull/79).

## [1.1.0] 2021-03-18

The major change for this release is that it no longer supports Python 2.7, and only supports Python versions
greater than 3.5.

- Add a function to download and parse the [globular cluster pulsar](http://www.naic.edu/~pfreire/GCpsr.html) table. See [#77](https://github.com/mattpitkin/psrqpy/pull/77).

## [1.0.11] 2020-10-19

- Allow `psrqpy` query to work using ATNF catalogue v1.64, which contains a position typo. See [#74](https://github.com/mattpitkin/psrqpy/pull/74).
- More fixes to enable parsing of more references using NASA ADS and fixes changes that occurred with v1.64. See [#70](https://github.com/mattpitkin/psrqpy/pull/70).

## [1.0.10] 2020-10-08

- Allow the `Pulsar` objects to access parameter references. See [#71](https://github.com/mattpitkin/psrqpy/pull/71).

## [1.0.9] 2020-06-08

- A minor fix to get ADS references to work. See [#66](https://github.com/mattpitkin/psrqpy/pull/66).

## [1.0.8] 2020-06-05

Changes for this release:

- The way references are stored in the ATNF pulsar catalogue has changed, so the `get_references` function has been changed to work with this new format. This breaks that function for older cached catalogue, but a warning is provided. See [#64](https://github.com/mattpitkin/psrqpy/pull/64).
- The astropy Galactocentric coordinates have been set to use default values from the pre-v4.0 release. See [#65](https://github.com/mattpitkin/psrqpy/pull/65).

## [1.0.7] 2020-03-29

Changes for this release:

- Fix an issue with Tables in Astropy v4.0 (see [#61](https://github.com/mattpitkin/psrqpy/pull/61)).

## [1.0.6] 2020-03-19

Changes for this release:

- Adds position errors on the right ascension and declination when converted into degrees (see [#60](https://github.com/mattpitkin/psrqpy/pull/60)).

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
