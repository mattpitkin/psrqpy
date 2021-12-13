"""
This submodule sets up common constants for use, such as the allowed pulsar parameters and various
URLs used for queries.
"""

import itertools


#: The ATNF pulsar catalogue base URL.
ATNF_BASE_URL = r"http://www.atnf.csiro.au/people/pulsar/psrcat/"

#: The name of the tarball containing the entire catalogue database.
ATNF_TARBALL = ATNF_BASE_URL + r"downloads/psrcat_pkg.tar.gz"

# The name of the tarball containing the entire catalogue database (allowing version string to
# be added).
ATNF_VERSION_TARBALL = ATNF_BASE_URL + r"downloads/psrcat_pkg.v{}.tar.gz"

#: The Jodrell Bank glitch catalogue table URL.
GLITCH_URL = r"http://www.jb.man.ac.uk/pulsar/glitches/gTable.html"

#: Paolo Freire's globular cluster pulsar table URL
GC_URL = r"http://www.naic.edu/~pfreire/GCpsr.txt"

#: Dunc Lorimer's MSP table URL
MSP_URL = r"http://astro.phys.wvu.edu/GalacticMSPs/GalacticMSPs.txt"

# Pulsar parameters (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html) that can be
# queried. For each parameter there is a dictionary giving:
#  - 'ref': True if the parameter can have an associated reference in the ATNF catalogue
#  - 'err': True if the parameter can have an associated error value
#  - 'unit': a string giving the units for the parameter (to be used if generating an astropy table)
PSR_GENERAL = dict()
PSR_GENERAL['NAME'] =     {'ref': True,  'err': False, 'units': None}  # Pulsar name. Default: B name
PSR_GENERAL['JNAME'] =    {'ref': True,  'err': False, 'units': None}  # Pulsar name (J2000)
PSR_GENERAL['BNAME'] =    {'ref': True,  'err': False, 'units': None}  # Pulsar Besselian name
PSR_GENERAL['PSRJ'] =     {'ref': True,  'err': False, 'units': None}  # Pulsar name (J2000)
PSR_GENERAL['PSRB'] =     {'ref': True,  'err': False, 'units': None}  # Pulsar Besselian name
PSR_GENERAL['RAJ'] =      {'ref': True,  'err': True,  'units': None}  # Right ascension (J2000) (hh:mm:ss.s)
PSR_GENERAL['DECJ'] =     {'ref': True,  'err': True,  'units': None}  # Declination (J2000) (+dd:mm:ss)
PSR_GENERAL['PMRA'] =     {'ref': True,  'err': True,  'units': 'mas/yr'}   # Proper motion (right ascension)
PSR_GENERAL['PMDEC'] =    {'ref': True,  'err': True,  'units': 'mas/yr'}   # Proper motion (declination)
PSR_GENERAL['PX'] =       {'ref': True,  'err': True,  'units': 'mas'}   # Annual parallax
PSR_GENERAL['POSEPOCH'] = {'ref': True,  'err': False, 'units': 'd'}   # Epoch of position. Default: PEpoch
PSR_GENERAL['ELONG'] =    {'ref': True,  'err': True,  'units': 'deg'}   # Ecliptic longitude
PSR_GENERAL['ELAT'] =     {'ref': True,  'err': True,  'units': 'deg'}   # Ecliptic latitude
PSR_GENERAL['PMELONG'] =  {'ref': True,  'err': True,  'units': 'mas/yr'}   # Proper motion in ecliptic long.
PSR_GENERAL['PMELAT'] =   {'ref': True,  'err': True,  'units': 'mas/yr'}   # Proper motion in ecliptic lat.
PSR_GENERAL['GL'] =       {'ref': False, 'err': False, 'units': 'deg'}   # Galactic longitude
PSR_GENERAL['GB'] =       {'ref': False, 'err': False, 'units': 'deg'}   # Galactic latitude
PSR_GENERAL['RAJD'] =     {'ref': False, 'err': True,  'units': 'deg'}   # Right ascension (J2000)
PSR_GENERAL['DECJD'] =    {'ref': False, 'err': True,  'units': 'deg'}   # Declination (J2000)
PSR_GENERAL['TYPE'] =     {'ref': False, 'err': False, 'units': None}  # Type codes for the pulsar
PSR_GENERAL['PML'] =      {'ref': False, 'err': False,  'units': 'mas/yr'}  # Proper motion in Galactic long.
PSR_GENERAL['PMB'] =      {'ref': False, 'err': False,  'units': 'mas/yr'}  # Proper motion in Galactic lat.
# DIST: Best estimate of the pulsar distance using the YMW16 DM-based distance as default
PSR_GENERAL['DIST'] =     {'ref': False, 'err': False, 'units': 'kpc'}
# DIST_DM: Distance based on the YMW16 electron density model
PSR_GENERAL['DIST_DM'] =  {'ref': True,  'err': False, 'units': 'kpc'}
# DIST_DM1: Distance based on NE2001 electron density model
PSR_GENERAL['DIST_DM1'] = {'ref': True,  'err': False, 'units': 'kpc'}
# DIST1: Best estimate of the pulsar distance using the NE2001 DM-based distance as default
PSR_GENERAL['DIST1'] =    {'ref': False, 'err': False, 'units': 'kpc'}
# DIST_AMN: Lower limit on independent distance estimate (NOTE: 'error' column is always zero)
PSR_GENERAL['DIST_AMN'] = {'ref': True,  'err': True,  'units': 'kpc'}
# DIST_AMX: Upper limit on independent distance estimate
PSR_GENERAL['DIST_AMX'] = {'ref': True,  'err': False, 'units': 'kpc'}
# DIST_A: Independent distance estimate
PSR_GENERAL['DIST_A'] =   {'ref': True,  'err': True,  'units': 'kpc'}
PSR_GENERAL['DMSINB'] =   {'ref': False, 'err': False, 'units': 'cm^-3 pc'}   # DM x sin(b)
# ZZ: Distance from the Galactic plane, based on Dist
PSR_GENERAL['ZZ'] =       {'ref': False, 'err': False, 'units': 'kpc'}
# XX: X-Distance in X-Y-Z Galactic coordinate system
PSR_GENERAL['XX'] =       {'ref': False, 'err': False, 'units': 'kpc'}
# YY: Y-Distance in X-Y-Z Galactic coordinate system
PSR_GENERAL['YY'] =       {'ref': False, 'err': False, 'units': 'kpc'}
# ASSOC: Name(s) of other objects associated with the pulsar
PSR_GENERAL['ASSOC'] =    {'ref': False, 'err': False, 'units': None}
PSR_GENERAL['SURVEY'] =   {'ref': False, 'err': False, 'units': None}   # Surveys that detected the pulsar
# OSURVEY: Surveys that detected the pulsar, encoded as bits in integer
PSR_GENERAL['OSURVEY'] =  {'ref': False, 'err': False, 'units': None}
PSR_GENERAL['DATE'] =     {'ref': False, 'err': False, 'units': 'yr'}   # Date of discovery publication
# NGLT: Number of glitches observed for the pulsar
PSR_GENERAL['NGLT'] =     {'ref': False, 'err': False, 'units': None}
# GLEP: Epoch of glitch (MJD)
PSR_GENERAL['GLEP'] =     {'ref': False, 'err': False, 'units': 'd'}
# GLPH: Phase increment at glitch
PSR_GENERAL['GLPH'] =     {'ref': False, 'err': False, 'units': None}
# GLF0: Permanent pulse frequency increment at glitch
PSR_GENERAL['GLF0'] =     {'ref': False, 'err': False, 'units': 'Hz'}
# GLF1: Permanent frequency derivative increment at glitch
PSR_GENERAL['GLF1'] =     {'ref': False, 'err': False, 'units': 'Hz/s'}
# GLF0D: Decaying frequency increment at glitch
PSR_GENERAL['GLF0D'] =    {'ref': False, 'err': False, 'units': 'Hz'}
# GLTD: Time constant for decaying frequency increment
PSR_GENERAL['GLTD'] =     {'ref': False, 'err': False, 'units': 'd'}
# CLK: Reference clock used for timing solution
PSR_GENERAL['CLK'] =      {'ref': True,  'err': False, 'units': None}
# EPHEM: Solar-system ephemeris used for timing solution
PSR_GENERAL['EPHEM'] =    {'ref': True,  'err': False, 'units': None}

# just return the parameter names
PSR_GENERAL_PARS = list(PSR_GENERAL.keys())

# timing solution and profile parameters
PSR_TIMING = dict()
PSR_TIMING['P0'] =      {'ref': True,  'err': True,  'units': 's'}  # Barycentric period
# P1: Time derivative of barcycentric period
PSR_TIMING['P1'] =      {'ref': True,  'err': True,  'units': None}
# F0: Barycentric rotation frequency
PSR_TIMING['F0'] =      {'ref': True,  'err': True,  'units': 'Hz'}
# F1: Time derivative of barycentric rotation frequency
PSR_TIMING['F1'] =      {'ref': True,  'err': True,  'units': 's^-2'}
# F2: Second time derivative of barycentric rotation frequency
PSR_TIMING['F2'] =      {'ref': True,  'err': True,  'units': 's^-3'}
# F3: Third time derivative of barycentric rotation frequency
PSR_TIMING['F3'] =      {'ref': True,  'err': True,  'units': 's^-4'}
# F4: Fourth time derivative of barycentric rotation frequency
PSR_TIMING['F4'] =      {'ref': True,  'err': True,  'units': 's^-5'}
# F5: Fifth time derivative of barycentric rotation frequency
PSR_TIMING['F5'] =      {'ref': True,  'err': True,  'units': 's^-6'}
# PEPOCH: Epoch of period or frequency (MJD)
PSR_TIMING['PEPOCH'] =  {'ref': True,  'err': False, 'units': 'd'}
# DMEPOCH: Reference epoch for DM, defaults to PEpoch (MJD)
PSR_TIMING['DMEPOCH'] = {'ref': True,  'err': False, 'units': 'd'}
PSR_TIMING['DM'] =      {'ref': True,  'err': True,  'units': 'cm^-3 pc'}  # Dispersion measure
# DM1: First time derivative of dispersion measure
PSR_TIMING['DM1'] =     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-1'}
# DM2: Second time derivative of dispersion measure
PSR_TIMING['DM2'] =     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-2'}
# DM3: Third time derivative of dispersion measure
PSR_TIMING['DM3'] =     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-3'}
PSR_TIMING['RM'] =      {'ref': True,  'err': True,  'units': 'rad m^-2'}  # Rotation measure
PSR_TIMING['W50'] =     {'ref': True,  'err': True,  'units': 'ms'}  # Width of pulse at 50% of peak
PSR_TIMING['W10'] =     {'ref': True,  'err': True,  'units': 'ms'}  # Width of pulse at 10% of peak
# UNITS: Timescale for period/frequency and epoch data: TCB or TDB
PSR_TIMING['UNITS'] =   {'ref': True,  'err': False, 'units': None}
# TAU_SC: Temporal broadening of pulses at 1 GHz due to interestellar scattering
PSR_TIMING['TAU_SC'] =  {'ref': True,  'err': True,  'units': 's'}
# SI414: Spectral index derived from S400 and S1400
PSR_TIMING['SI414'] =   {'ref': False, 'err': False, 'units': None}
PSR_TIMING['S400'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 400 MHz
PSR_TIMING['S1400'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 1400 MHz
PSR_TIMING['S2000'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 2000 MHz
PSR_TIMING['S40'] =     {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 40 MHz
PSR_TIMING['S50'] =     {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 50 MHz
PSR_TIMING['S60'] =     {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 60 MHz
PSR_TIMING['S80'] =     {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 80 MHz
PSR_TIMING['S100'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 100 MHz
PSR_TIMING['S150'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 150 MHz
PSR_TIMING['S200'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 200 MHz
PSR_TIMING['S300'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 300 MHz
PSR_TIMING['S600'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 600 MHz
PSR_TIMING['S700'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 700 MHz
PSR_TIMING['S800'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 800 MHz
PSR_TIMING['S900'] =    {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 900 MHz
PSR_TIMING['S1600'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 1600 MHz
PSR_TIMING['S3000'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 3000 MHz
PSR_TIMING['S4000'] =   {'ref': True,  'err': False, 'units': 'mJy'}  # Mean flux at 4000 MHz
PSR_TIMING['S5000'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 5000 MHz
PSR_TIMING['S6000'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 6000 MHz
PSR_TIMING['S8000'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 8000 MHz
PSR_TIMING['S100G'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 100 GHz
PSR_TIMING['S150G'] =   {'ref': True,  'err': True,  'units': 'mJy'}  # Mean flux at 150 GHz
PSR_TIMING['SPINDX'] =  {'ref': True,  'err': True,  'units': None}   # Radio spectral index

PSR_TIMING_PARS = list(PSR_TIMING.keys())

# binary system parameters
PSR_BINARY = dict()
PSR_BINARY['BINARY'] =   {'ref': True,  'err': False, 'units': None}   # Binary model
PSR_BINARY['T0'] =       {'ref': True,  'err': True,  'units': 'd'}   # Epoch of periastron (MJD)
PSR_BINARY['PB'] =       {'ref': True,  'err': True,  'units': 'd'}   # Binary period of pulsar
PSR_BINARY['A1'] =       {'ref': True,  'err': True,  'units': 's'}   # Projected semi-major axis
PSR_BINARY['OM'] =       {'ref': True,  'err': True,  'units': 'deg'}   # Longitude of periastron
PSR_BINARY['ECC'] =      {'ref': True,  'err': True,  'units': None}   # Eccentricity
PSR_BINARY['TASC'] =     {'ref': True,  'err': True,  'units': 'd'}   # Epoch of ascending node (MJD)
PSR_BINARY['EPS1'] =     {'ref': True,  'err': True,  'units': None}   # ECC x sin(OM)
PSR_BINARY['EPS2'] =     {'ref': True,  'err': True,  'units': None}   # ECC x cos(OM)
PSR_BINARY['BINCOMP'] =  {'ref': True,  'err': False, 'units': None}    # Companion type
PSR_BINARY['FB0'] =      {'ref': True,  'err': True,  'units': 'Hz'}    # Orbital frequency
# FB1: 1st time derivative of orbital frequency
PSR_BINARY['FB1'] =      {'ref': True,  'err': True,  'units': 'Hz s^-1'}
# FB2: 2nd time derivative of orbital frequency
PSR_BINARY['FB2'] =      {'ref': True,  'err': True,  'units': 'Hz s^-2'}
# OMDOT: 1st time derivative of periastron longitude
PSR_BINARY['OMDOT'] =    {'ref': True,  'err': True,  'units': 'deg/yr'}
# OM2DOT: 2nd time derivative of periastron longitude
PSR_BINARY['OM2DOT'] =   {'ref': True,  'err': True,  'units': 'deg/yr^2'}
# A1DOT: 1st time derivative of projected semi-major axis
PSR_BINARY['A1DOT'] =    {'ref': True,  'err': True,  'units': 's s^-1'}
# A12DOT: 2nd time derivative of projected semi-major axis
PSR_BINARY['A12DOT'] =   {'ref': True,  'err': True,  'units': 's s^-2'}
# ECCDOT: 1st time derivative of eccentricity
PSR_BINARY['ECCDOT'] =   {'ref': True,  'err': True,  'units': 's^-1'}
# PBDOT: 1st time derivative of binary period
PSR_BINARY['PBDOT'] =    {'ref': True,  'err': True,  'units': None}
PSR_BINARY['GAMMA'] =    {'ref': True,  'err': True,  'units': 's'}   # Post-Keplerian time-dilation term
PSR_BINARY['T0_2'] =     {'ref': True,  'err': True,  'units': 'd'}   # Epoch of periastron (2nd orbit; MJD)
PSR_BINARY['PB_2'] =     {'ref': True,  'err': True,  'units': 'd'}   # Binary period of pulsar (2nd orbit)
# A1_2: Projected semi-major axis of orbit (2nd orbit)
PSR_BINARY['A1_2'] =     {'ref': True,  'err': True,  'units': 's'}
# OM_2: Longitude of periastron (2nd orbit)
PSR_BINARY['OM_2'] =     {'ref': True,  'err': True,  'units': 'deg'}
PSR_BINARY['ECC_2'] =    {'ref': True,  'err': True,  'units': None}   # Eccentricity (2nd orbit)
PSR_BINARY['EPS1_2'] =   {'ref': True,  'err': True,  'units': None}   # ECC_2 x sin(OM_2) (2nd orbit)
PSR_BINARY['EPS2_2'] =   {'ref': True,  'err': True,  'units': None}   # ECC_2 x cos(OM_2) (2nd orbit)
# TASC_2: Epoch of ascending node (2nd orbit; MJD)
PSR_BINARY['TASC_2'] =   {'ref': True,  'err': True,  'units': 'd'}
PSR_BINARY['T0_3'] =     {'ref': True,  'err': True,  'units': 'd'}   # Epoch of periastron (3rd orbit; MJD)
PSR_BINARY['PB_3'] =     {'ref': True,  'err': True,  'units': 'd'}   # Binary period of pulsar (3rd orbit)
# A1_3: Projected semi-major axis of orbit (3rd orbit)
PSR_BINARY['A1_3'] =     {'ref': True,  'err': True,  'units': 's'}
# OM_3: Longitude of periastron (3rd orbit)
PSR_BINARY['OM_3'] =     {'ref': True,  'err': True,  'units': 'deg'}
PSR_BINARY['ECC_3'] =    {'ref': True,  'err': True,  'units': None}   # Eccentricity (3rd orbit)
PSR_BINARY['SINI'] =     {'ref': True,  'err': True,  'units': None}   # Sine of inclination angle
# SINI_2: Sine of inclination angle (2nd orbit)
PSR_BINARY['SINI_2'] =   {'ref': True,  'err': True,  'units': None}
# SINI_3: Sine of inclination angle (3rd orbit)
PSR_BINARY['SINI_3'] =   {'ref': True,  'err': True,  'units': None}
# KOM: Long. on sky of asc. node (N toward E) from ann. orbital parallax
PSR_BINARY['KOM'] =      {'ref': True,  'err': True,  'units': 'deg'}
# KIN: Orbit inclination from annual orbital parallax
PSR_BINARY['KIN'] =      {'ref': True,  'err': True,  'units': 'deg'}
PSR_BINARY['M2'] =       {'ref': True,  'err': True,  'units': 'M_sun'}   # Companion mass
PSR_BINARY['M2_2'] =     {'ref': True,  'err': True,  'units': 'M_sun'}   # Companion mass (2nd orbit)
PSR_BINARY['M2_3'] =     {'ref': True,  'err': True,  'units': 'M_sun'}   # Companion mass (3rd orbit)
PSR_BINARY['MASS_Q'] =   {'ref': True,  'err': True,  'units': None}   # Mass ratio for binary: M1/M2
# OM_ASC: Longitude on sky of ascending node (from N toward E)
PSR_BINARY['OM_ASC'] =   {'ref': True,  'err': True,  'units': 'deg'}
# DTHETA: Relativistic deformation of the orbit
PSR_BINARY['DTHETA'] =   {'ref': True,  'err': True,  'units': '1e-6'}
# XOMDOT: Rate of periastron advance minus GR prediction
PSR_BINARY['XOMDOT'] =   {'ref': False, 'err': False, 'units': 'deg/yr'}
# H3: Amplitude of 3rd Shapiro-delay harmonic
PSR_BINARY['H3'] =       {'ref': True,  'err': True,  'units': 's'}
# H4: Amplitude of 4th Shapiro-delay harmonic
PSR_BINARY['H4'] =       {'ref': True,  'err': True,  'units': 's'}
# STIG: Ratio of successive Shapiro-delay harmonics
PSR_BINARY['STIG'] =     {'ref': True,  'err': True,  'units': None}
PSR_BINARY['MASSFN'] =   {'ref': False, 'err': True,  'units': 'M_sun'}    # Pulsar mass function
PSR_BINARY['MINMASS'] =  {'ref': False, 'err': False, 'units': 'M_sun'}    # Minimum companion mass
PSR_BINARY['MEDMASS'] =  {'ref': False, 'err': False, 'units': 'M_sun'}    # Median companion mass
# UPRMASS: 90% confidence upper companion mass limit
PSR_BINARY['UPRMASS'] =  {'ref': False, 'err': False, 'units': 'M_sun'}
# MINOMDOT: Minimum omega dot, assuming sin i = 1 and M_ns = 1.4M_sol
PSR_BINARY['MINOMDOT'] = {'ref': False, 'err': False, 'units': 'deg/yr'}

PSR_BINARY_PARS = list(PSR_BINARY.keys())

# derived parameters
PSR_DERIVED = dict()
PSR_DERIVED['R_LUM'] =   {'ref': False, 'err': False, 'units': 'mJy kpc^2'}   # Radio luminosity at 400 MHz
PSR_DERIVED['R_LUM14'] = {'ref': False, 'err': False, 'units': 'mJy kpc^2'}   # Radio luminosity at 1400 MHz
PSR_DERIVED['AGE'] =     {'ref': False, 'err': False, 'units': 'yr'}   # Spin down age
PSR_DERIVED['BSURF'] =   {'ref': False, 'err': False, 'units': 'G'}   # Surface magnetic flux density
PSR_DERIVED['EDOT'] =    {'ref': False, 'err': False, 'units': 'erg/s'}   # Spin down energy loss rate
PSR_DERIVED['EDOTD2'] =  {'ref': False, 'err': False, 'units': 'erg s^-1/kpc^2'}   # Energy flux at the Sun
PSR_DERIVED['PMTOT'] =   {'ref': False, 'err': False,  'units': 'mas/yr'}   # Total proper motion
# VTRANS: Transverse velocity - based on DIST
PSR_DERIVED['VTRANS'] =  {'ref': False, 'err': False, 'units': 'km/s'}
# P1_I: Period derivative corrected for Shklovskii effect
PSR_DERIVED['P1_I'] =    {'ref': False, 'err': False, 'units': None}
PSR_DERIVED['AGE_I'] =   {'ref': False, 'err': False, 'units': 'yr'}   # Spin down age from P1_i
PSR_DERIVED['BSURF_I'] = {'ref': False, 'err': False, 'units': 'G'}   # Surface magnetic dipole from P1_i
PSR_DERIVED['B_LC'] =    {'ref': False, 'err': False, 'units': 'G'}   # Magnetic field at light cylinder
PSR_DERIVED['H0_SD'] =   {'ref': False, 'err': False, 'units': None}  # GW spin-down limit

PSR_DERIVED_PARS = list(PSR_DERIVED.keys())

# a list of all allowed parameters for querying
PSR_ALL = dict(itertools.chain(PSR_GENERAL.items(), PSR_TIMING.items(),
                               PSR_BINARY.items(), PSR_DERIVED.items()))
""": A dictionary of allowed pulsars parameters (e.g., name, position,
distance...)

Each parameter name key gives a dictionary containing the keys:

* ``ref`` (bool) - True if the parameter has an associated reference with it
* ``err`` (bool) - True if the parameter has an associated error value
* ``units`` (str) - a string with the parameters units that can be parsed by
  :class:`~astropy.units.Unit`

The allowed parameters and their units are given
`here <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#par_list>`_.
"""

PSR_ALL_PARS = PSR_GENERAL_PARS + PSR_TIMING_PARS + PSR_BINARY_PARS + PSR_DERIVED_PARS

PSR_TYPE = [
    'AXP',   # Anomalous X-ray Pulsar or Soft Gamma-ray Repeater with detected pulsations
    'BINARY',   # Pulsar has one or more stellar companion(s)
    'HE',   # Spin-powered pulsar with pulsed emission from radio to infrared or higher frequencies
    'NRAD',   # Spin-powered pulsar with pulsed emission only at infrared or higher frequencies
    'RADIO',   # Pulsars with pulsed emission in the radio band
    'RRAT',   # Pulsars with intermittently pulsed radio emission
    'XINS'  # Isolated neutron stars with pulsed X-ray emission but no detectable radio emission
    ]
""": Allowed
`pulsar types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_
for use in ``type()`` when setting logical conditions.
"""

PSR_BINARY_TYPE = [
    'MS',   # Main-sequence star
    'NS',   # Neutron star
    'CO',   # CO or ONeMg White Dwarf
    'He',   # Helium White Dwarf
    'UL'   # Ultra-light companion or planet (mass < 0.08 solar masses)
    ]
""": `Binary companion
types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type>`_
for use in ``bincomp()`` when setting logical conditions.
"""

PSR_ASSOC_TYPE = [
    'AXP',  # Anomolous X-ray pulsar
    'EXGAL',  # Extra-galactic source
    'GC',   # Globular cluster
    'GRS',   # Gamma-ray source
    'OPT',   # Optical source
    'PWN',   # Pulsar wind nebula
    'SNR',   # Supernova remnant
    'XRS'   # X-ray source
    ]
""": Other sources/objects associated with the pulsar.
"""

#: The URL for the NASA ADS.
ADS_URL = 'https://ui.adsabs.harvard.edu/#abs/{}/'
