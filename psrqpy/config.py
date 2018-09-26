"""
This submodules sets up common constants for use, such as the allowed pulsar parameters and various
URLs used for queries.
"""

import itertools

#: The default ATNF catalogue version.
ATNF_VERSION = '1.59'

#: The ATNF pulsar catalogue base URL.
ATNF_BASE_URL = r'http://www.atnf.csiro.au/people/pulsar/psrcat/'

#: The name of the tarball containing the entire catalogue database.
ATNF_TARBALL = ATNF_BASE_URL + r'downloads/psrcat_pkg.tar.gz'

#: The Jodrell Bank glitch catalogue table URL.
GLITCH_URL = r'http://www.jb.man.ac.uk/pulsar/glitches/gTable.html'

# pulsar parameters (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html) that can be
# queried. For each parameter there is a dictionary giving:
#  - 'ref': True if the parameter can have an associated reference in the ATNF catalogue
#  - 'err': True if the parameter can have an associated error value
#  - 'unit': a string giving the units for the parameter (to be used if generating an astropy table)
PSR_GENERAL = {
    'NAME':     {'ref': True,  'err': False, 'units': None},  # Pulsar name. Default: B name
    'JNAME':    {'ref': True,  'err': False, 'units': None},  # Pulsar name (J2000)
    'BNAME':    {'ref': True,  'err': False, 'units': None},  # Pulsar Besselian name
    'PSRJ':     {'ref': True,  'err': False, 'units': None},  # Pulsar name (J2000)
    'PSRB':     {'ref': True,  'err': False, 'units': None},  # Pulsar Besselian name
    'RAJ':      {'ref': True,  'err': True,  'units': None},  # Right ascension (J2000) (hh:mm:ss.s)
    'DECJ':     {'ref': True,  'err': True,  'units': None},  # Declination (J2000) (+dd:mm:ss)
    'PMRA':     {'ref': True,  'err': True,  'units': 'mas/yr'},   # Proper motion in right ascension
    'PMDEC':    {'ref': True,  'err': True,  'units': 'mas/yr'},   # Proper motion in declination
    'PX':       {'ref': True,  'err': True,  'units': 'mas'},   # Annual parallax 
    'POSEPOCH': {'ref': True,  'err': False, 'units': 'd'},   # Epoch of position, defaults to PEpoch (MJD)
    'ELONG':    {'ref': True,  'err': True,  'units': 'deg'},   # Ecliptic longitude
    'ELAT':     {'ref': True,  'err': True,  'units': 'deg'},   # Ecliptic latitude
    'PMELONG':  {'ref': True,  'err': True,  'units': 'mas/yr'},   # Proper motion in ecliptic long.
    'PMELAT':   {'ref': True,  'err': True,  'units': 'mas/yr'},   # Proper motion in ecliptic lat.
    'GL':       {'ref': False, 'err': False, 'units': 'deg'},   # Galactic longitude
    'GB':       {'ref': False, 'err': False, 'units': 'deg'},   # Galactic latitude
    'RAJD':     {'ref': False, 'err': False, 'units': 'deg'},   # Right ascension (J2000)
    'DECJD':    {'ref': False, 'err': False, 'units': 'deg'},   # Declination (J2000)
    'TYPE':     {'ref': False, 'err': False, 'units': None},  # Type codes for the pulsar
    'PML':      {'ref': False, 'err': False,  'units': 'mas/yr'},  # Proper motion in Galactic long.
    'PMB':      {'ref': False, 'err': False,  'units': 'mas/yr'},  # Proper motion in Galactic lat.
    'DIST':     {'ref': False, 'err': False, 'units': 'kpc'},   # Best estimate of the pulsar distance using the YMW16 DM-based distance as default
    'DIST_DM':  {'ref': True,  'err': False, 'units': 'kpc'},   # Distance based on the YMW16 electron density model
    'DIST_DM1': {'ref': True,  'err': False, 'units': 'kpc'},   # Distance based on NE2001 model
    'DIST1':    {'ref': False, 'err': False, 'units': 'kpc'},   # Best estimate of the pulsar distance using the NE2001 DM-based distance as default
    'DIST_AMN': {'ref': True,  'err': True,  'units': 'kpc'},   # Lower limit on independent distance estimate (NOTE: although there is an 'error' column it is always zero)
    'DIST_AMX': {'ref': True,  'err': False, 'units': 'kpc'},   # Upper limit on independent distance estimate
    'DIST_A':   {'ref': True,  'err': True,  'units': 'kpc'},   # Independent distance estimate
    'DMSINB':   {'ref': False, 'err': False, 'units': 'cm^-3 pc'},   # DM x sin(b)
    'ZZ':       {'ref': False, 'err': False, 'units': 'kpc'},   # Distance from the Galactic plane, based on Dist
    'XX':       {'ref': False, 'err': False, 'units': 'kpc'},   # X-Distance in X-Y-Z Galactic coordinate system
    'YY':       {'ref': False, 'err': False, 'units': 'kpc'},   # Y-Distance in X-Y-Z Galactic coordinate system
    'ASSOC':    {'ref': False, 'err': False, 'units': None},   # Name(s) of other objects associated with the pulsar
    'SURVEY':   {'ref': False, 'err': False, 'units': None},   # Surveys that detected the pulsar
    'OSURVEY':  {'ref': False, 'err': False, 'units': None},   # Surveys that detected the pulsar encoded as bits in integer
    'DATE':     {'ref': False, 'err': False, 'units': 'yr'},   # Date of discovery publication
    'NGLT':     {'ref': False, 'err': False, 'units': None},   # Number of glitches observed for the pulsar
    'GLEP':     {'ref': False, 'err': False, 'units': 'd'},   # Epoch of glitch
    'GLPH':     {'ref': False, 'err': False, 'units': None},   # Phase increment at glitch
    'GLF0':     {'ref': False, 'err': False, 'units': 'Hz'},   # Permanent pulse frequency increment at glitch
    'GLF1':     {'ref': False, 'err': False, 'units': 'Hz/s'},   # Permanent frequency derivative increment at glitch
    'GLF0D':    {'ref': False, 'err': False, 'units': 'Hz'},   # Decaying frequency increment at glitch
    'GLTD':     {'ref': False, 'err': False, 'units': 'd'},   # Time constant for decaying frequency increment
    'CLK':      {'ref': True,  'err': False, 'units': None},   # Reference clock used for timing solution
    'EPHEM':    {'ref': True,  'err': False, 'units': None}   # Solar-system ephemeris used for timing solution
    }

# just return the parameter names
PSR_GENERAL_PARS = list(PSR_GENERAL.keys())

# timing solution and profile parameters
PSR_TIMING = {
    'P0':      {'ref': True,  'err': True,  'units': 's'},      # Barycentric period
    'P1':      {'ref': True,  'err': True,  'units': None},     # Time derivative of barcycentric period
    'F0':      {'ref': True,  'err': True,  'units': 'Hz'},     # Barycentric rotation frequency
    'F1':      {'ref': True,  'err': True,  'units': 's^-2'},   # Time derivative of barycentric rotation frequency
    'F2':      {'ref': True,  'err': True,  'units': 's^-3'},   # Second time derivative of barycentric rotation frequency
    'F3':      {'ref': True,  'err': True,  'units': 's^-4'},   # Third time derivative of barycentric rotation frequency
    'F4':      {'ref': True,  'err': True,  'units': 's^-5'},   # Fourth time derivative of barycentric rotation frequency
    'F5':      {'ref': True,  'err': True,  'units': 's^-6'},   # Fifth time derivative of barycentric rotation frequency
    'PEPOCH':  {'ref': True,  'err': False, 'units': 'd'},      # Epoch of period or frequency (MJD)
    'DM':      {'ref': True,  'err': True,  'units': 'cm^-3 pc'},        # Dispersion measure
    'DM1':     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-1'},     # First time derivative of dispersion measure
    'DMEPOCH': {'ref': True,  'err': False, 'units': 'd'},  # Reference epoch for DM, defaults to PEpoch (MJD)
    'DM2':     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-2'},  # Second time derivative of dispersion measure
    'DM3':     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-3'},  # Third time derivative of dispersion measure
    'RM':      {'ref': True,  'err': True,  'units': 'rad m^-2'},  # Rotation measure
    'W50':     {'ref': True,  'err': True,  'units': 'ms'},  # Width of pulse at 50% of peak
    'W10':     {'ref': True,  'err': True,  'units': 'ms'},  # Width of pulse at 10% of peak
    'UNITS':   {'ref': True,  'err': False, 'units': None},  # Timescale for period/frequency and epoch data: TCB or TDB
    'TAU_SC':  {'ref': True,  'err': True,  'units': 's'},  # Temporal broadening of pulses at 1 GHz due to interestellar scattering
    'SI414':   {'ref': False, 'err': False, 'units': None},  # Spectral index between 400 and 1400 MHz
    'S400':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 400 MHz
    'S1400':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 1400 MHz
    'S2000':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 2000 MHz
    'S40':     {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 40 MHz
    'S50':     {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 50 MHz
    'S60':     {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 60 MHz
    'S80':     {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 80 MHz
    'S100':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 100 MHz
    'S150':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 150 MHz
    'S200':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 200 MHz
    'S300':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 300 MHz
    'S600':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 600 MHz
    'S700':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 700 MHz
    'S800':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 800 MHz
    'S900':    {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 900 MHz
    'S1600':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 1600 MHz
    'S3000':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 3000 MHz
    'S4000':   {'ref': True,  'err': False, 'units': 'mJy'},  # Mean flux at 4000 MHz
    'S5000':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 5000 MHz
    'S6000':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 6000 MHz
    'S8000':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 8000 MHz
    'S100G':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 100 GHz
    'S150G':   {'ref': True,  'err': True,  'units': 'mJy'},  # Mean flux at 150 GHz
    'SPINDX':  {'ref': True,  'err': True,  'units': None}   # Radio spectral index
    }

PSR_TIMING_PARS = list(PSR_TIMING.keys())

# binary system parameters
PSR_BINARY = {
    'BINARY':   {'ref': True,  'err': False, 'units': None},   # Binary model
    'T0':       {'ref': True,  'err': True,  'units': 'd'},   # Epoch of periastron (MJD)
    'PB':       {'ref': True,  'err': True,  'units': 'd'},   # Binary period of pulsar
    'A1':       {'ref': True,  'err': True,  'units': 's'},   # Projected semi-major axis
    'OM':       {'ref': True,  'err': True,  'units': 'deg'},   # Longitude of periastron
    'ECC':      {'ref': True,  'err': True,  'units': None},   # Eccentricity
    'TASC':     {'ref': True,  'err': True,  'units': 'd'},   # Epoch of ascending node (MJD)
    'EPS1':     {'ref': True,  'err': True,  'units': None},   # ECC x sin(OM)
    'EPS2':     {'ref': True,  'err': True,  'units': None},   # ECC x cos(OM)
    'MINMASS':  {'ref': False, 'err': False, 'units': 'M_sun'},    # Minimum companion mass
    'MEDMASS':  {'ref': False, 'err': False, 'units': 'M_sun'},    # Median companion mass
    'BINCOMP':  {'ref': True,  'err': False, 'units': None},    # Companion type
    'FB0':      {'ref': True,  'err': True,  'units': 'Hz'},    # Orbital frequency
    'FB1':      {'ref': True,  'err': True,  'units': 'Hz s^-1'},   # 1st time derivative of orbital frequency
    'FB2':      {'ref': True,  'err': True,  'units': 'Hz s^-2'},   # 2nd time derivative of orbital frequency
    'OMDOT':    {'ref': True,  'err': True,  'units': 'deg/yr'},   # 1st time derivative of periastron longitude
    'OM2DOT':   {'ref': True,  'err': True,  'units': 'deg/yr^2'},   # 2nd time derivative of periastron longitude
    'A1DOT':    {'ref': True,  'err': True,  'units': 's s^-1'},   # 1st time derivative of projected semi-major axis
    'A12DOT':   {'ref': True,  'err': True,  'units': 's s^-2'},   # 2nd time derivative of projected semi-major axis
    'ECCDOT':   {'ref': True,  'err': True,  'units': 's^-1'},   # 1st time derivative of eccentricity
    'PBDOT':    {'ref': True,  'err': True,  'units': None},   # 1st time derivative of binary period
    'GAMMA':    {'ref': True,  'err': True,  'units': 's'},   # Post-Keplerian time-dilation term
    'T0_2':     {'ref': True,  'err': True,  'units': 'd'},   # Epoch of periastron (2nd orbit; MJD)
    'PB_2':     {'ref': True,  'err': True,  'units': 'd'},   # Binary period of pulsar (2nd orbit)
    'A1_2':     {'ref': True,  'err': True,  'units': 's'},   # Projected semi-major axis of orbit (2nd orbit)
    'OM_2':     {'ref': True,  'err': True,  'units': 'deg'},   # Longitude of periastron (2nd orbit)
    'ECC_2':    {'ref': True,  'err': True,  'units': None},   # Eccentricity (2nd orbit)
    'EPS1_2':   {'ref': True,  'err': True,  'units': None},   # ECC_2 x sin(OM_2) (2nd orbit)
    'EPS2_2':   {'ref': True,  'err': True,  'units': None},   # ECC_2 x cos(OM_2) (2nd orbit)
    'TASC_2':   {'ref': True,  'err': True,  'units': 'd'},   # Epoch of ascending node (2nd orbit; MJD)
    'T0_3':     {'ref': True,  'err': True,  'units': 'd'},   # Epoch of periastron (3rd orbit; MJD)
    'PB_3':     {'ref': True,  'err': True,  'units': 'd'},   # Binary period of pulsar (3rd orbit)
    'A1_3':     {'ref': True,  'err': True,  'units': 's'},   # Projected semi-major axis of orbit (3rd orbit)
    'OM_3':     {'ref': True,  'err': True,  'units': 'deg'},   # Longitude of periastron (3rd orbit)
    'ECC_3':    {'ref': True,  'err': True,  'units': None},   # Eccentricity (3rd orbit)
    'SINI':     {'ref': True,  'err': True,  'units': None},   # Sine of inclination angle
    'SINI_2':   {'ref': True,  'err': True,  'units': None},   # Sine of inclination angle (2nd orbit)
    'SINI_3':   {'ref': True,  'err': True,  'units': None},   # Sine of inclination angle (3rd orbit)
    'KOM':      {'ref': True,  'err': True,  'units': 'deg'},   # Long. on sky of asc. node (N toward E) from ann. orbital parallax
    'KIN':      {'ref': True,  'err': True,  'units': 'deg'},   # Orbit inclination from annual orbital parallax
    'M2':       {'ref': True,  'err': True,  'units': 'M_sun'},   # Companion mass
    'M2_2':     {'ref': True,  'err': True,  'units': 'M_sun'},   # Companion mass (2nd orbit)
    'M2_3':     {'ref': True,  'err': True,  'units': 'M_sun'},   # Companion mass (3rd orbit)
    'MASS_Q':   {'ref': True,  'err': True,  'units': None},   # Mass ratio for binary: M1/M2
    'OM_ASC':   {'ref': True,  'err': True,  'units': 'deg'},    # Longitude on sky of ascending node (from N toward E)
    'DTHETA':   {'ref': True,  'err': True,  'units': '1e-6'},    # Relativistic deformation of the orbit
    'XOMDOT':   {'ref': False, 'err': False, 'units': 'deg/yr'},   # Rate of periastron advance minus GR prediction
    'H3':       {'ref': True,  'err': True,  'units': 's'},   # Amplitude of 3rd Shapiro-delay harmonic
    'H4':       {'ref': True,  'err': True,  'units': 's'},   # Amplitude of 4th Shapiro-delay harmonic
    'STIG':     {'ref': True,  'err': True,  'units': None},   # Ratio of successive Shapiro-delay harmonics
    'MASSFN':   {'ref': False, 'err': False, 'units': 'M_sun'},    # Pulsar mass function
    'UPRMASS':  {'ref': False, 'err': False, 'units': 'M_sun'},    # 90% confidence upper companion mass limit
    'MINOMDOT': {'ref': False, 'err': False, 'units': 'deg/yr'}    # Minimum omega dot, assuming sin i = 1 and M_ns = 1.4M_sol
    }

PSR_BINARY_PARS = list(PSR_BINARY.keys())

# derived parameters
PSR_DERIVED = {
    'R_LUM':    {'ref': False, 'err': False, 'units': 'mJy kpc^2'},   # Radio luminosity at 400 MHz
    'R_LUM14':  {'ref': False, 'err': False, 'units': 'mJy kpc^2'},   # Radio luminosity at 1400 MHz
    'AGE':      {'ref': False, 'err': False, 'units': 'yr'},   # Spin down age
    'BSURF':    {'ref': False, 'err': False, 'units': 'G'},   # Surface magnetic flux density
    'EDOT':     {'ref': False, 'err': False, 'units': 'erg/s'},   # Spin down energy loss rate
    'EDOTD2':   {'ref': False, 'err': False, 'units': 'erg s^-1/kpc^2'},   # Energy flux at the Sun
    'PMTOT':    {'ref': False, 'err': False,  'units': 'mas/yr'},   # Total proper motion
    'VTRANS':   {'ref': False, 'err': False, 'units': 'km/s'},    # Transverse velocity - based on DIST
    'P1_I':     {'ref': False, 'err': False, 'units': None},   # Period derivative corrected for Shklovskii effect
    'AGE_I':    {'ref': False, 'err': False, 'units': 'yr'},   # Spin down age from P1_i
    'BSURF_I':  {'ref': False, 'err': False, 'units': 'G'},   # Surface magnetic dipole from P1_i
    'B_LC':     {'ref': False, 'err': False, 'units': 'G'}   # Magnetic field at light cylinder
    }

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
  :class:`~astropy.units.core.Unit`

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
""": Allowed pulsar
`types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_
for use in ``type()`` when setting logical conditions.
"""

PSR_BINARY_TYPE = [
    'MS',   # Main-sequence star
    'NS',   # Neutron star
    'CO',   # CO or ONeMg White Dwarf
    'He',   # Helium White Dwarf
    'UL'   # Ultra-light companion or planet (mass < 0.08 solar masses)
    ]
""": Binary companion
`types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type>`_ 
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
""": Other objects associated with the pulsar.
"""

#: The URL for the NASA ADS.
ADS_URL = 'https://ui.adsabs.harvard.edu/#abs/{}/'
