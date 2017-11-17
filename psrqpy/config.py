"""
Configuration information
"""

import itertools

ATNF_VERSION = '1.57' # default ATNF version
ATNF_BASE_URL = r'http://www.atnf.csiro.au/people/pulsar/psrcat/'
ATNF_URL = ATNF_BASE_URL + r'proc_form.php?version={version}'

PARAMS_QUERY = '{params}'
USERDEFINED_QUERY = '&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val='
CONDITION_QUERY = '&condition={condition}'
PSRNAMES_QUERY = '&pulsar_names={psrnames}'
SORT_QUERY = '&sort_attr={sortattr}&sort_order={sortorder}'
EPHEMERIS_QUERY = '&submit_ephemeris={getephemeris}'
QUERY_FLUFF = '&ephemeris=long&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=*&nohead=nohead&state=query&table_bottom.x=30&table_bottom.y=22'

QUERY_URL = ATNF_URL + PARAMS_QUERY + USERDEFINED_QUERY + SORT_QUERY + CONDITION_QUERY + PSRNAMES_QUERY + EPHEMERIS_QUERY + QUERY_FLUFF

# pulsar parameters (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html) that can be queried
# For each parameter there is a dictionary giving:
#  - 'ref': True if the parameter can have an associated reference in the ATNF catalogue
#  - 'err': True if the parameter can have an associated error value
#  - 'unit': a string giving the units for the parameter (to be used if generating an astropy table)
#  - 'format': a string giving the parameter format (to be used if generating an astropy table)

# general parameters (e.g., name, position, distance...)
PSR_GENERAL = {'Name':     {'ref': True, 'err': False, 'units': None, 'format': 'S32'},        # Pulsar name.  The B name if exists, otherwise the J name.
               'JName':    {'ref': True, 'err': False, 'units': None, 'format': 'S32'},        # Pulsar name based on J2000 coordinates
               'RAJ':      {'ref': True, 'err': True, 'units': None, 'format': 'S32'},         # Right ascension (J2000) (hh:mm:ss.s)
               'DecJ':     {'ref': True, 'err': True, 'units': None, 'format': 'S32'},         # Declination (J2000) (+dd:mm:ss)
               'PMRA':     {'ref': True, 'err': True, 'units': 'mas/yr', 'format': 'f8'},      # Proper motion in the right ascension direction (mas/yr)
               'PMDec':    {'ref': True, 'err': True, 'units': 'mas/yr', 'format': 'f8'},      # Proper motion in declination (mas/yr)
               'PX':       {'ref': True, 'err': True, 'units': 'mas', 'format': 'f8'},         # Annual parallax (mas)
               'PosEpoch': {'ref': True, 'err': True, 'units': 'd', 'format': 'f8'},           # Epoch of position, defaults to PEpoch (MJD)
               'ELong':    {'ref': True, 'err': True, 'units': 'deg', 'format': 'f8'},         # Ecliptic longitude (degrees)
               'ELat':     {'ref': True, 'err': True, 'units': 'deg', 'format': 'f8'},         # Ecliptic latitude (degrees)
               'PMElong':  {'ref': True, 'err': True, 'units': 'mas/yr', 'format': 'f8'},      # Proper motion in the ecliptic longitude direction (mas/yr)
               'PMElat':   {'ref': True, 'err': True, 'units': 'mas/yr', 'format': 'f8'},      # Proper motion in ecliptic latitude (mas/yr)
               'GL':       {'ref': False, 'err': False, 'units': 'deg', 'format': 'f8'},       # Galactic longitude (degrees)
               'GB':       {'ref': False, 'err': False, 'units': 'deg', 'format': 'f8'},       # Galactic latitude (degrees)
               'RAJD':     {'ref': False, 'err': False, 'units': 'deg', 'format': 'f8'},       # Right ascension (J2000) (degrees)
               'DecJD':    {'ref': False, 'err': False, 'units': 'deg', 'format': 'f8'},       # Declination (J2000) (degrees)
               'Type':     {'ref': True, 'err': False, 'units': None, 'format': 'S32'},        # Type codes for the pulsar (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types)
               'Dist':     {'ref': False, 'err': False, 'units': 'kpc', 'format': 'f8'},       # Best estimate of the pulsar distance using the YMW16 DM-based distance as default (kpc)
               'Dist_DM':  {'ref': True, 'err': False, 'units': 'kpc', 'format': 'f8'},        # Distance based on the YMW16 electron density model. In 'LONG' or 'PUBLICATION QUALITY' modes, lower limits from the distance model are preceded by a '+' sign.
               'DMsinb':   {'ref': False, 'err': False, 'units': 'cm^-3 pc', 'format': 'f8'},  # DM x sin(b) (cm-3 pc)
               'ZZ':       {'ref': False, 'err': False, 'units': 'kpc', 'format': 'f8'},       # Distance from the Galactic plane, based on Dist
               'XX':       {'ref': False, 'err': False, 'units': 'kpc', 'format': 'f8'},       # X-Distance in X-Y-Z Galactic coordinate system (kpc)
               'YY':       {'ref': False, 'err': False, 'units': 'kpc', 'format': 'f8'},       # Y-Distance in X-Y-Z Galactic coordinate system (kpc)
               'Assoc':    {'ref': False, 'err': False, 'units': None, 'format': 'S32'},       # Names of other objects, e.g., supernova remnant, globular cluster or gamma-ray source associated with the pulsar
               'Survey':   {'ref': False, 'err': False, 'units': None, 'format': 'S32'},       # Surveys that detected the pulsar (discovery survey first) (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#surveys)
               'OSurvey':  {'ref': False, 'err': False, 'units': None, 'format': 'S32'},       # Surveys that detected the pulsar encoded as bits in integer
               'Date':     {'ref': False, 'err': False, 'units': 'yr', 'format': 'i4'},        # Date of discovery publication.
               'NGlt':     {'ref': False, 'err': False, 'units': None, 'format': 'i4'}         # Number of glitches observed for the pulsar
              }

# just return the parameter names
PSR_GENERAL_PARS = PSR_GENERAL.keys()

# timing solution and profile parameters
PSR_TIMING = {'P0':        {'ref': True, 'err': True, 'units': 's', 'format': 'f8'},           # Barycentric period of the pulsar (s)
              'P1':        {'ref': True, 'err': True, 'units': None, 'format': 'f8'},          # Time derivative of barcycentric period (dimensionless)
              'F0':        {'ref': True, 'err': True, 'units': 'Hz', 'format': 'f8'},          # Barycentric rotation frequency (Hz)
              'F1':        {'ref': True, 'err': True, 'units': 's^-2', 'format': 'f8'},        # Time derivative of barycentric rotation frequency (s-2)
              'F2':        {'ref': True, 'err': True, 'units': 's^-3', 'format': 'f8'},        # Second time derivative of barycentric rotation frequency (s-3)
              'F3':        {'ref': True, 'err': True, 'units': 's^-4', 'format': 'f8'},        # Third time derivative of barycentric rotation frequency (s-4)
              'PEpoch':    {'ref': True, 'err': False, 'units': 'd', 'format': 'f8'},          # Epoch of period or frequency (MJD)
              'DM':        {'ref': True, 'err': True, 'units': 'cm^-3 pc', 'format': 'f8'},    # Dispersion measure (cm-3 pc)
              'DM1':       {'ref': True, 'err': True, 'units': 'cm^-3 pc/yr', 'format': 'f8'}, # First time derivative of dispersion measure (cm-3 pc yr-1)
              'RM':        {'ref': True, 'err': True, 'units': 'rad m^-2', 'format': 'f8'},    # Rotation measure (rad m-2)
              'W50':       {'ref': True, 'err': True, 'units': 'ms', 'format': 'f8'},          # Width of pulse at 50% of peak (ms). Note, pulse widths are a function of both observing frequency and observational time resolution,so quoted widths are indicative only. Refer to the original reference for details.
              'W10':       {'ref': True, 'err': True, 'units': 'ms', 'format': 'f8'},          # Width of pulse at 10% (ms). Note the comments above for W50.
              'Units':     {'ref': True, 'err': False, 'units': None, 'format': 'S4'},         # Timescale for period/frequency and epoch data: TCB or TDB. See Hobbs, Edwards & Manchester (2006) for a discussion of the relationship between TCB and TDB.
              'Tau_sc':    {'ref': True, 'err': True, 'units': 's', 'format': 'f8'},           # Temporal broadening of pulses at 1 GHz due to interestellar scattering (s)
              'S400':      {'ref': True, 'err': True, 'units': 'mJy', 'format': 'f8'},         # Mean flux density at 400 MHz (mJy)
              'S1400':     {'ref': True, 'err': True, 'units': 'mJy', 'format': 'f8'},         # Mean flux density at 1400 MHz (mJy)
              'S2000':     {'ref': True, 'err': True, 'units': 'mJy', 'format': 'f8'}          # Mean flux density at 2000 MHz (mJy)
             }

PSR_TIMING_PARS = PSR_TIMING.keys()

# binary system parameters
PSR_BINARY = {'Binary':    {'ref': True, 'err': False, 'units': None, 'format': 'S5'},         # Binary model (usually one of several recognised by the pulsar timing programs TEMPO or TEMPO2). Modified versions of standard models are often used - refer to the source paper for details of the binary model used.
              'T0':        {'ref': True, 'err': True, 'units': 'd', 'format': 'f8'},           # Epoch of periastron (MJD)
              'PB':        {'ref': True, 'err': True, 'units': 'd', 'format': 'f8'},           # Binary period of pulsar (days)
              'A1':        {'ref': True, 'err': True, 'units': 's', 'format': 'f8'},           # Projected semi-major axis of orbit (lt s)
              'OM':        {'ref': True, 'err': True, 'units': 'deg', 'format': 'f8'},         # Longitude of periastron (degrees)
              'ECC':       {'ref': True, 'err': True, 'units': None, 'format': 'f8'},          # Eccentricity
              'TASC':      {'ref': True, 'err': True, 'units': 'd', 'format': 'f8'},           # Epoch of ascending node(MJD) - ELL1 binary model
              'EPS1':      {'ref': True, 'err': True, 'units': None, 'format': 'f8'},          # ECC x sin(OM) - ELL1 binary model
              'EPS2':      {'ref': True, 'err': True, 'units': None, 'format': 'f8'},          # ECC x cos(OM) - ELL1 binary model
              'MinMass':   {'ref': False, 'err': False, 'units': 'M_sun', 'format': 'f8'},     # Minimum companion mass assuming i=90 degrees and neutron star mass is 1.35 Mo
              'MedMass':   {'ref': False, 'err': False, 'units': 'M_sun', 'format': 'f8'},     # Median companion mass assuming i=60 degrees
              'BinComp':   {'ref': True, 'err': False, 'units': None, 'format': 'S4'}          # Companion type
             }

PSR_BINARY_PARS = PSR_BINARY.keys()

# derived parameters
PSR_DERIVED = {'R_Lum':    {'ref': False, 'err': False, 'units': 'mJy kpc^2', 'format': 'f8'},      # Radio luminosity at 400 MHz (mJy kpc2)
               'R_Lum14':  {'ref': False, 'err': False, 'units': 'mJy kpc^2', 'format': 'f8'},      # Radio luminosity at 1400 MHz (mJy kpc2)
               'Age':      {'ref': False, 'err': False, 'units': 'yr', 'format': 'f8'},             # Spin down age (yr) [tau = P/(2 Pdot)]
               'Bsurf':    {'ref': False, 'err': False, 'units': 'G', 'format': 'f8'},              # Surface magnetic flux density (Gauss) [B = 3.2e19 sqrt(P * Pdot)]
               'Edot':     {'ref': False, 'err': False, 'units': 'erg/s', 'format': 'f8'},          # Spin down energy loss rate (ergs/s)
               'Edotd2':   {'ref': False, 'err': False, 'units': 'erg s^-1/kpc^2', 'format': 'f8'}, # Energy flux at the Sun (ergs/kpc2/s)
               'PMTot':    {'ref': False, 'err': False, 'units': 'mas/yr', 'format': 'f8'},         # Total proper motion (mas/yr)
               'VTrans':   {'ref': False, 'err': False, 'units': 'km/s', 'format': 'f8'},           # Transverse velocity - based on DIST (km/s)
               'P1_i':     {'ref': False, 'err': False, 'units': None, 'format': 'f8'},             # Period derivative corrected for Shklovskii (proper motion) effect
               'Age_i':    {'ref': False, 'err': False, 'units': 'yr', 'format': 'f8'},             # Spin down age from P1_i (yr)
               'BSurf_i':  {'ref': False, 'err': False, 'units': 'G', 'format': 'f8'},              # Surface magnetic dipole from P1_i (gauss)
               'B_LC':     {'ref': False, 'err': False, 'units': 'G', 'format': 'f8'}               # Magnetic field at light cylinder (gauss)
              }

PSR_DERIVED_PARS = PSR_DERIVED.keys()

# a list of all allowed parameters for querying
PSR_ALL = dict(itertools.chain(PSR_GENERAL.items(), PSR_TIMING.items(), PSR_BINARY.items(), PSR_DERIVED.items()))
PSR_ALL_PARS = PSR_GENERAL_PARS + PSR_TIMING_PARS + PSR_BINARY_PARS + PSR_DERIVED_PARS

# "types" of pulsar http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types
# for use in 'type()' when setting logical conditions
PSR_TYPES = ['AXP',           # Anomalous X-ray Pulsar or Soft Gamma-ray Repeater with detected pulsations
             'BINARY',        # Pulsar has one or more stellar companion(s)
             'HE',            # Spin-powered pulsar with pulsed emission from radio to infrared or higher frequencies
             'NRAD',          # Spin-powered pulsar with pulsed emission only at infrared or higher frequencies
             'RADIO',         # Pulsars with pulsed emission in the radio band
             'RRAT',          # Pulsars with intermittently pulsed radio emission
             'XINS'           # Isolated neutron stars with pulsed thermal X-ray emission but no detectable radio emission
             ]

# binary companion types http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type
# for use in 'bincomp()' when setting logical conditions
PSR_BINARY_TYPE = ['MS',      # Main-sequence star
                   'NS',      # Neutron star
                   'CO',      # CO or ONeMg White Dwarf
                   'He',      # Helium White Dwarf
                   'UL'       # Ultra-light companion or planet (mass < 0.08 solar masses)
                  ]

# other objects associated with the pulsar (this is not an exhaustive list
# # for use in 'assoc()' when setting logical conditions)
PSR_ASSOC_TYPE = ['GC',  # globular cluster
                  'SNR'  # supernova remnant
                 ]

# URL for the NASA ADS
ADS_URL = 'https://ui.adsabs.harvard.edu/#abs/{}/'