"""
This submodules sets up common constants for use, such as the allowed pulsar parameters and various
URLs used for queries.
"""

import itertools

ATNF_VERSION = '1.58'  #: the default ATNF catalogue version
ATNF_BASE_URL = r'http://www.atnf.csiro.au/people/pulsar/psrcat/'  #: the ATNF pulsar catalogue base URL
ATNF_URL = ATNF_BASE_URL + r'proc_form.php?version={version}'  #: the ATNF pulsar catalogue base URL for queries
ATNF_TARBALL = ATNF_BASE_URL + r'downloads/psrcat_pkg.tar.gz'  #: name of the tarball containing the entire catalogue database

GLITCH_URL = r'http://www.jb.man.ac.uk/pulsar/glitches/gTable.html'  #: Jodrell Bank glitch catalogue table URL 

PARAMS_QUERY = r'{params}'
USERDEFINED_QUERY = r'&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val='
CONDITION_QUERY = r'&condition={condition}'
PSRNAMES_QUERY = r'&pulsar_names={psrnames}'
SORT_QUERY = r'&sort_attr={sortattr}&sort_order={sortorder}'
EPHEMERIS_QUERY = r'&submit_ephemeris={getephemeris}'
COORD1_QUERY = r'&coords_1={coord1}'
COORD2_QUERY = r'&coords_2={coord2}'
RADIUS_QUERY = r'&radius={radius}' 
QUERY_FLUFF = r'&ephemeris=long&coords_unit=raj%2Fdecj&style=Long+with+errors&no_value=*&nohead=nohead&state=query&table_bottom.x=30&table_bottom.y=22'

#: the full ATNF catalogue query URL
QUERY_URL = ATNF_URL + PARAMS_QUERY + USERDEFINED_QUERY + SORT_QUERY + CONDITION_QUERY + PSRNAMES_QUERY + EPHEMERIS_QUERY + RADIUS_QUERY + COORD1_QUERY + COORD2_QUERY + QUERY_FLUFF

# pulsar parameters (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html) that can be
# queried. For each parameter there is a dictionary giving:
#  - 'ref': True if the parameter can have an associated reference in the ATNF catalogue
#  - 'err': True if the parameter can have an associated error value
#  - 'unit': a string giving the units for the parameter (to be used if generating an astropy table)
#  - 'format': a string giving the parameter format (to be used if generating an astropy table)
PSR_GENERAL = {'NAME':     {'ref': True,  'err': False, 'units': None,       'format': 'U32'},  # Pulsar name.  The B name if exists, otherwise the J name.
               'JNAME':    {'ref': True,  'err': False, 'units': None,       'format': 'U32'},  # Pulsar name based on J2000 coordinates
               'BNAME':    {'ref': True,  'err': False, 'units': None,       'format': 'U32'},  # Pulsar Besselian name
               'RAJ':      {'ref': True,  'err': True,  'units': None,       'format': 'U64'},  # Right ascension (J2000) (hh:mm:ss.s)
               'DECJ':     {'ref': True,  'err': True,  'units': None,       'format': 'U64'},  # Declination (J2000) (+dd:mm:ss)
               'PMRA':     {'ref': True,  'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in the right ascension direction (mas/yr)
               'PMDEC':    {'ref': True,  'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in declination (mas/yr)
               'PX':       {'ref': True,  'err': True,  'units': 'mas',      'format': 'f8'},   # Annual parallax (mas)
               'POSEPOCH': {'ref': True,  'err': False, 'units': 'd',        'format': 'f8'},   # Epoch of position, defaults to PEpoch (MJD)
               'ELONG':    {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},   # Ecliptic longitude (degrees)
               'ELAT':     {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},   # Ecliptic latitude (degrees)
               'PMELONG':  {'ref': True,  'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in the ecliptic longitude direction (mas/yr)
               'PMELAT':   {'ref': True,  'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in ecliptic latitude (mas/yr)
               'GL':       {'ref': False, 'err': False, 'units': 'deg',      'format': 'f8'},   # Galactic longitude (degrees)
               'GB':       {'ref': False, 'err': False, 'units': 'deg',      'format': 'f8'},   # Galactic latitude (degrees)
               'RAJD':     {'ref': False, 'err': False, 'units': 'deg',      'format': 'f8'},   # Right ascension (J2000) (degrees)
               'DECJD':    {'ref': False, 'err': False, 'units': 'deg',      'format': 'f8'},   # Declination (J2000) (degrees)
               'TYPE':     {'ref': True,  'err': False, 'units': None,       'format': 'U64'},  # Type codes for the pulsar (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types)
               'PML':      {'ref': False, 'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in Galactic longitude (mas/yr)
               'PMB':      {'ref': False, 'err': True,  'units': 'mas/yr',   'format': 'f8'},   # Proper motion in Galactic latitude (mas/yr)
               'DIST':     {'ref': False, 'err': False, 'units': 'kpc',      'format': 'f8'},   # Best estimate of the pulsar distance using the YMW16 DM-based distance as default (kpc)
               'DIST_DM':  {'ref': True,  'err': False, 'units': 'kpc',      'format': 'f8'},   # Distance based on the YMW16 electron density model. In 'LONG' or 'PUBLICATION QUALITY' modes, lower limits from the distance model are preceded by a '+' sign.
               'DIST_DM1': {'ref': True,  'err': False, 'units': 'kpc',      'format': 'f8'},   # Distance based on NE2001 model (kpc)
               'DIST1':    {'ref': False, 'err': False, 'units': 'kpc',      'format': 'f8'},   # Best estimate of the pulsar distance using the NE2001 DM-based distance as default (kpc)
               'DIST_AMN': {'ref': True,  'err': True,  'units': 'kpc',      'format': 'f8'},   # Lower limit on independent distance estimate (kpc) (NOTE: although there is an 'error' column it is always zero)
               'DIST_AMX': {'ref': True,  'err': False, 'units': 'kpc',      'format': 'f8'},   # Upper limit on independent distance estimate (kpc)
               'DIST_A':   {'ref': True,  'err': True,  'units': 'kpc',      'format': 'f8'},   # Independent distance estimate - takes precedence over other distance estimates for DIST and DIST1 (kpc)
               'DMSINB':   {'ref': False, 'err': False, 'units': 'cm^-3 pc', 'format': 'f8'},   # DM x sin(b) (cm-3 pc)
               'ZZ':       {'ref': False, 'err': False, 'units': 'kpc',      'format': 'f8'},   # Distance from the Galactic plane, based on Dist
               'XX':       {'ref': False, 'err': False, 'units': 'kpc',      'format': 'f8'},   # X-Distance in X-Y-Z Galactic coordinate system (kpc)
               'YY':       {'ref': False, 'err': False, 'units': 'kpc',      'format': 'f8'},   # Y-Distance in X-Y-Z Galactic coordinate system (kpc)
               'ASSOC':    {'ref': False, 'err': False, 'units': None,       'format': 'U128'}, # Names of other objects, e.g., supernova remnant, globular cluster or gamma-ray source associated with the pulsar
               'SURVEY':   {'ref': False, 'err': False, 'units': None,       'format': 'U64'},  # Surveys that detected the pulsar (discovery survey first) (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#surveys)
               'OSURVEY':  {'ref': False, 'err': False, 'units': None,       'format': 'U32'},  # Surveys that detected the pulsar encoded as bits in integer
               'DATE':     {'ref': False, 'err': False, 'units': 'yr',       'format': 'i4'},   # Date of discovery publication
               'NGLT':     {'ref': False, 'err': False, 'units': None,       'format': 'i4'},   # Number of glitches observed for the pulsar
               'GLEP':     {'ref': False, 'err': False, 'units': 'd',        'format': 'f8'},   # Epoch of glitch
               'GLPH':     {'ref': False, 'err': False, 'units': None,       'format': 'f8'},   # Phase increment at glitch
               'GLF0':     {'ref': False, 'err': False, 'units': 'Hz',       'format': 'f8'},   # Permanent pulse frequency increment at glitch
               'GLF1':     {'ref': False, 'err': False, 'units': 'Hz/s',     'format': 'f8'},   # Permanent frequency derivative increment at glitch
               'GLF0D':    {'ref': False, 'err': False, 'units': 'Hz',       'format': 'f8'},   # Decaying frequency increment at glitch
               'GLTD':     {'ref': False, 'err': False, 'units': 'd',        'format': 'f8'},   # Time constant for decaying frequency increment
               'CLK':      {'ref': True,  'err': False, 'units': None,       'format': 'U32'},  # Reference clock used for timing solution
               'EPHEM':    {'ref': True,  'err': False, 'units': None,       'format': 'U32'}  # Solar-system ephemeris used for timing solution
              }

# just return the parameter names
PSR_GENERAL_PARS = list(PSR_GENERAL.keys())

# timing solution and profile parameters
PSR_TIMING = {'P0':      {'ref': True,  'err': True,  'units': 's',              'format': 'f8'}, # Barycentric period of the pulsar (s)
              'P1':      {'ref': True,  'err': True,  'units': None,             'format': 'f8'}, # Time derivative of barcycentric period (dimensionless)
              'F0':      {'ref': True,  'err': True,  'units': 'Hz',             'format': 'f8'}, # Barycentric rotation frequency (Hz)
              'F1':      {'ref': True,  'err': True,  'units': 's^-2',           'format': 'f8'}, # Time derivative of barycentric rotation frequency (s-2)
              'F2':      {'ref': True,  'err': True,  'units': 's^-3',           'format': 'f8'}, # Second time derivative of barycentric rotation frequency (s-3)
              'F3':      {'ref': True,  'err': True,  'units': 's^-4',           'format': 'f8'}, # Third time derivative of barycentric rotation frequency (s-4)
              'F4':      {'ref': True,  'err': True,  'units': 's^-5',           'format': 'f8'}, # Fourth time derivative of barycentric rotation frequency (s^-5)
              'F5':      {'ref': True,  'err': True,  'units': 's^-6',           'format': 'f8'}, # Fifth time derivative of barycentric rotation frequency (s^-6)
              'F6':      {'ref': False, 'err': False, 'units': 's^-7',           'format': 'f8'}, # Sixth time derivative of barycentric rotation frequency (s^-7)
              'F7':      {'ref': False, 'err': False, 'units': 's^-8',           'format': 'f8'}, # Seventh time derivative of barycentric rotation frequency (s^-7)
              'F8':      {'ref': False, 'err': False, 'units': 's^-9',           'format': 'f8'}, # Eighth time derivative of barycentric rotation frequency (s^-7)
              'F9':      {'ref': False, 'err': False, 'units': 's^-10',          'format': 'f8'}, # Ninth time derivative of barycentric rotation frequency (s^-7)
              'PEPOCH':  {'ref': True,  'err': False, 'units': 'd',              'format': 'f8'}, # Epoch of period or frequency (MJD)
              'DM':      {'ref': True,  'err': True,  'units': 'cm^-3 pc',       'format': 'f8'}, # Dispersion measure (cm-3 pc)
              'DM1':     {'ref': True,  'err': True,  'units': 'cm^-3 pc/yr',    'format': 'f8'}, # First time derivative of dispersion measure (cm-3 pc yr-1)
              'DMEPOCH': {'ref': True,  'err': False, 'units': 'd',              'format': 'f8'}, # Reference epoch for DM, defaults to PEpoch (MJD)
              'DM2':     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-2', 'format': 'f8'}, # Second time derivative of dispersion measure (cm-3pc yr-2)
              'DM3':     {'ref': True,  'err': True,  'units': 'cm^-3 pc yr^-3', 'format': 'f8'}, # Third time derivative of dispersion measure (cm-3pc yr-3)
              'DM4':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-4', 'format': 'f8'}, # Fourth time derivative of dispersion measure (cm-3pc yr-4)
              'DM5':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-5', 'format': 'f8'}, # Fifth time derivative of dispersion measure (cm-3pc yr-5)
              'DM6':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-6', 'format': 'f8'}, # Sixth time derivative of dispersion measure (cm-3pc yr-6)
              'DM7':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-7', 'format': 'f8'}, # Seventh time derivative of dispersion measure (cm-3pc yr-7)
              'DM8':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-8', 'format': 'f8'}, # Eighth time derivative of dispersion measure (cm-3pc yr-8)
              'DM9':     {'ref': False, 'err': False, 'units': 'cm^-3 pc yr^-9', 'format': 'f8'}, # Ninth time derivative of dispersion measure (cm-3pc yr-9)
              'RM':      {'ref': True,  'err': True,  'units': 'rad m^-2',       'format': 'f8'}, # Rotation measure (rad m-2)
              'W50':     {'ref': True,  'err': True,  'units': 'ms',             'format': 'f8'}, # Width of pulse at 50% of peak (ms). Note, pulse widths are a function of both observing frequency and observational time resolution,so quoted widths are indicative only. Refer to the original reference for details.
              'W10':     {'ref': True,  'err': True,  'units': 'ms',             'format': 'f8'}, # Width of pulse at 10% (ms). Note the comments above for W50.
              'UNITS':   {'ref': True,  'err': False, 'units': None,             'format': 'U4'}, # Timescale for period/frequency and epoch data: TCB or TDB. See Hobbs, Edwards & Manchester (2006) for a discussion of the relationship between TCB and TDB.
              'TAU_SC':  {'ref': True,  'err': True,  'units': 's',              'format': 'f8'}, # Temporal broadening of pulses at 1 GHz due to interestellar scattering (s)
              'S400':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Mean flux density at 400 MHz (mJy)
              'S1400':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Mean flux density at 1400 MHz (mJy)
              'S2000':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Mean flux density at 2000 MHz (mJy)
              'S30':     {'ref': False, 'err': False, 'units': 'mJy',            'format': 'f8'}, # Flux at 30 MHz (mJy)
              'S40':     {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 40 MHz (mJy)
              'S50':     {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 50 MHz (mJy)
              'S60':     {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 60 MHz (mJy)
              'S80':     {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 80 MHz (mJy)
              'S100':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 100 MHz (mJy)
              'S200':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 200 MHz (mJy)
              'S300':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 300 MHz (mJy)
              'S600':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 600 MHz (mJy)
              'S700':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 700 MHz (mJy)
              'S800':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 800 MHz (mJy)
              'S900':    {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 900 MHz (mJy)
              'S1600':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 1600 MHz (mJy)
              'S3000':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 3000 MHz (mJy)
              'S4000':   {'ref': True,  'err': False, 'units': 'mJy',            'format': 'f8'}, # Flux at 4000 MHz (mJy)
              'S6000':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 6000 MHz (mJy)
              'S8000':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 8000 MHz (mJy)
              'S10G':    {'ref': False, 'err': False, 'units': 'mJy',            'format': 'f8'}, # Flux at 10 GHz (mJy)
              'S20G':    {'ref': False, 'err': False, 'units': 'mJy',            'format': 'f8'}, # Flux at 20 GHz (mJy)
              'S50G':    {'ref': False, 'err': False, 'units': 'mJy',            'format': 'f8'}, # Flux at 50 GHz (mJy)
              'S100G':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 100 GHz (mJy)
              'S150G':   {'ref': True,  'err': True,  'units': 'mJy',            'format': 'f8'}, # Flux at 150 GHz (mJy)
              'SPINDX':  {'ref': True,  'err': True,  'units': None,             'format': 'f8'}  # Radio spectral index
             }

PSR_TIMING_PARS = list(PSR_TIMING.keys())

# binary system parameters
PSR_BINARY = {'BINARY':   {'ref': True,  'err': False, 'units': None,       'format': 'U8'},  # Binary model (usually one of several recognised by the pulsar timing programs TEMPO or TEMPO2). Modified versions of standard models are often used - refer to the source paper for details of the binary model used.
              'T0':       {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Epoch of periastron (MJD)
              'PB':       {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Binary period of pulsar (days)
              'A1':       {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Projected semi-major axis of orbit (lt s)
              'OM':       {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Longitude of periastron (degrees)
              'ECC':      {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Eccentricity
              'TASC':     {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Epoch of ascending node(MJD) - ELL1 binary model
              'EPS1':     {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # ECC x sin(OM) - ELL1 binary model
              'EPS2':     {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # ECC x cos(OM) - ELL1 binary model
              'MINMASS':  {'ref': False, 'err': False, 'units': 'M_sun',    'format': 'f8'},  # Minimum companion mass assuming i=90 degrees and neutron star mass is 1.35 Mo
              'MEDMASS':  {'ref': False, 'err': False, 'units': 'M_sun',    'format': 'f8'},  # Median companion mass assuming i=60 degrees
              'BINCOMP':  {'ref': True,  'err': False, 'units': None,       'format': 'U32'}, # Companion type
              'FB0':      {'ref': True,  'err': True,  'units': 'Hz',       'format': 'f8'},  # Orbital frequency (Hz)
              'FB1':      {'ref': True,  'err': True,  'units': 'Hz s^-1',  'format': 'f8'},  # 1st time derivative of orbital frequency (Hz s^-1)
              'FB2':      {'ref': True,  'err': True,  'units': 'Hz s^-2',  'format': 'f8'},  # 2nd time derivative of orbital frequency (Hz s^-2)
              'OMDOT':    {'ref': True,  'err': True,  'units': 'deg/yr',   'format': 'f8'},  # 1st time derivative of periastron longitude (periastron advance) (deg yr^-1)
              'OM2DOT':   {'ref': True,  'err': True,  'units': 'deg/yr^2', 'format': 'f8'},  # 2nd time derivative of periastron longitude (deg yr^-2)
              'A1DOT':    {'ref': True,  'err': True,  'units': 's s^-1',   'format': 'f8'},  # 1st time derivative of projected semi-major axis (lt-s s^-1)
              'A12DOT':   {'ref': True,  'err': True,  'units': 's s^-2',   'format': 'f8'},  # 2nd time derivative of projected semi-major axis (lt-s s^-2)
              'ECCDOT':   {'ref': True,  'err': True,  'units': 's^-1',     'format': 'f8'},  # 1st time derivative of eccentricity (s^-1)
              'ECC2DOT':  {'ref': False, 'err': False, 'units': 's^-2',     'format': 'f8'},  # 2nd time derivative of eccentricity (s^-2)
              'PBDOT':    {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # 1st time derivative of binary period (dimensionless)
              'GAMMA':    {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Post-Keplerian time-dilation term (s)
              'T0_2':     {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Epoch of periastron [2nd orbit (where innermost orbit is 1st)] (MJD)
              'PB_2':     {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Binary period of pulsar [2nd orbit] (days)
              'A1_2':     {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Projected semi-major axis of orbit [2nd orbit] (s)
              'OM_2':     {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Longitude of periastron [2nd orbit] (deg)
              'ECC_2':    {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Eccentricity [2nd orbit]
              'OMDOT_2':  {'ref': False, 'err': False, 'units': 'deg/yr',   'format': 'f8'},  # Periastron advance [2nd orbit] (deg/yr) 
              'PBDOT_2':  {'ref': False, 'err': False, 'units': None,       'format': 'f8'},  # 1st time derivative of binary period [2nd orbit]
              'EPS1_2':   {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # ECC_2 x sin(OM_2) [2nd orbit]
              'EPS2_2':   {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # ECC_2 x cos(OM_2) [2nd orbit]
              'TASC_2':   {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Epoch of ascending node (MJD) [2nd orbit] 
              'T0_3':     {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Epoch of periastron [3rd orbit] (MJD)
              'PB_3':     {'ref': True,  'err': True,  'units': 'd',        'format': 'f8'},  # Binary period of pulsar [3rd orbit] (days)
              'A1_3':     {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Projected semi-major axis of orbit [3rd orbit] (s)
              'OM_3':     {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Longitude of periastron [3rd orbit] (deg)
              'ECC_3':    {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Eccentricity [3rd orbit]
              'OMDOT_3':  {'ref': False, 'err': False, 'units': 'deg/yr',   'format': 'f8'},  # Periastron advance [3rd orbit] (deg/yr)
              'PBDOT_3':  {'ref': False, 'err': False, 'units': None,       'format': 'f8'},  # 1st time derivative of binary period [3rd orbit]
              'PPNGAMMA': {'ref': False, 'err': False, 'units': 's',        'format': 'f8'},  # PPN parameter gamma (s)
              'SINI':     {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Sine of inclination angle
              'SINI_2':   {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Sine of inclination angle [2nd orbit]
              'SINI_3':   {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Sine of inclination angle [3rd orbit]
              'XPBDOT':   {'ref': False, 'err': False, 'units': None,       'format': 'f8'},  # Rate of change of orbital period minus GR prediction
              'KOM':      {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Long. on sky of asc. node (N toward E) from ann. orbital parallax  (deg)
              'KIN':      {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Orbit inclination from annual orbital parallax  (deg)
              'M2':       {'ref': True,  'err': True,  'units': 'M_sun',    'format': 'f8'},  # Companion mass (solar masses)
              'M2_2':     {'ref': True,  'err': True,  'units': 'M_sun',    'format': 'f8'},  # Companion mass [2nd orbit] (solar masses)
              'M2_3':     {'ref': True,  'err': True,  'units': 'M_sun',    'format': 'f8'},  # Companion mass [3rd orbit] (solar masses)
              'MASS_Q':   {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Mass ratio for binary: M1/M2
              'MASS_Q_2': {'ref': False, 'err': False, 'units': None,       'format': 'f8'},  # Mass ratio for 2nd binary companion: M1/M2_2
              'OM_ASC':   {'ref': True,  'err': True,  'units': 'deg',      'format': 'f8'},  # Longitude on sky of ascending node (from N toward E) (deg)
              'OM_ASC_2': {'ref': False, 'err': False, 'units': 'deg',      'format': 'f8'},  # Longitude on sky of ascending node (2nd orbit) (deg)
              'DTHETA':   {'ref': True,  'err': True,  'units': '1e-6',     'format': 'f8'},  # Relativistic deformation of the orbit (1e-6)
              'XOMDOT':   {'ref': False, 'err': False, 'units': 'deg/yr',   'format': 'f8'},  # Rate of periastron advance minus GR prediction (deg/yr)
              'H3':       {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Amplitude of 3rd Shapiro-delay harmonic (s)
              'H4':       {'ref': True,  'err': True,  'units': 's',        'format': 'f8'},  # Amplitude of 4th Shapiro-delay harmonic (s)
              'STIG':     {'ref': True,  'err': True,  'units': None,       'format': 'f8'},  # Ratio of successive Shapiro-delay harmonics (e.g., H4/H3)
              'MASSFN':   {'ref': False, 'err': False, 'units': 'M_sun',    'format': 'f8'},  # The pulsar mass function (solar mass)
              'UPRMASS':  {'ref': False, 'err': False, 'units': 'M_sun',    'format': 'f8'},  # 90% confidence upper companion mass limit, i=26deg (solar mass)
              'MINOMDOT': {'ref': False, 'err': False, 'units': 'deg/yr',   'format': 'f8'}   # Minimum omega dot, assuming sin i = 1 and M_ns = 1.4Mo (deg/yr)
             }

PSR_BINARY_PARS = list(PSR_BINARY.keys())

# derived parameters
PSR_DERIVED = {'R_LUM':    {'ref': False, 'err': False, 'units': 'mJy kpc^2',      'format': 'f8'}, # Radio luminosity at 400 MHz (mJy kpc2)
               'R_LUM14':  {'ref': False, 'err': False, 'units': 'mJy kpc^2',      'format': 'f8'}, # Radio luminosity at 1400 MHz (mJy kpc2)
               'AGE':      {'ref': False, 'err': False, 'units': 'yr',             'format': 'f8'}, # Spin down age (yr) [tau = P/(2 Pdot)]
               'BSURF':    {'ref': False, 'err': False, 'units': 'G',              'format': 'f8'}, # Surface magnetic flux density (Gauss) [B = 3.2e19 sqrt(P * Pdot)]
               'EDOT':     {'ref': False, 'err': False, 'units': 'erg/s',          'format': 'f8'}, # Spin down energy loss rate (ergs/s)
               'EDOTD2':   {'ref': False, 'err': False, 'units': 'erg s^-1/kpc^2', 'format': 'f8'}, # Energy flux at the Sun (ergs/kpc2/s)
               'PMTOT':    {'ref': False, 'err': True,  'units': 'mas/yr',         'format': 'f8'}, # Total proper motion (mas/yr)
               'VTRANS':   {'ref': False, 'err': False, 'units': 'km/s',           'format': 'f8'}, # Transverse velocity - based on DIST (km/s)
               'P1_I':     {'ref': False, 'err': False, 'units': None,             'format': 'f8'}, # Period derivative corrected for Shklovskii (proper motion) effect
               'AGE_I':    {'ref': False, 'err': False, 'units': 'yr',             'format': 'f8'}, # Spin down age from P1_i (yr)
               'BSURF_I':  {'ref': False, 'err': False, 'units': 'G',              'format': 'f8'}, # Surface magnetic dipole from P1_i (gauss)
               'B_LC':     {'ref': False, 'err': False, 'units': 'G',              'format': 'f8'}  # Magnetic field at light cylinder (gauss)
              }

PSR_DERIVED_PARS = list(PSR_DERIVED.keys())


# a list of all allowed parameters for querying
PSR_ALL = dict(itertools.chain(PSR_GENERAL.items(), PSR_TIMING.items(), PSR_BINARY.items(), PSR_DERIVED.items()))
""": a dict of allowed pulsars parameters (e.g., name, position, distance...)

Each parameter name key gives a dictionary containing the keys:

* ``ref`` (bool) - True if the parameter has an associated reference with it
* ``err`` (bool) - True if the parameter has an associated error value
* ``units`` (str) - a string with the parameters units that can be parsed by
  :class:`~astropy.units.core.Unit`
* ``format`` (str) - a string with a :class:`numpy.dtype` for storing the parameter in an
  :class:`~astropy.table.Table`

The allowed parameters and their units are given
`here <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#par_list>`_.
"""

PSR_ALL_PARS = PSR_GENERAL_PARS + PSR_TIMING_PARS + PSR_BINARY_PARS + PSR_DERIVED_PARS

PSR_TYPES = ['AXP',           # Anomalous X-ray Pulsar or Soft Gamma-ray Repeater with detected pulsations
             'BINARY',        # Pulsar has one or more stellar companion(s)
             'HE',            # Spin-powered pulsar with pulsed emission from radio to infrared or higher frequencies
             'NRAD',          # Spin-powered pulsar with pulsed emission only at infrared or higher frequencies
             'RADIO',         # Pulsars with pulsed emission in the radio band
             'RRAT',          # Pulsars with intermittently pulsed radio emission
             'XINS'           # Isolated neutron stars with pulsed thermal X-ray emission but no detectable radio emission
             ]
""": `types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types>`_ of
pulsar for use in ``type()`` when setting logical conditions.
"""

PSR_BINARY_TYPE = ['MS',      # Main-sequence star
                   'NS',      # Neutron star
                   'CO',      # CO or ONeMg White Dwarf
                   'He',      # Helium White Dwarf
                   'UL'       # Ultra-light companion or planet (mass < 0.08 solar masses)
                  ]
""": binary companion
`types <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#bincomp_type>`_
for use in ``bincomp()`` when setting logical conditions.
"""

#: other objects associated with the pulsar (this is not an exhaustive list for use in ``assoc()`` when setting logical conditions)
PSR_ASSOC_TYPE = ['GC',  # globular cluster
                  'SNR'  # supernova remnant
                 ]

#: URL for the NASA ADS
ADS_URL = 'https://ui.adsabs.harvard.edu/#abs/{}/'
