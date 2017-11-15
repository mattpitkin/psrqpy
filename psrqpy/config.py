"""
Configuration information
"""

ATNF_VERSION = '1.56'
ATNF_URL = r'http://www.atnf.csiro.au/people/pulsar/psrcat'

# pulsar parameters (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html) that can be queried

# general parameters (e.g., name, position, distance...)
PSR_GENERAL = ['Name',        # Pulsar name.  The B name if exists, otherwise the J name.
               'JName',       # Pulsar name based on J2000 coordinates
               'RAJ',         # Right ascension (J2000) (hh:mm:ss.s)
               'DecJ',        # Declination (J2000) (+dd:mm:ss)
               'PMRA',        # Proper motion in the right ascension direction (mas/yr)
               'PMDec',       # Proper motion in declination (mas/yr)
               'PX',          # Annual parallax (mas)
               'PosEpoch',    # Epoch of position, defaults to PEpoch (MJD)
               'ELong',       # Ecliptic longitude (degrees)
               'ELat',        # Ecliptic latitude (degrees)
               'PMElong',     # Proper motion in the ecliptic longitude direction (mas/yr)
               'PMElat',      # Proper motion in ecliptic latitude (mas/yr)
               'GL',          # Galactic longitude (degrees)
               'GB',          # Galactic latitude (degrees)
               'RAJD',        # Right ascension (J2000) (degrees)
               'DecJD',       # Declination (J2000) (degrees)
               'Type',        # Type codes for the pulsar (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#psr_types)
               'Dist',        # Best estimate of the pulsar distance using the YMW16 DM-based distance as default (kpc)
               'Dist_DM',     # Distance based on the YMW16 electron density model. In 'LONG' or 'PUBLICATION QUALITY' modes, lower limits from the distance model are preceded by a '+' sign.
               'DMsinb',      # DM x sin(b) (cm-3 pc)
               'ZZ',          # Distance from the Galactic plane, based on Dist
               'XX',          # X-Distance in X-Y-Z Galactic coordinate system (kpc)
               'YY',          # Y-Distance in X-Y-Z Galactic coordinate system (kpc)
               'Assoc',       # Names of other objects, e.g., supernova remnant, globular cluster or gamma-ray source associated with the pulsar
               'Survey',      # Surveys that detected the pulsar (discovery survey first) (http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html#surveys)
               'OSurvey',     # Surveys that detected the pulsar encoded as bits in integer
               'Date',        # Date of discovery publication.
               'NGlt'         # Number of glitches observed for the pulsar
               ]

# timing solution and profile parameters
PSR_TIMING = ['P0',           # Barycentric period of the pulsar (s)
              'P1',           # Time derivative of barcycentric period (dimensionless)
              'F0',           # Barycentric rotation frequency (Hz)
              'F1',           # Time derivative of barycentric rotation frequency (s-2)
              'F2',           # Second time derivative of barycentric rotation frequency (s-3)
              'F3',           # Third time derivative of barycentric rotation frequency (s-4)
              'PEpoch',       # Epoch of period or frequency (MJD)
              'DM',           # Dispersion measure (cm-3 pc)
              'DM1',          # First time derivative of dispersion measure (cm-3 pc yr-1)
              'RM',           # Rotation measure (rad m-2)
              'W50',          # Width of pulse at 50% of peak (ms). Note, pulse widths are a function of both observing frequency and observational time resolution,so quoted widths are indicative only. Refer to the original reference for details.
              'W10',          # Width of pulse at 10% (ms). Note the comments above for W50.
              'Units',        # Timescale for period/frequency and epoch data: TCB or TDB. See Hobbs, Edwards & Manchester (2006) for a discussion of the relationship between TCB and TDB.
              'Tau_sc',       # Temporal broadening of pulses at 1 GHz due to interestellar scattering (s)
              'S400',         # Mean flux density at 400 MHz (mJy)
              'S1400',        # Mean flux density at 1400 MHz (mJy)
              'S2000'         # Mean flux density at 2000 MHz (mJy)
             ]

# binary system parameters
PSR_BINARY = ['Binary',       # Binary model (usually one of several recognised by the pulsar timing programs TEMPO or TEMPO2). Modified versions of standard models are often used - refer to the source paper for details of the binary model used.
              'T0',           # Epoch of periastron (MJD)
              'PB',           # Binary period of pulsar (days)
              'A1',           # Projected semi-major axis of orbit (lt s)
              'OM',           # Longitude of periastron (degrees)
              'ECC',          # Eccentricity
              'TASC',         # Epoch of ascending node(MJD) - ELL1 binary model
              'EPS1',         # ECC x sin(OM) - ELL1 binary model
              'EPS2',         # ECC x cos(OM) - ELL1 binary model
              'MinMass',      # Minimum companion mass assuming i=90 degrees and neutron star mass is 1.35 Mo
              'MedMass',      # Median companion mass assuming i=60 degrees
              'BinComp'       # Companion type
             ]

# derived parameters
PSR_DERIVED = ['R_Lum',       # Radio luminosity at 400 MHz (mJy kpc2)
               'R_Lum14',     # Radio luminosity at 1400 MHz (mJy kpc2)
               'Age',         # Spin down age (yr) [￼tau = P/(2 Pdot)]
               'Bsurf',       # Surface magnetic flux density (Gauss) [B = 3.2e19 sqrt(P * Pdot)]
               'Edot',        # Spin down energy loss rate (ergs/s)
               'Edotd2',      # Energy flux at the Sun (ergs/kpc2/s)
               'PMTot',       # Total proper motion (mas/yr)
               'VTrans',      # Transverse velocity - based on DIST (km/s)
               'P1_i',        # Period derivative corrected for Shklovskii (proper motion) effect
               'Age_i',       # Spin down age from P1_i (yr)
               'BSurf_i',     # Surface magnetic dipole from P1_i (gauss)
               'B_LC'         # Magnetic field at light cylinder
              ]

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