"""
Test script.
"""

import pytest
import os
from psrqpy import QueryATNF
import numpy as np
from pandas import Series
import pytest_socket
from astropy.table.column import MaskedColumn


def sf_scale(value):
    """
    Calculate the base-10 scale of the final significant figure for a given
    number.

    E.g. for a value of 12000.0 you would get 1000.0 as the scale of the final
    significant figure. Or, for 1.2345e-8 you would get 0.0001.
    """

    # base-10 exponent to which the value is raise
    valexp = np.floor(np.log10(np.abs(value)))

    # value that is raise to 10**valexp
    val = (value/10**valexp).astype('float32')  # type as float32 to avoid numerical noise

    valstr = str(val)

    # get the number of decimal places
    numdp = len(valstr) - valstr.find('.') - 1

    if valstr[-1] == '0':
        numdp -= 1

    # return the scale of the final significant figure
    return 10**(valexp - numdp)


def round_err(errvalue, atnferrvalue):
    """
    Round the derived error to the same number of significant figures at the
    error produced by `psrcat` and used for the ATNF pulsar catalogue, noting
    that `psrcat` rounds error up. Return True if the derived error and
    equivalent ATNF are the same
    """

    # get ATNF derived error value
    errval = (atnferrvalue/sf_scale(atnferrvalue)).astype('float32')

    # ATNF derived errors are always rounded up
    derval = np.ceil(errvalue/sf_scale(atnferrvalue))

    return derval == errval


def test_crab(query):
    """
    Test that the Crab pulsar is present and the frequency is as expected, i.e.
    the frequency rounds down to 29 Hz (should be OK for another ~80 years!)
    """

    f0 = query.get_pulsar('J0534+2200')['F0'][0]

    assert np.floor(f0) == 29.0

    # try Crab's B-name
    f0B = query.get_pulsar('B0531+21')['F0'][0]

    assert f0 == f0B

    # check reference and error are not None
    assert query.get_pulsar('B0531+21')['F0_ERR'][0] is not None
    assert query.get_pulsar('B0531+21')['F0_REF'][0] is not None


def test_catalogue_shape(query):
    """
    Test the catalogue for shape consistency
    """

    length = query.catalogue_len
    shape = query.catalogue_shape
    rows = query.catalogue_nrows
    cols = query.catalogue_ncols
    colnames = query.columns

    assert length == rows and length == shape[0]
    assert cols == len(colnames) and cols == shape[1]


def test_get_pulsars(query):
    """
    Test the 'Pulsars' class.
    """

    psrs = query.get_pulsars()

    assert len(psrs) == query.num_pulsars

    # check Crab frequency
    f01 = query.get_pulsar('J0534+2200')['F0'][0]
    f02 = psrs['J0534+2200'].F0  # frequency attribute

    assert f01 == f02

    # test removing a pulsar
    crab = psrs.pop('J0534+2200')

    f03 = crab.F0
    assert f03 == f01
    assert len(psrs) == (query.num_pulsars - 1)

    # get the ephemeris string for the Crab pulsar
    crabeph = query.get_ephemeris('J0534+2200AB')  # wrong name
    assert crabeph is None

    crabeph = query.get_ephemeris('J0534+2200')

    assert isinstance(crabeph, str)

    for line in crabeph.split('\n'):
        if line.split()[0].strip() == 'F0':
            f0str = line.split()[1].strip()
            break

    assert f01 == float(f0str)


def test_getitem(query):
    """
    Test using the __getitem__ method.
    """

    with pytest.raises(KeyError):
        query["BLah"]

    # get a pulsar
    crabeph = query["J0534+2200"]

    assert np.floor(crabeph["F0"][0]) == 29.0

    # get a column
    psrjs = query["PSRJ"]

    assert "J0534+2200" in psrjs


def test_save_load_file(tmp_path, query):
    """
    Test saving and reloading a query as a pickle file.
    """

    # test exception handling
    testfilebad = '/jkshfdjfd/jkgsdfjkj/kgskfd.jhfd'

    with pytest.raises(IOError):
        query.save(testfilebad)

    # test exception handling
    with pytest.raises(IOError):
        querynew = QueryATNF(loadquery=testfilebad)

    testfile = tmp_path / 'query.pkl'
    query.save(testfile)

    # re-load in as a new query
    querynew = QueryATNF(loadquery=testfile)

    assert query.num_pulsars == querynew.num_pulsars


def test_condition(query):
    """
    Test the parsing of logical conditions.
    """

    with pytest.raises(TypeError):
        # test for error if condition is not a string
        query.condition = 2.3

    # test that we only return pulsars with F0 > 100 Hz
    query.condition = 'F0 > 100'

    psrs = query.table
    f0s = psrs['F0']

    assert not np.any(f0s < 100.)

    # test that we only return pulsars with F0 > 100 Hz in binary systems
    query.condition = 'F0 > 100 && type(binary)'

    psrs = query.table
    f0s = psrs['F0']
    binary = psrs['BINARY']

    if type(binary) == MaskedColumn:
        assert not np.any(f0s < 100.) and not np.any(binary.mask)
    else:
        assert not np.any(f0s < 100.)

    # test 'OR'
    query.condition = 'F0 > 100 || type(binary)'

    psrs = query.table
    f0s = psrs['F0']
    binary = psrs['BINARY']

    if type(binary) == MaskedColumn:
        assert np.all(~binary[f0s < 100.].mask)

    # test 'NOT'
    query.condition = 'F0 > 100 and not type(binary)'

    psrs = query.table
    f0s = psrs['F0']
    binary = psrs['BINARY']

    if type(binary) == MaskedColumn:
        assert not np.any(f0s < 100) and np.all(binary.mask)
    else:
        assert not np.any(f0s < 100)

    # test type (not binary)
    query.condition = 'type(HE)'

    psrs = query.table
    types = psrs['TYPE']

    assert np.all(['HE' in ttype for ttype in types])

    # test associations
    query.condition = 'assoc(GC)'

    psrs = query.table
    assocs = psrs['ASSOC']

    assert np.all(['GC' in assoc for assoc in assocs])

    # test binary companion
    query.condition = 'bincomp(MS)'

    psrs = query.table
    bincomps = psrs['BINCOMP']

    assert np.all(['MS' in bincomp for bincomp in bincomps])

    # test exists
    query.condition = None
    allpsrs = query.table
    query.condition = 'exist(PMRA)'

    psrs = query.table
    pmras = psrs['PMRA']

    if type(allpsrs['PMRA']) == MaskedColumn:
        assert len(pmras) == np.sum(~allpsrs['PMRA'].mask)

    # test survey
    query.condition = "survey(gb350)"
    psrs = query.table
    surveys = psrs["SURVEY"]

    assert np.all(["gb350" in survey for survey in surveys])

    # test discovery
    query.condition = "discovery(htru_eff)"
    psrs = query.table
    discoveries = psrs["SURVEY"]

    assert np.all(["htru_eff" == disc.split(",")[0] for disc in discoveries])

    # reset condition
    query.condition = None


def test_num_pulsars(query):
    """
    Test that the number of pulsars returned is as expected.
    """

    query.psrs = 'J9999+9999'  # bad pulsar

    # length should be zero
    assert len(query) == 0

    query.psrs = 'J0534+2200'  # Crab pulsar

    # length should be one
    assert len(query) == 1

    query.psrs = ['J0534+2200', 'J0537-6910']

    # length should be two
    assert len(query) == 2


def test_num_columns(query):
    """
    Test that the number of columns if correct.
    """

    query.query_params = []
    assert len(query.query_params) == 0

    with pytest.raises(TypeError):
        # test error for non-string parameter
        query.query_params = 1.2

    with pytest.raises(TypeError):
        # test error for non-string parameter in list
        query.query_params = ['F0', 1.3]

    query.query_params = 'F0'
    query.include_errs = False

    # number of columns should be 1
    assert len(query.table.columns) == 1

    # name of column should be 'F0'
    assert query.table.keys()[0] == 'F0'

    query.include_errs = True

    # number of columns should be 2
    assert len(query.table.columns) == 2
    assert 'F0' in query.table.keys() and 'F0_ERR' in query.table.keys()

    query.query_params = ['F0', 'F1']

    # number of columns should be 4
    assert len(query.table.columns) == 4


def test_get_references(query):
    """
    Test getting the references (without wanting ADS urls).
    """

    query.get_references()

    # test parsing a reference
    ref = query.parse_ref('ksm+06')

    assert isinstance(ref, str)
    assert "Kramer" in ref and "Science" in ref and "314" in ref and "2006" in ref

    # test parsing two references
    ref = query.parse_ref(['ksm+06', 'abb+18'])

    assert len(ref) == 2
    assert "Kramer" in ref[0] and "Science" in ref[0] and "314" in ref[0] and "2006" in ref[0]
    assert "2018" in ref[1] and "ApJS" in ref[1] and "235" in ref[1] and"37" in ref[1]

    # test parsing two references (with the second one being gibberish)
    ref = query.parse_ref(['ksm+06', 'wlihlacljkblf'])

    assert len(ref) == 2
    assert "Kramer" in ref[0] and "Science" in ref[0] and "314" in ref[0] and "2006" in ref[0]
    assert ref[1] is None


def test_update(query):
    """
    Test update method.
    """

    # no name set (and column not a Series)
    column = 1
    with pytest.raises(ValueError):
        query.update(column)

    with pytest.raises(ValueError):
        query.update('dummystring', name='F0')

    # add an additional column using a numpy array
    newcol = np.ones(query.catalogue_len)
    newname = "TEST1"
    numcols = len(query.columns)

    assert newname not in query.columns

    query.update(newcol, name=newname)

    assert newname in query.columns
    assert len(query.columns) == (numcols + 1)
    assert np.all(query.catalogue[newname] == 1.)

    # add an additional columns using a Series
    numcols = len(query.columns)
    newname = "TEST2"
    newseries = Series(np.full(query.catalogue_len, 0.5), name=newname)

    assert newname not in query.columns

    query.update(newseries)

    assert newname in query.columns
    assert len(query.columns) == (numcols + 1)
    assert np.all(query.catalogue[newname] == 0.5)

    # make sure only NaNs get updated
    newname = "TEST3"
    newseries = Series(np.full(query.catalogue_len, np.nan), name=newname)
    newseries[0] = 1.
    query.update(newseries)

    assert query.catalogue[newname][0] == 1
    assert np.all(np.isnan(query.catalogue[newname][1:]))

    newseries = np.full(query.catalogue_len, 2.)
    query.update(newseries, name=newname)

    assert query.catalogue[newname][0] == 1
    assert np.all(query.catalogue[newname][1:] == 2.)

    # check that overwrite works
    query.update(newseries, name=newname, overwrite=True)

    assert np.all(query.catalogue[newname] == 2.)


def test_ppdot_diagram(query):
    """
    Test the creation of the P-Pdot diagram
    """

    from matplotlib.figure import Figure

    fig = query.ppdot(showtypes='ALL', showGCs=True, intrinsicpdot=True)
    assert isinstance(fig, Figure)


def test_glitch_table():
    """
    Try downloading the glitch table for the Crab.
    """

    from psrqpy.utils import get_glitch_catalogue

    table = get_glitch_catalogue(psr='J0534+2200')

    assert len(table) > 1


def test_gc_table():
    """
    Try downloading the globular cluster pulsar table.
    """

    from psrqpy.utils import get_gc_catalogue

    table = get_gc_catalogue()

    assert "47 Tuc" in table.clusters()
    assert "J0023-7204C" in table.cluster_pulsars("47 Tuc")["Pulsar"]
    assert "J0023-7204C" in table.cluster_pulsars("NGC 104")["Pulsar"]
    assert "J0023-7204C" in table.cluster_pulsars("104")["Pulsar"]
    assert "J0023-7204C" in table.cluster_pulsars("47Tuc")["Pulsar"]


def test_msp_table(query):
    """
    Try downloading the MSP table.
    """

    from psrqpy.utils import get_msp_catalogue

    table = get_msp_catalogue()

    # when added there were 400 pulsars in the table
    assert len(table) >= 400
    psrcheck = "J0023+0923"
    assert psrcheck in table["NAME"]

    # compare values to ATNF values to with 1%
    psr = query.get_pulsar(psrcheck)
    msppsr = table[table["NAME"] == psrcheck]

    assert 0.99 < psr["GB"] / msppsr["GB"] < 1.01
    assert 0.99 < psr["GL"] / msppsr["GL"] < 1.01
    assert 0.99 < psr["DM"] / msppsr["DM"] < 1.01
    assert 0.99 < psr["A1"] / msppsr["A1"] < 1.01
    assert 0.99 < psr["PB"] / msppsr["PB"] < 1.01


# TEST DERIVED PARAMETERS #
def test_derived_p0_p1(query_derived, query_atnf):
    """
    Test the derived period and period derivative values against the values
    from the ATNF Pulsar Catalogue.
    """

    p0 = query_derived.get_pulsar('TEST1')['P0'][0]
    p0atnf = query_atnf.get_pulsar('TEST1')['P0'][0]
    p0err = query_derived.get_pulsar('TEST1')['P0_ERR'][0]
    p0atnferr = query_atnf.get_pulsar('TEST1')['P0_ERR'][0]

    assert abs(p0-p0atnf) < sf_scale(p0atnf)
    assert round_err(p0err, p0atnferr)

    p1 = query_derived.get_pulsar('TEST1')['P1'][0]
    p1atnf = query_atnf.get_pulsar('TEST1')['P1'][0]
    p1err = query_derived.get_pulsar('TEST1')['P1_ERR'][0]
    p1atnferr = query_atnf.get_pulsar('TEST1')['P1_ERR'][0]

    assert abs(p1-p1atnf) < sf_scale(p1atnf)
    assert round_err(p1err, p1atnferr)


def test_derived_f0_f1(query_derived, query_atnf):
    """
    Test the derived frequency and frequency derivative values against the
    values from the ATNF Pulsar Catalogue.
    """

    f0 = query_derived.get_pulsar('TEST2')['F0'][0]
    f0atnf = query_atnf.get_pulsar('TEST2')['F0'][0]
    f0err = query_derived.get_pulsar('TEST2')['F0_ERR'][0]
    f0atnferr = query_atnf.get_pulsar('TEST2')['F0_ERR'][0]

    assert abs(f0-f0atnf) < sf_scale(f0atnf)
    assert round_err(f0err, f0atnferr)

    f1 = query_derived.get_pulsar('TEST2')['F1'][0]
    f1atnf = query_atnf.get_pulsar('TEST2')['F1'][0]
    f1err = query_derived.get_pulsar('TEST2')['F1_ERR'][0]
    f1atnferr = query_atnf.get_pulsar('TEST2')['F1_ERR'][0]

    assert abs(f1-f1atnf) < sf_scale(f1atnf)
    assert round_err(f1err, f1atnferr)


def test_derived_galactic(query_derived, query_atnf):
    """
    Test the derived galactic longitude and latitude. The galactic
    cartesian coordinates are defined differently between astropy
    and psrcat, so these cannot be directly compared.
    """

    gb = query_derived.get_pulsar('TEST1')['GB'][0]
    gbatnf = query_atnf.get_pulsar('TEST1')['GB'][0]

    gl = query_derived.get_pulsar('TEST1')['GL'][0]
    glatnf = query_atnf.get_pulsar('TEST1')['GL'][0]

    assert abs(gb-gbatnf) < sf_scale(gbatnf)
    assert abs(gl-glatnf) < sf_scale(glatnf)

    dmsinb = query_derived.get_pulsar('TEST1')['DMSINB'][0]
    dmsinbatnf = query_atnf.get_pulsar('TEST1')['DMSINB'][0]

    assert abs(dmsinb-dmsinbatnf) < sf_scale(dmsinbatnf)


def test_derived_ecliptic(query_derived, query_atnf):
    """
    Test the derived ecliptic longitude and latitude.
    """

    elong = query_derived.get_pulsar('TEST1')['ELONG'][0]
    elongatnf = query_atnf.get_pulsar('TEST1')['ELONG'][0]

    elat = query_derived.get_pulsar('TEST1')['ELAT'][0]
    elatatnf = query_atnf.get_pulsar('TEST1')['ELAT'][0]

    assert abs(elong-elongatnf) < 1e-5
    assert abs(elat-elatatnf) < 1e-5


def test_derived_binary_mass(query_derived, query_atnf):
    """
    Test the derived binary mass function and associated quantities.
    """

    # mass function
    massfn = query_derived.get_pulsar('TEST1')['MASSFN'][0]
    massfnatnf = query_atnf.get_pulsar('TEST1')['MASSFN'][0]

    # minimum companion mass
    minmass = query_derived.get_pulsar('TEST1')['MINMASS'][0]
    minmassatnf = query_atnf.get_pulsar('TEST1')['MINMASS'][0]

    # median companion mass
    medmass = query_derived.get_pulsar('TEST1')['MEDMASS'][0]
    medmassatnf = query_atnf.get_pulsar('TEST1')['MEDMASS'][0]

    # 90% upper bound in companion mass
    uprmass = query_derived.get_pulsar('TEST1')['UPRMASS'][0]
    uprmassatnf = query_atnf.get_pulsar('TEST1')['UPRMASS'][0]

    # minimum change in angle of peristron
    minomdot = query_derived.get_pulsar('TEST1')['MINOMDOT'][0]
    minomdotatnf = query_atnf.get_pulsar('TEST1')['MINOMDOT'][0]

    assert abs(massfn - massfnatnf) < sf_scale(massfnatnf)
    assert abs(minmass - minmassatnf) < sf_scale(minmassatnf)
    assert abs(medmass - medmassatnf) < sf_scale(medmassatnf)
    assert abs(uprmass - uprmassatnf) < sf_scale(uprmassatnf)
    assert abs(minomdot - minomdotatnf) < sf_scale(minomdotatnf)


def test_derived_binary_om_ecc(query_derived, query_atnf):
    """
    Test the values of the binary system angle of periastron and
    eccentricity derived from EPS1 and EPS2.
    """

    # eccentricity
    ecc = query_derived.get_pulsar('TEST1')['ECC'][0]
    eccatnf = query_atnf.get_pulsar('TEST1')['ECC'][0]
    eccerr = query_derived.get_pulsar('TEST1')['ECC_ERR'][0]
    eccerratnf = query_atnf.get_pulsar('TEST1')['ECC_ERR'][0]

    assert abs(ecc - eccatnf) < sf_scale(eccatnf)
    assert round_err(eccerr, eccerratnf)

    # angle of periastron
    om = query_derived.get_pulsar('TEST1')['OM'][0]
    omatnf = query_atnf.get_pulsar('TEST1')['OM'][0]
    omerr = query_derived.get_pulsar('TEST1')['OM_ERR'][0]
    omerratnf = query_atnf.get_pulsar('TEST1')['OM_ERR'][0]

    assert abs(om - omatnf) < sf_scale(omatnf)
    assert round_err(omerr, omerratnf)


def test_derived_age(query_derived, query_atnf):
    """
    Test the derived characteristic age.
    """

    age = query_derived.get_pulsar('TEST1')['AGE'][0]
    ageatnf = query_atnf.get_pulsar('TEST1')['AGE'][0]

    assert abs(age - ageatnf) < sf_scale(ageatnf)


def test_derived_bsurf(query_derived, query_atnf):
    """
    Test the derived surface magnetic field.
    """

    bsurf = query_derived.get_pulsar('TEST1')['BSURF'][0]
    bsurfatnf = query_atnf.get_pulsar('TEST1')['BSURF'][0]

    assert abs(bsurf - bsurfatnf) < sf_scale(bsurfatnf)


def test_instrinsic_pdot(query_derived, query_atnf):
    """
    Test parameters derived from the instrinsic period.
    """

    # instrinsic period derivative
    p1i = query_derived.get_pulsar('TEST1')['P1_I'][0]
    p1iatnf = query_atnf.get_pulsar('TEST1')['P1_I'][0]

    assert abs(p1i - p1iatnf) < sf_scale(p1iatnf)

    # instrinsic charateristic age
    agei = query_derived.get_pulsar('TEST1')['AGE_I'][0]
    ageiatnf = query_atnf.get_pulsar('TEST1')['AGE_I'][0]

    assert abs(agei - ageiatnf) < sf_scale(ageiatnf)

    # intrinsic surface magnetic field strength
    bsurfi = query_derived.get_pulsar('TEST1')['BSURF_I'][0]
    bsurfiatnf = query_atnf.get_pulsar('TEST1')['BSURF_I'][0]

    assert abs(bsurfi - bsurfiatnf) < sf_scale(bsurfiatnf)


def test_distance(query_derived, query_atnf):
    """
    Test the derived distance.
    """

    # distance derived from parallax for TEST1
    dist = query_derived.get_pulsar('TEST1')['DIST'][0]
    distatnf = query_atnf.get_pulsar('TEST1')['DIST'][0]

    assert abs(dist - distatnf) < sf_scale(distatnf)

    dist1 = query_derived.get_pulsar('TEST1')['DIST1'][0]

    assert dist == dist1

    # distance obtained from DM and DM1
    dist = query_derived.get_pulsar('TEST2')['DIST'][0]
    distatnf = query_atnf.get_pulsar('TEST2')['DIST'][0]

    assert abs(dist - distatnf) < sf_scale(distatnf)

    distdm = query_derived.get_pulsar('TEST2')['DIST_DM'][0]

    assert dist == distdm

    dist = query_derived.get_pulsar('TEST2')['DIST1'][0]
    distatnf = query_atnf.get_pulsar('TEST2')['DIST1'][0]

    assert abs(dist - distatnf) < sf_scale(distatnf)

    distdm = query_derived.get_pulsar('TEST2')['DIST_DM1'][0]

    assert dist == distdm

    # distance obtained from DIST_AMN and DIST_AMX
    dist = query_derived.get_pulsar('TEST3')['DIST'][0]
    distatnf = query_atnf.get_pulsar('TEST3')['DIST'][0]

    assert abs(dist - distatnf) < sf_scale(distatnf)

    dist = query_derived.get_pulsar('TEST3')['DIST1'][0]
    distatnf = query_atnf.get_pulsar('TEST3')['DIST1'][0]

    assert abs(dist - distatnf) < sf_scale(distatnf)


def test_luminosity(query_derived, query_atnf):
    """
    Test derived spin-down luminosity parameters.
    """

    edot = query_derived.get_pulsar('TEST1')['EDOT'][0]
    edotatnf = query_atnf.get_pulsar('TEST1')['EDOT'][0]

    assert abs(edot - edotatnf) < sf_scale(edotatnf)

    edotd2 = query_derived.get_pulsar('TEST1')['EDOTD2'][0]
    edotd2atnf = query_atnf.get_pulsar('TEST1')['EDOTD2'][0]

    assert abs(edotd2 - edotd2atnf) < sf_scale(edotd2atnf)


def test_derived_radio_luminosity(query_derived, query_atnf):
    """
    Test the derived radio luminosity.
    """

    rlum = query_derived.get_pulsar('TEST1')['R_LUM'][0]
    rlumatnf = query_atnf.get_pulsar('TEST1')['R_LUM'][0]

    assert abs(rlum - rlumatnf) < sf_scale(rlumatnf)

    rlum14 = query_derived.get_pulsar('TEST1')['R_LUM14'][0]
    rlum14atnf = query_atnf.get_pulsar('TEST1')['R_LUM14'][0]

    assert abs(rlum14 - rlum14atnf) < sf_scale(rlum14atnf)


def test_derived_proper_motion(query_derived, query_atnf):
    """
    Test the derived proper motion values.
    """

    pmtot = query_derived.get_pulsar('TEST1')['PMTOT'][0]
    pmtotatnf = query_atnf.get_pulsar('TEST1')['PMTOT'][0]
    pmtoterr = query_derived.get_pulsar('TEST1')['PMTOT_ERR'][0]
    pmtoterratnf = query_atnf.get_pulsar('TEST1')['PMTOT_ERR'][0]

    assert abs(pmtot - pmtotatnf) < sf_scale(pmtotatnf)
    assert round_err(pmtoterr, pmtoterratnf)

    vtrans = query_derived.get_pulsar('TEST1')['VTRANS'][0]
    vtransatnf = query_atnf.get_pulsar('TEST1')['VTRANS'][0]

    assert abs(vtrans - vtransatnf) < sf_scale(vtransatnf)

    # proper motion in galactic coordinates
    # NOTE: these are currently left out of the test as psrcat performs an
    # additional velocity correction to their local standard of rest by
    # removing both the solar system velocity and their local galactic
    # rotation
    # pml = query_derived.get_pulsar('TEST1')['PML'][0]
    # pmlatnf = query_atnf.get_pulsar('TEST1')['PML'][0]

    # assert abs(pml - pmlatnf) < sf_scale(pmlatnf)

    # pmb = query_derived.get_pulsar('TEST1')['PMB'][0]
    # pmbatnf = query_atnf.get_pulsar('TEST1')['PMB'][0]

    # assert abs(pmb - pmbatnf) < sf_scale(pmbatnf)


def test_derived_pb_pbdot(query_derived, query_atnf):
    """
    Test binary period and period derivative from orbital frequency.
    """

    pb = query_derived.get_pulsar('TEST4')['PB'][0]
    pbatnf = query_atnf.get_pulsar('TEST4')['PB'][0]
    pberr = query_derived.get_pulsar('TEST4')['PB_ERR'][0]
    pberratnf = query_atnf.get_pulsar('TEST4')['PB_ERR'][0]

    assert abs(pb - pbatnf) < sf_scale(pbatnf)
    assert round_err(pberr, pberratnf)

    pbdot = query_derived.get_pulsar('TEST4')['PBDOT'][0]
    pbdotatnf = query_atnf.get_pulsar('TEST4')['PBDOT'][0]
    pbdoterr = query_derived.get_pulsar('TEST4')['PBDOT_ERR'][0]
    pbdoterratnf = query_atnf.get_pulsar('TEST4')['PBDOT_ERR'][0]

    assert abs(pbdot - pbdotatnf) < sf_scale(pbdotatnf)
    assert round_err(pbdoterr, pbdoterratnf)


def test_pdot_to_fdot(query):
    """
    Test the conversion functions from Pdot to fdot and vice versa.
    """

    from psrqpy.utils import pdot_to_fdot, fdot_to_pdot

    vela = query.get_pulsar('J0835-4510')
    crab = query.get_pulsar('J0534+2200')

    F1v = vela["F1"].data[0]
    P1v = vela["P1"].data[0]

    F0v = vela["F0"].data[0]
    P0v = vela["P0"].data[0]

    F1c = crab["F1"].data[0]
    P1c = crab["P1"].data[0]

    F0c = crab["F0"].data[0]
    P0c = crab["P0"].data[0]

    # test pdot to fdot
    with pytest.raises(ValueError):
        pdot_to_fdot(P1v)

    # test individual value
    f1calcv1 = pdot_to_fdot(P1v, period=P0v)
    f1calcv2 = pdot_to_fdot(P1v, frequency=F0v)

    assert isinstance(f1calcv1, float)
    assert isinstance(f1calcv2, float)
    assert f1calcv1 == f1calcv2
    assert np.allclose([f1calcv1], [F1v])

    # test pair of values
    f1calc1 = pdot_to_fdot([P1v, P1c], period=[P0v, P0c])
    f1calc2 = pdot_to_fdot([P1v, P1c], frequency=[F0v, F0c])

    assert len(f1calc1) == 2
    assert len(f1calc2) == 2
    assert np.allclose(f1calc1, f1calc2)

    # test fdot to pdot
    with pytest.raises(ValueError):
        fdot_to_pdot(F1v)

    # test individual value
    p1calcv1 = fdot_to_pdot(F1v, period=P0v)
    p1calcv2 = fdot_to_pdot(F1v, frequency=F0v)

    assert isinstance(p1calcv1, float)
    assert isinstance(p1calcv2, float)
    assert p1calcv1 == p1calcv2
    assert np.allclose([p1calcv1], [P1v])

    # test pair of values
    p1calc1 = fdot_to_pdot([F1v, F1c], period=[P0v, P0c])
    p1calc2 = fdot_to_pdot([F1v, F1c], frequency=[F0v, F0c])

    assert len(p1calc1) == 2
    assert len(p1calc2) == 2
    assert np.allclose(p1calc1, p1calc2)


def test_derived_gw_parameters(query):
    """
    Test derived gravitational-wave parameters.
    """

    from psrqpy.utils import h0_to_ellipticity, ellipticity_to_q22, gw_luminosity

    crab = query.get_pulsar('J0534+2200')

    # check derived h0 spin-down upper limit
    expectedh0 = 1.4e-24  # e.g., Table 3 of https://arxiv.org/abs/2007.14251
    derivedh0 = crab["H0_SD"].data[0]

    # make sure value is within 5% of expected value
    assert 0.95 < derivedh0 / expectedh0 < 1.05

    # check derived luminosity
    expectedL = 4.5e31  # e.g., Table 2 of https://arxiv.org/abs/2007.14251
    derivedL = crab["EDOT"].to("W").data[0]

    # make sure value is within 5% of expected value
    assert 0.95 < derivedL / expectedL < 1.05

    # check ellipticity conversion
    expectedell = 1.8e-4 * 4.1  # see Sec 3 of https://arxiv.org/abs/0805.4758
    ellipticity = h0_to_ellipticity(crab["H0_SD"], crab["F0"], crab["DIST"])

    # make sure value is within 5% of expected value
    assert 0.95 < ellipticity / expectedell < 1.05

    # check ellipticity to q22
    expectedQ22 = 12.6e32 / 0.021  # see Table 3 of https://arxiv.org/abs/2007.14251
    q22 = ellipticity_to_q22(ellipticity)

    # make sure value is within 5% of expected value
    assert 0.95 < q22 / expectedQ22 < 1.05

    # check GW luminosity
    h0ul = 1.5e-26  # see Table 3 of https://arxiv.org/abs/2007.14251
    gwlum = gw_luminosity(h0ul, crab["F0"], crab["DIST"])
    lumrat = gwlum / crab["EDOT"].to("W").data[0]
    expectedrat = 0.01 ** 2

    # make sure value is within 10% of expected value
    assert 0.9 < lumrat / expectedrat < 1.1

    # use query table functions
    ell = query.gw_ellipticity()
    ellipticity = ell[ell["PSRJ"] == "J0534+2200"]["ELL"]
    assert 0.95 < ellipticity / expectedell < 1.05

    # use query table functions
    q22t = query.gw_mass_quadrupole()
    q22 = q22t[q22t["PSRJ"] == "J0534+2200"]["Q22"]
    assert 0.95 < q22 / expectedQ22 < 1.05


# TEST EXCEPTIONS #
def test_bad_database():
    """
    Try loading in random file.
    """

    baddbfile = 'sdhfjjdf'  # bad database file
    with pytest.raises(RuntimeError):
        _ = QueryATNF(loadfromdb=baddbfile)


@pytest.mark.disable_socket
def test_download_db():
    """
    Try downloading the database with the socket disabled.
    """

    with pytest.raises(RuntimeError):
        _ = QueryATNF(checkupdate=True)


@pytest.mark.disable_socket
def test_download_glitch_table():
    """
    Try downloading the glitch table with the socket disabled.
    """

    from psrqpy.utils import get_glitch_catalogue

    with pytest.raises(RuntimeError):
        _ = get_glitch_catalogue()


@pytest.mark.disable_socket
def test_download_gc_table():
    """
    Try downloading the globular cluster table with the socket disabled.
    """

    from psrqpy.utils import get_gc_catalogue

    with pytest.raises(RuntimeError):
        _ = get_gc_catalogue()


@pytest.mark.disable_socket
def test_download_msp_table():
    """
    Try downloading the MSP table with the socket disabled.
    """

    from psrqpy.utils import get_msp_catalogue

    with pytest.raises(RuntimeError):
        _ = get_msp_catalogue()


def test_sort_exception(query):
    """
    Test exception in sort method.
    """

    curkey = query.sort_key
    sortval = 'kgsdkfkfd'  # random sort parameter
    with pytest.raises(KeyError):
        query.sort(sort_attr=sortval)

    sort_key = 1.2  # non-string sort key
    with pytest.raises(ValueError):
        query.sort_key = sort_key

    # reset sort key
    query.sort_key = curkey
