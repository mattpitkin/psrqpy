"""
Test script.
"""

import pytest
from psrqpy import QueryATNF
import numpy as np
from pandas import Series
import pytest_socket


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


def test_update(query):
    """
    Test update method.
    """

    # no name set (and column not a Series)
    column = 1
    with pytest.raises(ValueError):
        query.update(column)

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

    fig = query.ppdot(showtypes='BINARY', showGCs=True)


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
    Test the derived galactic longitude and latitude and galactic
    cartesian coordinates.
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

    # xx = query_derived.get_pulsar('TEST1')['XX'][0]
    # xxatnf = query_atnf.get_pulsar('TEST1')['XX'][0]

    # yy = query_derived.get_pulsar('TEST1')['YY'][0]
    # yyatnf = query_atnf.get_pulsar('TEST1')['YY'][0]

    # zz = query_derived.get_pulsar('TEST1')['ZZ'][0]
    # zzatnf = query_atnf.get_pulsar('TEST1')['ZZ'][0]

    # assert abs(xx-xxatnf) < sf_scale(xxatnf)
    # assert abs(yy-yyatnf) < sf_scale(yyatnf)
    # assert abs(zz-zzatnf) < sf_scale(zzatnf)


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


# TEST EXCEPTIONS #
def test_bad_database():
    """
    Try loading in random file.
    """

    baddbfile = 'sdhfjjdf'  # bad database file
    with pytest.raises(IOError):
        query = QueryATNF(loadfromdb=baddbfile)


@pytest.mark.disable_socket
def test_download_db():
    """
    Try downloading the database with the socket disabled.
    """

    with pytest.raises(RuntimeError):
        query = QueryATNF(checkupdate=True)


def test_sort_exception(query):
    """
    Test exception in sort method.
    """

    curkey = query.sort_key
    sortval = 'kgsdkfkfd'  # random sort parameter
    with pytest.raises(KeyError):
        query.sort(sort_attr=sortval)

    # reset sort key
    query.sort_key = curkey
