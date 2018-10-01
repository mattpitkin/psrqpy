"""
Test script.
"""

import pytest
from psrqpy import QueryATNF
import numpy as np


@pytest.mark.enable_socket
def test_crab(query):
    """
    Test that the Crab pulsar is present and the frequency is as expected, i.e.
    the frequency rounds down to 29 Hz (should be OK for another ~80 years!)
    """

    f0 = query.get_pulsar('J0534+2200')['F0'][0]

    assert np.floor(f0) == 29.0

@pytest.mark.enable_socket
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

@pytest.mark.enable_socket
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

### Test derived parameters ###
@pytest.mark.enable_socket
def test_derived_p0(query_derived, query_atnf):
    """
    Test the derived period value against the values from the ATNF Pulsar
    Catalogue.
    """

    p0 = query_derived.get_pulsar('TEST1')['P0'][0]
    p0atnf = query_atnf.get_pulsar('TEST1')['P0'][0]

    p0err = query_derived.get_pulsar('TEST1')['P0_ERR'][0]
    p0atnferr = query_atnf.get_pulsar('TEST1')['P0_ERR'][0]

    assert abs(p0-p0atnf) < p0atnferr

    # get ATNF derived error value
    errexp = np.floor(np.log10(p0atnferr))    # exponent
    errval = np.round(p0atnferr/10**errexp)   # value before exponent

    # ATNF derived errors are always rounded up
    derval = np.ceil(p0err/10**errexp)

    assert errval == derval


### Test exceptions ###
def test_bad_database():
    """
    Try loading in random file.
    """

    baddbfile = 'sdhfjjdf'  # bad database file
    with pytest.raises(IOError):
        query = QueryATNF(loadfromdb=baddbfile)

def test_download_db():
    """
    Try downloading the database without the socket being enabled.
    """

    with pytest.raises(RuntimeError):
        query = QueryATNF(checkupdate=True)

@pytest.mark.enable_socket
def test_sort_exception(query):
    """
    Test exception in sort method.
    """

    sortval = 'kgsdkfkfd'  # random sort parameter
    with pytest.raises(KeyError):
        query.sort(sort_attr=sortval)
