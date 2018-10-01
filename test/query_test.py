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
