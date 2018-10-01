import pytest
from psrqpy import QueryATNF
from pytest_socket import disable_socket

def pytest_runtest_setup():
    """
    Disable sockets by default, so that we can test correct errors when
    trying to access the database.
    """

    disable_socket()

@pytest.fixture(scope="module")
def query():
    return QueryATNF()

@pytest.fixture(scope="module")
def query_atnf():
    return QueryATNF(loadfromdb='test/derived_catalogue.db')

@pytest.fixture(scope="session")
def query_derived():
    return QueryATNF(loadfromdb='test/test_catalogue.db')