import pytest
from psrqpy import QueryATNF


@pytest.fixture(scope="module")
def query():
    return QueryATNF()


@pytest.fixture(scope="module")
def query_atnf():
    return QueryATNF(loadfromdb='test/derived_catalogue.db')


@pytest.fixture(scope="session")
def query_derived():
    return QueryATNF(loadfromdb='test/test_catalogue.db')
