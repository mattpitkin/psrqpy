from pathlib import Path

import pytest
from psrqpy import QueryATNF

TESTS_DIR = Path(__file__).parent


@pytest.fixture(scope="module")
def query():
    return QueryATNF()


@pytest.fixture(scope="module")
def query_atnf():
    return QueryATNF(loadfromdb=TESTS_DIR / 'derived_catalogue.db')


@pytest.fixture(scope="session")
def query_derived():
    return QueryATNF(loadfromdb=TESTS_DIR / 'test_catalogue.db')
