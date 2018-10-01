"""
Test script.
"""

from __future__ import absolute_import

import unittest
from psrqpy import QueryATNF
from psrqpy.config import *
import numpy as np

class TestCrab(unittest.TestCase):
    """
    Test that the Crab pulsar is present and the frequency is as expected, i.e.
    the frequency rounds down to 29 Hz (should be OK for another ~80 years!)
    """

    def setUp(self):
        self.query_crab = QueryATNF(psrs='J0534+2200', params='F0')

    def tearDown(self):
        del self.query_crab

    def test_f0(self):
        f0 = self.query_crab.table['F0'][0]

        self.assertTrue(np.floor(f0) == 29.0)
    

class TestDerived(unittest.TestCase):
    def setUp(self):
        self.query_derived = QueryATNF(loadfromdb='test/test_catalogue.db')
        self.query_atnf = QueryATNF(loadfromdb='test/derived_catalogue.db')

    def tearDown(self):
        del self.query_atnf
        del self.query_derived

    def test_derived_p0(self):
        """
        Test the derived period value against the values from the ATNF Pulsar
        Catalogue.
        """

        p0 = self.query_derived.get_pulsar('TEST1')['P0'][0]
        p0atnf = self.query_atnf.get_pulsar('TEST1')['P0'][0]

        p0err = self.query_derived.get_pulsar('TEST1')['P0_ERR'][0]
        p0atnferr = self.query_atnf.get_pulsar('TEST1')['P0_ERR'][0]

        self.assertTrue(abs(p0-p0atnf) < p0atnferr)

        # get ATNF derived error value
        errexp = np.floor(np.log10(p0atnferr))    # exponent
        errval = np.round(p0atnferr/10**errexp)   # value before exponent

        # ATNF derived errors are always rounded up
        derval = np.ceil(p0err/10**errexp)

        self.assertTrue(errval == derval)


if __name__ == '__main__':
    unittest.main()
