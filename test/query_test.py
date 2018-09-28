"""
Test script.
"""

from __future__ import absolute_import

import unittest
from psrqpy import QueryATNF
from psrqpy.config import *


class TestDerived(unittest.TestCase):
    def setUp(self):
        self.query_derived = QueryATNF(loadfromdb='test_catalogue.db')
        self.query_atnf = QueryATNF(loadfromdb='derived_catalogue.db')

    def tearDown(self):
        del self.query_known
        del self.query_derived

    def test_derived_p0(self):
        """
        Test the derived period value against the values from the ATNF Pulsar
        Catalogue.
        """

        p0 = self.query_derived.get_pulsar('TEST1')['P0'].value
        p0atnf = self.query_atnf.get_pulsar('TEST1')['P0'].value

        tol = 1e-9  # tolerance of how close values should be

        self.assertTrue(abs(p0-p0atnf) < tol)

        #p0err =

if __name__ == '__main__':
    unittest.main()
