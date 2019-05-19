# -*- coding: utf-8 -*-
"""
very basic tests on numericalunits
"""
import unittest
import numericalunits as nu

class TestStuff(unittest.TestCase):

    def setUp(self):
        pass

    def assert_almost_equal(self, a, b, rtol):
        """helper function to check if two floats are approximately equal,
        allowing for rounding errors etc. Similar to math.isclose in py3."""
        self.assertTrue(abs(a-b) <= rtol * (abs(a) + abs(b)))

    def test_everything(self):
        """just some very basic smoke tests"""
		# example from README
        x = 5 * nu.mL
        self.assert_almost_equal(x, 5e21 * nu.nm**3, rtol=1e-9)

        # example from README
        Efield = 1e5 * (nu.V / nu.cm)
        force = nu.e * Efield
        accel = force / nu.me
        self.assert_almost_equal(accel, 1.75882002e18 * nu.m / nu.s**2, rtol=1e-6)

        # make sure reset_units('SI') works
        nu.reset_units('SI')
        self.assert_almost_equal(nu.G, 1e-4, rtol=1e-9)

if __name__ == '__main__':
    unittest.main()
