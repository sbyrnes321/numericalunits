# -*- coding: utf-8 -*-
"""
very basic tests on numericalunits
"""
import unittest
import numericalunits as nu
from math import isclose

class TestStuff(unittest.TestCase):

    def setUp(self):
        pass

    def test_everything(self):
        """just some very basic smoke tests"""
        # example from README
        self.assertTrue(isclose(5 * nu.mL, 5e21 * nu.nm**3, rel_tol=1e-9))

        # example from README
        Efield = 1e5 * (nu.V / nu.cm)
        force = nu.e * Efield
        accel = force / nu.me
        self.assertTrue(isclose(accel, 1.75882002e18 * nu.m / nu.s**2, rel_tol=1e-6))

        # check nu_eval()
        self.assertTrue(isclose(nu.nu_eval('kg'), nu.kg, rel_tol=1e-9))
        self.assertTrue(isclose(nu.nu_eval('kg * m / s**2'), nu.kg * nu.m / nu.s**2, rel_tol=1e-9))
        self.assertTrue(isclose(nu.nu_eval('kg**-3.6'), nu.kg**-3.6, rel_tol=1e-9))

        # make sure reset_units('SI') works
        nu.reset_units('SI')
        self.assertTrue(isclose(nu.G, 1e-4, rel_tol=1e-9))


if __name__ == '__main__':
    unittest.main()
