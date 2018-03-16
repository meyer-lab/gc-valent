"""
Unit test file.
"""
import unittest
import numpy as np
from hypothesis import given, settings
from hypothesis.strategies import floats
from hypothesis.extra.numpy import arrays as harrays
from ..model import dy_dt, fullModel, solveAutocrine, getTotalActiveCytokine, solveAutocrineComplete, runCkine

settings.register_profile("ci", max_examples=1000)
settings.load_profile("ci")

class TestModel(unittest.TestCase):
    def assertPosEquilibrium(self, X, func):
        """Assert that all species came to equilibrium."""
        # All the species abundances should be above zero
        self.assertGreater(np.min(X), -1.0E-7)

        # Test that it came to equilirbium
        self.assertLess(np.linalg.norm(func(X)) / (1.0 + np.sum(X)), 1E-5)

    def assertConservation(self, y, y0, IDX):
        """Assert the conservation of species throughout the experiment."""
        species_delta = y - y0

        # Check for conservation of species sum
        self.assertAlmostEqual(np.sum(species_delta[IDX]), 0.0, places=5)

    def setUp(self):
        self.ts = np.array([0.0, 100000.0])
        self.y0 = np.random.lognormal(0., 1., 26)
        self.args = np.random.lognormal(0., 1., 17)
        self.tfargs = np.random.lognormal(0., 1., 11)
        # need to convert args from an array to a tuple of numbers

        if (self.tfargs[2] > 1.):
            self.tfargs[2] = self.tfargs[2] - np.floor(self.tfargs[2])

    def test_length(self):
        self.assertEqual(len(dy_dt(self.y0, 0, self.args)), self.y0.size)

    @given(y0=harrays(np.float, 26, elements=floats(0, 10)))
    def test_conservation(self, y0):
        """Check for the conservation of each of the initial receptors."""
        dy = dy_dt(y0, 0.0, self.args)
        #Check for conservation of gc
        self.assertConservation(dy, 0.0, np.array([2, 5, 7, 8, 9, 13, 15, 16, 17, 20, 24, 21, 25]))
        #Check for conservation of IL2Rb
        self.assertConservation(dy, 0.0, np.array([1, 4, 6, 8, 9, 12, 14, 16, 17]))
        #Check for conservation of IL2Ra
        self.assertConservation(dy, 0.0, np.array([0, 3, 6, 7, 9]))
        #Check for conservation of IL15Ra
        self.assertConservation(dy, 0.0, np.array([10, 11, 14, 15, 17]))
        #Check for conservation of IL7Ra
        self.assertConservation(dy, 0.0, np.array([18, 19, 21]))
        #Check for Conservation of IL9R
        self.assertConservation(dy, 0.0, np.array([22, 23, 25]))

    @given(y0=harrays(np.float, 2*26 + 4, elements=floats(0, 10)))
    def test_conservation_full(self, y0):
        """In the absence of trafficking, mass balance should hold in both compartments."""
        kw = np.zeros(11, dtype=np.float64)

        dy = fullModel(y0, 0.0, self.args, kw)

        #Check for conservation of gc
        self.assertConservation(dy, 0.0, np.array([2, 5, 7, 8, 9, 13, 15, 16, 17, 20, 24, 21, 25]))
        #Check for conservation of IL2Rb
        self.assertConservation(dy, 0.0, np.array([1, 4, 6, 8, 9, 12, 14, 16, 17]))
        #Check for conservation of IL2Ra
        self.assertConservation(dy, 0.0, np.array([0, 3, 6, 7, 9]))
        #Check for conservation of IL15Ra
        self.assertConservation(dy, 0.0, np.array([10, 11, 14, 15, 17]))
        #Check for conservation of IL7Ra
        self.assertConservation(dy, 0.0, np.array([18, 19, 21]))
        #Check for Conservation of IL9R
        self.assertConservation(dy, 0.0, np.array([22, 23, 25]))

        #Check for conservation of gc
        self.assertConservation(dy, 0.0, np.array([2, 5, 7, 8, 9, 13, 15, 16, 17, 20, 24, 21, 25]) + 26)
        #Check for conservation of IL2Rb
        self.assertConservation(dy, 0.0, np.array([1, 4, 6, 8, 9, 12, 14, 16, 17]) + 26)
        #Check for conservation of IL2Ra
        self.assertConservation(dy, 0.0, np.array([0, 3, 6, 7, 9]) + 26)
        #Check for conservation of IL15Ra
        self.assertConservation(dy, 0.0, np.array([10, 11, 14, 15, 17]) + 26)
        #Check for conservation of IL7Ra
        self.assertConservation(dy, 0.0, np.array([18, 19, 21]) + 26)
        #Check for Conservation of IL9R
        self.assertConservation(dy, 0.0, np.array([22, 23, 25]) + 26)

    def test_fullModel(self):
        """Assert the two functions solveAutocrine and solveAutocrine complete return the same values."""
        yOut = solveAutocrine(self.tfargs)

        yOut2 = solveAutocrineComplete(self.args, self.tfargs)

        kw = self.args.copy()

        kw[0:4] = 0.

        # Autocrine condition assumes no cytokine present, and so no activity
        self.assertAlmostEqual(getTotalActiveCytokine(0, yOut), 0.0, places=5)

        # Autocrine condition assumes no cytokine present, and so no activity
        self.assertAlmostEqual(getTotalActiveCytokine(0, yOut2), 0.0, places=5)

        self.assertPosEquilibrium(yOut, lambda y: fullModel(y, 0.0, kw, self.tfargs))

    @given(y0=harrays(np.float, 2*26 + 4, elements=floats(0, 10)))
    def test_reproducible(self, y0):

        dy1 = fullModel(y0, 0.0, self.args, self.tfargs)

        dy2 = fullModel(y0, 1.0, self.args, self.tfargs)

        dy3 = fullModel(y0, 2.0, self.args, self.tfargs)

        # Test that there's no difference
        self.assertLess(np.linalg.norm(dy1 - dy2), 1E-8)

        # Test that there's no difference
        self.assertLess(np.linalg.norm(dy1 - dy3), 1E-8)

    @given(vec=harrays(np.float, 27, elements=floats(0.001, 10.0)), sortF=floats(0.1, 0.9))
    def test_runCkine(self, vec, sortF):
        vec = np.insert(vec, 2, sortF)
        # 11 trafRates and 17 rxnRates
        trafRates = vec[0:11]
        rxnRates = vec[11:28]

        ys, retVal = runCkine(self.ts, rxnRates, trafRates)
        
        # test that return value of runCkine isn't negative (model run didn't fail)
        self.assertGreaterEqual(retVal, 0)
        
    def test_failed_runCkine(self):
        '''Internal CVode error in CVode
        At t = 0 and h = 0, the corrector convergence test failed repeatedly or with |h| = hmin.
        In module: CVODE
        Error code: -4
        CVode error in CVode. Code: -4
        Model run failed
        Endocytosis: 3.146053504525751e-304
        activeEndo: 0.10614172503188912
        sortF: 4.597500267118165e-237
        kRec: 1.168921624346571e-212
        kDeg: 6.175221406423857e-25
        Receptor expression: [0.e+00 0.e+00 1.e-07 0.e+00 0.e+00 0.e+00]
        .....Reaction rates.....
        IL2: 500.0
        IL15: 0.0
        IL7: 0.0
        IL9: 0.0
        kfwd: 3.1356349315053714e-172
        k5rev: 0.037367428484416224
        k6rev: 0.07733074217829684
        [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]        '''
        
        Endocytosis = 3.146053504525751e-304
        activeEndo = 0.10614172503188912
        sortF = 4.597500267118165e-237
        kRec = 1.168921624346571e-212
        kDeg = 6.175221406423857e-25
        kfwd = 3.1356349315053714e-172
        k5rev = 0.037367428484416224
        k6rev = 0.07733074217829684
        tfrRates = np.array([Endocytosis, activeEndo, sortF, kRec, kDeg, 0., 0., 1e-07, 0., 0., 0.])
        rxnRates = np.array([500., 0., 0., 0., kfwd, k5rev, k6rev, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        
        ys, retVal = runCkine(self.ts, rxnRates, trfRates)
        
        # test that return value of runCkine isn't negative (model run didn't fail)
        self.assertGreaterEqual(retVal, 0)

        '''Internal CVode error in CVode
        At t = 82.2539 and h = 1.78983e-09, the error test failed repeatedly or with |h| = hmin.
        In module: CVODE
        Error code: -3
        CVode error in CVode. Code: -3
        Model run failed
        Endocytosis: 0.0002682816512461031
        activeEndo: 0.08750230849382867
        sortF: 0.022692859699932316
        kRec: 0.01237493760188316
        kDeg: 5.048972752964616e-08
        Receptor expression: [0.00000000e+00 1.00000000e-07 8.71966236e-15 0.00000000e+00
         0.00000000e+00 0.00000000e+00]
        .....Reaction rates.....
        IL2: 500.0
        IL15: 0.0
        IL7: 0.0
        IL9: 0.0
        kfwd: 76.08505747988166
        k5rev: 0.13781451657116525
        k6rev: 0.14745955118122073
        [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]        '''