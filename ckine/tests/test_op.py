import unittest
import theano
import theano.tensor as T
from theano.tests import unittest_tools as utt
import numpy as np
from ..differencing_op import runCkineOp, runCkineKineticOp, runCkineDoseOp
from ..model import nSpecies, nParams, getTotalActiveSpecies


class TestOp(unittest.TestCase):
    def setUp(self):
        self.unkV = [np.full(nParams(), 0.9, dtype=np.float64)]
        self.cond = np.ones(nSpecies(), dtype=np.float64)
        self.conditions = np.ones((10, 6), dtype=np.float64)
        self.tol = {'abs_tol': 0.001, 'rel_tol': 0.001}

    def test_runCkineOp_T0(self):
        utt.verify_grad(runCkineOp(np.array([0.0])), self.unkV)

    def test_runCkineOp(self):
        utt.verify_grad(runCkineOp(np.array([100000.])), self.unkV)

    def test_runCkineKineticOp(self):
        ts = np.linspace(0, 1000, dtype=np.float64)

        utt.verify_grad(runCkineKineticOp(ts, self.cond), self.unkV, **self.tol)

    #def test_runCkineDoseOp(self):
    #    utt.verify_grad(runCkineDoseOp(np.array(1.0), self.cond, self.conditions), self.unkV, **self.tol)

    def test_runCkineDoseOp_noActivity(self):
        """ Test that in the absence of ligand most values and gradients are zero. """
        a = T.dvector()
        fexpr = runCkineDoseOp(np.array(10.0), getTotalActiveSpecies().astype(np.float64), np.zeros_like(self.conditions))(a)
        J = theano.scan(lambda i, y, x : T.grad(fexpr[i], a), sequences=T.arange(fexpr.shape[0]), non_sequences=[fexpr, a])[0]

        f = theano.function([a], fexpr)
        fprime = theano.function([a], J)

        self.assertAlmostEqual(np.max(f(self.unkV[0])), 0.0)
        print(fprime(self.unkV[0]))
