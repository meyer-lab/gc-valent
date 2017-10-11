import unittest
from ..model import dy_dt
import numpy as np
from scipy.integrate import odeint


class TestModel(unittest.TestCase):
    def setUp(self):
        print("Setup testing")
        self.ts = np.array([0.0, 100000.0])
        self.y0 = np.random.lognormal(0., 1., 10)
        self.args1 = list(np.random.lognormal(0., 1., 7))
        self.args = tuple(self.args1)
        # need to convert args from an array to a tuple of numbers
        
    def test_equilibrium(self):
        y, fullout = odeint(dy_dt, self.y0, self.ts, self.args,
                    full_output = True, mxstep = 5000)
        z = np.linalg.norm(dy_dt(y[1, :], 0, *self.args)) # have the sum of squares for all the dy_dt values be z
        self.assertTrue(z < 1E-5)