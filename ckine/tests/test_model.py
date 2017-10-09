import unittest
from ..model import dy_dt
import numpy as np
from scipy.integrate import odeint


class TestModel(unittest.TestCase):
    def setUp(self):
        print("Setup testing")

    def test_equilibrium(self):
        ts = np.array([0.0, 100000.0])
        y0 = np.ones((18, ), dtype=np.float64)
        args = (1., 1., 1., 1., 1., 1., 0.5, 1., 1., 1., 1., 1., 1., 1., 1.)  # these 8 float values represent the inputs IL2 through k11rev
        y = odeint(dy_dt, y0, ts, args, mxstep=5000)

        print(y)
        z = np.linalg.norm(dy_dt(y[1, :], 0, *args)) # have the sum of squares for all the dy_dt values be z
        self.assertLess(z, 1E-5)

