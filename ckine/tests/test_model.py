import unittest
from ..model import dy_dt
import numpy as np
from scipy.integrate import odeint


class TestModel(unittest.TestCase):
    def setUp(self):
        print("Setup testing")
        
    def test_conservation(self):
        ts = np.array([0.0, 100000.0])
        y0 = np.ones((10, ), dtype = np.float64)
        args = (1., 1., 1., 1., 1., 1., 0.5) # these 7 float values represent the inputs IL2 through k11rev
        y, fullout = odeint(dy_dt, y0, ts, args,
                    full_output = True, mxstep = 5000)
        