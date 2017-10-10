import unittest
from ..model import dy_dt
import numpy as np
from scipy.integrate import odeint


class TestModel(unittest.TestCase):
    def setUp(self):
        print("Setup testing")
        
    def test_length(self):
        y0 = np.ones((10, ), dtype = np.float64)
        args = (1., 1., 1., 1., 1., 1., 0.5) # these 7 float values represent the inputs IL2 through k11rev
                        
        self.assertEqual(len(dy_dt(y0, 0, *args)),10)
                
