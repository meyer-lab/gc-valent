from .model import solveAutocrine, fullModel, __active_species_IDX
from scipy.integrate import odeint
import numpy as np

_epsilon = np.sqrt(np.finfo(float).eps)

def approx_jacobian(x,func,epsilon,*args):
    """Approximate the Jacobian matrix of callable function func

       * Parameters
         x       - The state vector at which the Jacobian matrix is desired
         func    - A vector-valued function of the form f(x,*args)
         epsilon - The peturbation used to determine the partial derivatives
         *args   - Additional arguments passed to func

       * Returns
         An array of dimensions (lenf, lenx) where lenf is the length
         of the outputs of func, and lenx is the number of

       * Notes
         The approximation is done using forward differences

    """
    x0 = np.asfarray(x)
    f0 = func(*((x0,)+args))
    jac = np.zeros([len(x0),len(f0)])
    dx = np.zeros(len(x0))
    for i in range(len(x0)):
       dx[i] = epsilon
       jac[i] = (func(*((x0+dx,)+args)) - f0)/epsilon
       dx[i] = 0.0
    return jac.transpose()


trafRates = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])
t = np.array([0.0, 100000.0])
rxnRates = np.array([1.,1.,1.,1.,1.])
y = np.ones(56)
approx_jacobian( solveAutocrine(trafRates) , fullModel, _epsilon, y, t, rxnRates, trafRates, __active_species_IDX)